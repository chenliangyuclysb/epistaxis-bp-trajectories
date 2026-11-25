"""Apply G1 covariate and landmark rules without changing the historical P1 reports."""

from importlib import import_module
from pathlib import Path
import sys

import duckdb


ROOT = Path(__file__).resolve().parents[1]
DERIVED = ROOT / "derived"
TABLES = ROOT / "outputs" / "tables"
DICTIONARY = ROOT / "outputs" / "dictionary_final"
SITES = ("mcmed", "mimic")
NEW_COLUMNS = (
    "prior_htn_hosp_diagnosis", "prior_htn_ed_diagnosis", "escalation_reversal_confirmed",
    "landmark_eligible_any_admission_diagnostic", "hemostasis_level3_by_60min_any_admission_diagnostic",
    "hemostasis_level3_after_60min_any_admission_diagnostic",
)
ESCALATION_FIELDS = (
    "hemostasis_level3", "hemostasis_level3_by_60min", "hemostasis_level3_after_60min", "landmark_eligible",
    "hemostasis_level3_any_admission_diagnostic", "first_escalation_minutes",
)


def prepare_indices(con, literal):
    for site in SITES:
        con.execute(f"""
            CREATE OR REPLACE TEMP TABLE p2_input_{site} AS
            SELECT * FROM read_parquet({literal(DERIVED / f'cohort_{site}.parquet')})
        """)
    con.execute("""
        CREATE OR REPLACE TEMP TABLE p1_mcmed_index AS
        SELECT p.visit_id AS CSN, p.patient_id AS MRN, p.arrival_ts, p.departure_ts, p.dx_match,
               v.Dx_ICD10, v.Dx_ICD9, v.ED_dispo, v.Dispo_time
        FROM p2_input_mcmed p JOIN mcmed_visits v ON p.visit_id = v.CSN AND p.patient_id = v.MRN
    """)
    con.execute("""
        CREATE OR REPLACE TEMP TABLE p1_mimic_index AS
        SELECT visit_id AS stay_id, patient_id AS subject_id, arrival_ts FROM p2_input_mimic
    """)


def stage_reversal_dictionary(con, dictionary):
    dictionary.map_labels(
        con, "mcmed_orders", ("Order_type", "Procedure_name"), "p1_dict_mcmed_orders",
        dictionary.ED_FIELDS + dictionary.ESCALATION_FIELDS,
        lambda kind, name: dictionary.classify_aht_ed(name, kind) | dictionary.classify_escalation(kind, name), None,
    )
    source_fields = ("Order_type", "Procedure_name", "reversal_family", "record_count")
    result = con.execute(f"""
        SELECT {', '.join(source_fields)} FROM p1_dict_mcmed_orders WHERE reversal_candidate
        ORDER BY Order_type, Procedure_name
    """)
    records = [dict(zip(source_fields, (kind, label, family, count)))
               | dictionary.classify_g1_reversal(kind, label, family)
               for kind, label, family, count in result.fetchall()]
    con.execute("""
        CREATE OR REPLACE TEMP TABLE p2_dict_mcmed_reversal (
            Order_type VARCHAR, Procedure_name VARCHAR, reversal_confirmed BOOLEAN
        )
    """)
    if records:
        con.executemany(
            "INSERT INTO p2_dict_mcmed_reversal VALUES (?, ?, ?)",
            [(row["Order_type"], row["Procedure_name"], row["reversal_confirmed"]) for row in records],
        )
    return (*source_fields, "reversal_confirmed", "g1_hemostasis_level", "g1_reversal_exclusion_reason"), records


def stage_updates(con, covariates):
    covariates.create_mcmed_g1_history(con)
    covariates.create_mimic_g1_history(con)
    covariates.create_mcmed_order_events(con, g1=True)
    covariates.create_mcmed_escalation(con)
    con.execute("CREATE OR REPLACE TEMP TABLE p2_mcmed_escalation_before AS SELECT * FROM p1_mcmed_hemostasis")
    covariates.create_mcmed_escalation(con, g1=True)

    for site in SITES:
        existing = {col[0] for col in con.execute(f"DESCRIBE p2_input_{site}").fetchall()}
        existing_additions = [field for field in NEW_COLUMNS if field in existing]
        base = "old.*"
        if existing_additions:
            base += " EXCLUDE (" + ", ".join(existing_additions) + ")"
        replacements = [f"h.{field} AS {field}" for field in ("ckd", "liver_disease", "hht")]
        additions = [f"NULL::BOOLEAN AS {field}" for field in NEW_COLUMNS]
        additional_join = ""
        if site == "mcmed":
            replacements.extend(f"s.{field} AS {field}" for field in ESCALATION_FIELDS)
            additions[2:] = [f"s.{field} AS {field}" for field in NEW_COLUMNS[2:]]
            additional_join = "JOIN p1_mcmed_hemostasis s USING (visit_id)"
        else:
            additions[:2] = [f"h.{field} AS {field}" for field in NEW_COLUMNS[:2]]
        con.execute(f"""
            CREATE OR REPLACE TEMP TABLE p2_updated_{site} AS
            SELECT {base} REPLACE ({', '.join(replacements)}), {', '.join(additions)}
            FROM p2_input_{site} old
            JOIN p2_{site}_g1_history h USING (visit_id)
            {additional_join}
        """)


def validate_updates(con):
    previous = "hemostasis_level3_by_60min"
    current = "hemostasis_level3_after_60min"
    eligible = "landmark_eligible"
    suffix = "_any_admission_diagnostic"
    for site in SITES:
        invalid = con.execute(f"""
            SELECT count(*) FROM p2_input_{site} a FULL JOIN p2_updated_{site} u USING (visit_id)
            WHERE a.visit_id IS NULL OR u.visit_id IS NULL
               OR a.prior_htn IS DISTINCT FROM u.prior_htn
               OR a.has_medication_list IS DISTINCT FROM u.has_medication_list
               OR a.is_c3 IS DISTINCT FROM u.is_c3
        """).fetchone()[0]
        output_counts = con.execute(f"SELECT count(*), count(DISTINCT visit_id) FROM p2_updated_{site}").fetchone()
        if invalid or output_counts[0] != output_counts[1]:
            raise RuntimeError("G1 update changed the index keys or unchanged HTN variables.")
    broken_components = con.execute("""
        SELECT count(*) FROM p2_updated_mimic
        WHERE prior_htn_diagnosis IS DISTINCT FROM (prior_htn_hosp_diagnosis OR prior_htn_ed_diagnosis)
           OR prior_htn IS DISTINCT FROM (prior_htn_hosp_diagnosis OR prior_htn_ed_diagnosis OR prior_htn_medication)
    """).fetchone()[0]
    bad_risksets = con.execute(f"""
        SELECT count(*) FROM p2_updated_mcmed
        WHERE ({eligible} AND {previous}) OR ({current} IS NOT NULL) IS DISTINCT FROM {eligible}
           OR ({eligible}{suffix} AND {previous}{suffix}) OR ({current}{suffix} IS NOT NULL) IS DISTINCT FROM {eligible}{suffix}
    """).fetchone()[0]
    if broken_components or bad_risksets:
        raise RuntimeError("G1 HTN-source or event-free landmark checks failed.")


def audit_counts(con, table, field):
    return con.execute(f"SELECT count(*) FILTER (WHERE {field}), count({field}) FROM {table}").fetchone()


def aggregate_audit(con):
    rows = []
    for site, label in (("mcmed", "MC-MED"), ("mimic", "MIMIC-IV-ED")):
        index_count = con.execute(f"SELECT count(*) FROM p2_updated_{site}").fetchone()[0]
        comparisons = (
            ("ckd", f"p2_{site}_g1_history", "ckd_before"),
            ("liver_disease", f"p2_{site}_g1_history", "liver_before"),
            ("hht", f"p2_{site}_g1_history", "hht_before"),
            ("prior_htn", f"p2_input_{site}", "prior_htn"),
            ("has_medication_list", f"p2_input_{site}", "has_medication_list"),
        )
        if site == "mcmed":
            comparisons += tuple((field, "p2_mcmed_escalation_before", field) for field in ESCALATION_FIELDS
                                 if field != "first_escalation_minutes")
            comparisons += tuple((field, None, None) for field in NEW_COLUMNS[2:])
        else:
            comparisons += tuple((field, None, None) for field in NEW_COLUMNS[:2])
        for field, before_table, before_field in comparisons:
            before = audit_counts(con, before_table, before_field) if before_table else (None, None)
            after = audit_counts(con, f"p2_updated_{site}", field)
            rows.append((label, field, index_count, *before, *after))
    return rows


def write_updates(con, literal):
    destinations = [DERIVED / f"cohort_{site}.parquet" for site in SITES]
    pending = [path.with_suffix(".parquet.g1.tmp") for path in destinations]
    try:
        for site, path in zip(SITES, pending):
            con.execute(f"""
                COPY (SELECT * FROM p2_updated_{site} ORDER BY patient_id, arrival_ts, visit_id)
                TO {literal(path)} (FORMAT PARQUET, COMPRESSION ZSTD)
            """)
        for staged, destination in zip(pending, destinations):
            staged.replace(destination)
    finally:
        for path in pending:
            path.unlink(missing_ok=True)


def main():
    p0 = import_module("01_p0_audit")
    setup = import_module("02_p1_setup")
    dictionaries = import_module("02_p1_dictionaries")
    covariates = import_module("02_p1_covariates")
    step = "project runtime and P1 inputs"
    try:
        setup.verify_project_runtime()
        with duckdb.connect(str(DERIVED / "epistaxis.duckdb"), read_only=True) as con:
            con.execute(f"SET temp_directory = {p0.literal(DERIVED / 'duckdb_tmp')}")
            con.execute("SET preserve_insertion_order = false")
            setup.register_hosp(con)
            prepare_indices(con, p0.literal)
            step = "G1 reversal order mappings"
            fields, records = stage_reversal_dictionary(con, dictionaries)
            step = "G1 revised diagnosis and landmark rules"
            stage_updates(con, covariates)
            validate_updates(con)
            counts = aggregate_audit(con)
            step = "updated cohort Parquets and aggregate audits"
            write_updates(con, p0.literal)
        dictionaries.write_table(DICTIONARY / "mcmed_reversal_g1_mapping.csv", fields, records)
        p0.write_csv(TABLES / "p2_g1_updates.csv",
                     ("site", "field", "index_visits", "before_positive", "before_observed",
                      "after_positive", "after_observed"), counts)
    except Exception as error:
        # Source parser exceptions can expose individual values.
        print(f"G1 updates stopped at {step} ({type(error).__name__}).", file=sys.stderr)
        return 1
    print("G1 updates completed and verified.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
