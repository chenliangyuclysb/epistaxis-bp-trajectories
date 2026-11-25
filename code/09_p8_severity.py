"""Run the G7 post hoc severity checks without changing frozen P4 artifacts."""

import argparse
import csv
from importlib import import_module
from pathlib import Path
import re
import sys

import duckdb

import config as cfg


TABLES = Path("outputs/tables")
SITES = (("mcmed", "MC-MED"), ("mimic", "MIMIC-IV-ED"))
PHENOTYPES = ("N", "T", "P")
OUTCOME = "return_7d_event"
AUDIT_HEADER = ("name", "decision", "reason", "start_to_discharge_records", "start_to_discharge_visits", "all_pyxis_records")
DESCRIPTIVE_HEADER = (
    "site", "phenotype", "visits", "patients", "local_hemostatic_visits", "local_hemostatic_percent", "local_wilson_ci95_low",
    "local_wilson_ci95_high", "hemoglobin_observed_visits", "hemoglobin_available_percent", "hemoglobin_median_g_dl",
    "hemoglobin_q1_g_dl", "hemoglobin_q3_g_dl",
)
MODEL_HEADER = (
    "analysis", "site", "model", "stratum", "contrast", "eligible_visits", "modeled_visits", "modeled_patients", "modeled_events",
    "n_visits", "n_events", "t_visits", "t_events", "p_visits", "p_events", "rr", "ci95_low", "ci95_high", "parameters",
    "events_per_parameter", "stability", "fit_status", "notes",
)


def stage_index(con):
    fields = (
        "visit_id", "patient_id", "arrival_ts", "departure_ts", "age_years", "sex", "any_antithrombotic_use", "prior_htn",
        "triage_acuity_group", "hemostasis_level2", "hemoglobin",
    )
    parts = [f"""
        SELECT '{site}' AS site, {', '.join('c.' + field for field in fields)}, p.bp_phenotype
        FROM read_parquet('derived/cohort_{short}.parquet') c
        LEFT JOIN read_parquet('derived/phenotypes_{short}.parquet') p USING (visit_id)
        WHERE c.is_c2
    """ for short, site in SITES]
    con.execute(f"CREATE OR REPLACE TEMP TABLE p8_return_index AS {' UNION ALL '.join(parts)}")
    total, distinct, invalid = con.execute("""
        SELECT count(*), count(DISTINCT (site, visit_id)), count(*) FILTER (
            WHERE visit_id IS NULL OR patient_id IS NULL OR arrival_ts IS NULL OR departure_ts IS NULL
               OR arrival_ts >= departure_ts OR NOT coalesce(bp_phenotype IN ('N','T','P'), false)
               OR (site = 'MC-MED' AND hemostasis_level2 IS NULL))
        FROM p8_return_index
    """).fetchone()
    if not total or total != distinct or invalid:
        raise ValueError("The frozen return-cohort index is incomplete or duplicated.")


def pyxis_decision(label):
    if re.search(cfg.P8_PYXIS_TXA_PATTERN, label):
        return "include", "Tranexamic acid is included across routes per G7"
    if re.search(cfg.P8_PYXIS_SYSTEMIC_THROMBIN_PATTERN, label):
        return "exclude", "Systemic antithrombin or prothrombin product"
    if re.search(cfg.P8_PYXIS_NONLOCAL_FORM_PATTERN, label):
        return "exclude", "Explicit systemic or ophthalmic formulation"
    if re.search(cfg.P8_PYXIS_PHENYLEPHRINE_PATTERN, label):
        if re.search(cfg.P8_PYXIS_LOCAL_FORM_PATTERN, label):
            return "include", "Phenylephrine explicitly labeled nasal or topical"
        return "exclude", "Phenylephrine label has no explicit nasal or topical formulation"
    return "include", "Candidate local hemostatic product or nasal vasoconstrictor"


def stage_pyxis_audit(con, write_csv):
    con.execute("""
        CREATE OR REPLACE TEMP TABLE p8_pyxis_candidates AS
        SELECT stay_id, name, try_cast(charttime AS TIMESTAMP) AS dispense_time
        FROM mimic_pyxis WHERE regexp_matches(name, ?)
    """, [cfg.P8_PYXIS_CANDIDATE_PATTERN])
    labels = con.execute("""
        SELECT p.name, count(c.visit_id), count(DISTINCT c.visit_id), count(*)
        FROM p8_pyxis_candidates p LEFT JOIN p8_return_index c
          ON c.site = 'MIMIC-IV-ED' AND c.visit_id = p.stay_id
             AND p.dispense_time BETWEEN c.arrival_ts AND c.departure_ts
        GROUP BY p.name ORDER BY p.name
    """).fetchall()
    if not labels:
        raise ValueError("No G7 Pyxis candidate labels were found.")
    records = [(name, *pyxis_decision(name), records, visits, all_records) for name, records, visits, all_records in labels]
    write_csv("p8_pyxis_local_hemostatic_audit.csv", AUDIT_HEADER, records)
    con.execute("CREATE OR REPLACE TEMP TABLE p8_pyxis_decisions (name VARCHAR, decision VARCHAR)")
    con.executemany("INSERT INTO p8_pyxis_decisions VALUES (?, ?)", [(row[0], row[1]) for row in records])
    con.execute("""
        CREATE OR REPLACE TEMP TABLE p8_mimic_local_hemostasis AS
        SELECT c.visit_id, count(*) AS local_dispensing_records
        FROM p8_pyxis_candidates p JOIN p8_pyxis_decisions d USING (name)
        JOIN p8_return_index c ON c.site = 'MIMIC-IV-ED' AND c.visit_id = p.stay_id
          AND p.dispense_time BETWEEN c.arrival_ts AND c.departure_ts
        WHERE d.decision = 'include' GROUP BY c.visit_id
    """)
    return len(records), sum(row[1] == "include" for row in records)


def stage_hemoglobin(con, write_csv):
    con.execute("""
        CREATE OR REPLACE TEMP TABLE p8_hemoglobin_dictionary AS
        SELECT itemid, label, fluid, category
        FROM mimic_hosp_d_labitems WHERE lower(trim(fluid)) = ? AND lower(trim(label)) = ?
    """, [cfg.P8_HEMOGLOBIN_FLUID, cfg.P8_HEMOGLOBIN_LABEL])
    n, distinct = con.execute("SELECT count(*), count(DISTINCT itemid) FROM p8_hemoglobin_dictionary").fetchone()
    if not n or n != distinct:
        raise ValueError("Blood hemoglobin dictionary codes are missing or ambiguous.")
    con.execute("""
        CREATE OR REPLACE TEMP TABLE p8_mimic_hb_candidates AS
        WITH selected AS (
            SELECT c.visit_id, c.departure_ts, l.itemid,
                   try_cast(l.charttime AS TIMESTAMP) AS chart_ts, try_cast(l.storetime AS TIMESTAMP) AS store_ts,
                   l.labevent_id, try_cast(l.valuenum AS DOUBLE) AS numeric_value,
                   upper(replace(trim(l.valueuom), ' ', '')) AS unit
            FROM mimic_hosp_labevents l JOIN p8_hemoglobin_dictionary d USING (itemid)
            JOIN p8_return_index c ON c.site = 'MIMIC-IV-ED' AND l.subject_id = c.patient_id
              AND try_cast(l.charttime AS TIMESTAMP) BETWEEN c.arrival_ts AND c.departure_ts
        )
        SELECT *, CASE WHEN isfinite(numeric_value) AND list_contains($gdl, unit) THEN numeric_value
                       WHEN isfinite(numeric_value) AND list_contains($gl, unit) THEN numeric_value / $unit_conversion
                  END AS hemoglobin_g_dl
        FROM selected
    """, {"gdl": list(cfg.P8_HEMOGLOBIN_UNITS_G_DL), "gl": list(cfg.P8_HEMOGLOBIN_UNITS_G_L),
          "unit_conversion": cfg.P8_HEMOGLOBIN_GL_PER_GDL})
    con.execute("""
        CREATE OR REPLACE TEMP TABLE p8_mimic_first_hb AS
        SELECT visit_id, departure_ts, itemid, chart_ts, store_ts, labevent_id, hemoglobin_g_dl
        FROM p8_mimic_hb_candidates WHERE hemoglobin_g_dl IS NOT NULL
        QUALIFY row_number() OVER (
            PARTITION BY visit_id ORDER BY chart_ts, store_ts NULLS LAST, try_cast(labevent_id AS BIGINT), labevent_id
        ) = 1
    """)
    rows = con.execute("""
        WITH candidates AS (
            SELECT itemid, count(*) AS records, count(DISTINCT visit_id) AS visits,
                   count(*) FILTER (WHERE hemoglobin_g_dl IS NOT NULL) AS eligible_results,
                   count(*) FILTER (WHERE NOT coalesce(isfinite(numeric_value), false)) AS invalid_numeric_results,
                   count(*) FILTER (WHERE isfinite(numeric_value) AND hemoglobin_g_dl IS NULL) AS unsupported_unit_results
            FROM p8_mimic_hb_candidates GROUP BY itemid
        ), chosen AS (
            SELECT itemid, count(*) AS selected_visits,
                   count(*) FILTER (WHERE store_ts > departure_ts) AS stored_after_ed,
                   count(*) FILTER (WHERE store_ts IS NULL) AS store_time_unavailable
            FROM p8_mimic_first_hb GROUP BY itemid
        )
        SELECT d.itemid, d.label, d.fluid, d.category, coalesce(c.records, 0), coalesce(c.visits, 0),
               coalesce(c.eligible_results, 0), coalesce(c.invalid_numeric_results, 0), coalesce(c.unsupported_unit_results, 0),
               coalesce(s.selected_visits, 0), coalesce(s.stored_after_ed, 0), coalesce(s.store_time_unavailable, 0)
        FROM p8_hemoglobin_dictionary d LEFT JOIN candidates c USING (itemid) LEFT JOIN chosen s USING (itemid)
        ORDER BY try_cast(d.itemid AS INTEGER), d.itemid
    """).fetchall()
    write_csv("p8_hemoglobin_audit.csv", (
        "item_code", "label", "fluid", "category", "ed_records", "ed_visits", "valid_ed_results", "invalid_numeric_results",
        "unsupported_unit_results", "selected_first_results", "selected_first_results_stored_after_ed", "selected_results_store_time_missing",
    ), rows)


def analysis_frame(con):
    sources = [f"""
        SELECT '{site}' AS site, visit_id, {OUTCOME}
        FROM read_parquet('derived/cohort_{short}.parquet') WHERE is_c2
    """ for short, site in SITES]
    frame = con.execute(f"""
        WITH frozen_outcomes AS ({' UNION ALL '.join(sources)})
        SELECT i.site, i.visit_id, i.patient_id, i.bp_phenotype, i.age_years, i.sex, i.any_antithrombotic_use,
               i.prior_htn, i.triage_acuity_group, o.{OUTCOME},
               CASE WHEN i.site = 'MC-MED' THEN i.hemostasis_level2
                    ELSE x.local_dispensing_records IS NOT NULL END AS local_hemostatic_treatment,
               CASE WHEN i.site = 'MC-MED' THEN i.hemoglobin ELSE h.hemoglobin_g_dl END AS hemoglobin_g_dl
        FROM p8_return_index i JOIN frozen_outcomes o USING (site, visit_id)
        LEFT JOIN p8_mimic_local_hemostasis x ON x.visit_id = i.visit_id AND i.site = 'MIMIC-IV-ED'
        LEFT JOIN p8_mimic_first_hb h ON h.visit_id = i.visit_id AND i.site = 'MIMIC-IV-ED'
        ORDER BY i.site, i.patient_id, i.visit_id
    """).df()
    if frame[["local_hemostatic_treatment", OUTCOME]].isna().any().any() or frame.duplicated(["site", "visit_id"]).any():
        raise ValueError("Post hoc proxy or outcome flags are incomplete or duplicated.")
    with (TABLES / "p4_main_group_counts.csv").open(encoding="utf-8-sig", newline="") as stream:
        original = {(row["site"], row["phenotype"]): (int(row["observed_outcomes"]), int(row["observed_events"]))
                    for row in csv.DictReader(stream) if row["endpoint"] == "7-day return" and row["site"] != "Pooled"}
    if original != {(site, label): (len(group), int(group[OUTCOME].sum()))
                    for (site, label), group in frame.groupby(["site", "bp_phenotype"])}:
        raise ValueError("Return-cohort counts differ from frozen P4 counts.")
    frame["hemoglobin_per_g_dl"] = frame.hemoglobin_g_dl / cfg.P8_HEMOGLOBIN_MODEL_INCREMENT_G_DL
    return frame


def descriptive(frame, wilson):
    rows, coverage = [], []
    for _, site in SITES:
        stratum = frame.loc[frame.site.eq(site)]
        observed = int(stratum.hemoglobin_g_dl.notna().sum())
        percentage = 100 * observed / len(stratum)
        gate = percentage >= cfg.P8_MIMIC_HEMOGLOBIN_MIN_AVAILABLE_PERCENT if site == SITES[1][1] else None
        coverage.append((site, len(stratum), observed, percentage, gate))
        for label in PHENOTYPES:
            group = stratum.loc[stratum.bp_phenotype.eq(label)]
            n = len(group)
            local = int(group.local_hemostatic_treatment.sum())
            hb = group.hemoglobin_g_dl.dropna()
            if len(hb) > 1:
                q1, median, q3 = hb.quantile([0.25, 0.5, 0.75])
            else:
                q1 = median = q3 = None
            rows.append((site, label, n, group.patient_id.nunique(), local, 100 * local / n, *wilson(local, n),
                         len(hb), 100 * len(hb) / n, median, q1, q3))
    return rows, coverage


def candidate_group_counts(group):
    return {label: (int(group.bp_phenotype.eq(label).sum()), int(group.loc[group.bp_phenotype.eq(label), OUTCOME].sum()))
            for label in PHENOTYPES}


def result_rows(models, group, analysis, site, title, stratum, adjust, additional_covariates,
                contrasts=("P", "N"), skip=None):
    if skip:
        fit = None
        counts = candidate_group_counts(group)
        base = (len(group), None, None, None)
        note = f"{skip}. Displayed phenotype counts are for the eligible visits."
    else:
        fit = models.fit_model(group, OUTCOME, adjustment=adjust, continuous=additional_covariates)
        counts = fit.group_counts
        base = (fit.n_eligible, fit.n, fit.patients, fit.events)
        note = "Post hoc. Local treatment during the visit may follow phenotype assessment."
        if fit.cross_check is not None:
            note += f" Triggered clustered Poisson check after GEE; Poisson status {fit.cross_check.status}."
        if fit.n_gt_one:
            note += " Some fitted GEE means exceeded one."
    n_terms = len(fit.columns) if fit is not None and fit.params is not None else None
    epv = fit.events / n_terms if n_terms else None
    stability = ("unstable" if epv < cfg.P8_MIN_EVENTS_PER_PARAMETER else "stable") if epv is not None else "not assessed"
    group_totals = [value for label in PHENOTYPES for value in counts.get(label, (0, 0))]
    for label in contrasts:
        if fit is None:
            effect = models.blank_effect("Not fitted")
        else:
            effect = fit.estimate(f"bp_phenotype={label}")
        yield (analysis, site, title, stratum, f"{label} vs T", *base, *group_totals,
               *(effect[key] for key in ("estimate", "ci95_low", "ci95_high")), n_terms, epv, stability, effect["status"], note)


def fit_all(frame, models, coverage):
    pooled = "Pooled"
    local = ("local_hemostatic_treatment",)
    both = (*local, "hemoglobin_per_g_dl")
    labels = ("P4 model with local hemostatic treatment", "P4 model with local hemostatic treatment and first ED hemoglobin")
    rows = list(result_rows(models, frame, "B pooled", pooled, labels[0], "All", "pooled", local))
    for site, adjustment in ((SITES[0][1], "mcmed"), (SITES[1][1], "mimic")):
        group = frame.loc[frame.site.eq(site)]
        rows.extend(result_rows(models, group, "C site-specific", site, labels[0], "All", adjustment, local))
        if adjustment == "mcmed":
            rows.extend(result_rows(models, group, "C site-specific with hemoglobin", site, labels[1], "All", adjustment, both))
    has_mcmed_hb = coverage[0][2] > 0
    accepts_mimic_hb = coverage[1][-1]
    no_pool = None if has_mcmed_hb and accepts_mimic_hb else "Hemoglobin did not meet the cross-center availability rule"
    rows.extend(result_rows(models, frame, "D pooled with hemoglobin", pooled, labels[1], "All", "pooled", both, skip=no_pool))
    for positive in (False, True):
        group = frame.loc[frame.local_hemostatic_treatment.eq(positive)]
        reason = "Stratum event count was below the minimum for estimation" if int(group[OUTCOME].sum()) < cfg.P8_MIN_STRATUM_EVENTS else None
        rows.extend(result_rows(models, group, "E stratified", pooled, "P4 model within the local hemostatic-treatment stratum",
                                "Local treatment" if positive else "No local treatment", "pooled", (), contrasts=("P",), skip=reason))
    return rows


def run(con, audit_only=False):
    setup = import_module("02_p1_setup")
    models = import_module("05_p4_models")
    stage_index(con)
    labels, included = stage_pyxis_audit(con, models.write_csv)
    print(f"P8 Pyxis audit written before model fitting: {labels} candidate labels, {included} included.", flush=True)
    if audit_only:
        return
    checksum_path = TABLES / "p8_hosp_labs_checksum.csv"
    with checksum_path.open(encoding="utf-8-sig", newline="") as stream:
        verified = {row["member"] for row in csv.DictReader(stream)
                    if row["status"] == "pass" and row["expected_sha256"] == row["actual_sha256"]}
    if verified != set(setup.P8_LAB_PAYLOAD):
        raise RuntimeError("Verify the opt-in P8 lab files before the severity analysis.")
    setup.register_hosp(con, include_p8_labs=True)
    stage_hemoglobin(con, models.write_csv)
    frame = analysis_frame(con)
    wilson = import_module("02_p1_tables").wilson
    descriptors, coverage = descriptive(frame, wilson)
    models.write_csv("p8_posthoc_severity_descriptive.csv", DESCRIPTIVE_HEADER, descriptors)
    models.write_csv("p8_hemoglobin_coverage.csv", ("site", "return_visits", "hemoglobin_available_visits",
                     "hemoglobin_available_percent", "satisfies_mimic_threshold"), coverage)
    models.write_csv("p8_posthoc_severity.csv", MODEL_HEADER, fit_all(frame, models, coverage))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit-pyxis", action="store_true", help="Write only the pre-model Pyxis candidate-label audit.")
    args = parser.parse_args()
    try:
        import_module("02_p1_setup").verify_project_runtime()
        with duckdb.connect("derived/epistaxis.duckdb", read_only=True) as con:
            con.execute("SET temp_directory = 'derived/duckdb_tmp'")
            con.execute("SET preserve_insertion_order = false")
            run(con, audit_only=args.audit_pyxis)
    except Exception as error:
        # Parser and modeling exception strings can contain protected input values.
        print(f"P8 8.1 stopped ({type(error).__name__}).", file=sys.stderr)
        return 1
    print("P8 8.1 aggregate outputs completed.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
