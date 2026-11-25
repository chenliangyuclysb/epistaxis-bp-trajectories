"""Run the P1 cohort build after the environment and hosp setup."""

from importlib import import_module
from pathlib import Path
import random
import sys

import duckdb

import config as cfg


ROOT = Path(__file__).resolve().parents[1]
DERIVED = ROOT / "derived"
OUTPUTS = ROOT / "outputs"
TABLES = OUTPUTS / "tables"
SITE_FIELDS = (
    ("mcmed", "CSN", "MRN", "ED_dispo", "NULL::VARCHAR"),
    ("mimic", "stay_id", "subject_id", "disposition", "i.hadm_id"),
)


def combine_cohorts(con):
    for site, visit_id, patient_id, disposition, hadm_id in SITE_FIELDS:
        if site == "mimic":
            inpatient_fields = "h.* EXCLUDE (visit_id)"
            inpatient_join = "JOIN p1_mimic_inpatient_outcome h ON h.visit_id = i.stay_id"
        else:
            inpatient_fields = ("NULL::BOOLEAN AS inpatient_ent_outcome_eligible, "
                                "NULL::BOOLEAN AS inpatient_ent_hemostasis_event")
            inpatient_join = ""

        con.execute(f"""
            CREATE OR REPLACE TEMP TABLE p1_{site}_cohort AS
            SELECT '{site}' AS site, i.{visit_id} AS visit_id, i.{patient_id} AS patient_id,
                {hadm_id} AS hadm_id, i.arrival_ts, i.departure_ts, i.age_years,
                i.{disposition} AS disposition, i.dx_match, i.cc_match,
                i.prior_epistaxis_365d, i.is_first_eligible_visit, i.incomplete_visit,
                i.is_c1, i.is_c2, i.bp_readings, i.bp_span_minutes,
                c.* EXCLUDE (visit_id), r.* EXCLUDE (visit_id), {inpatient_fields}
            FROM p1_{site}_index i
            JOIN p1_{site}_covariates c ON c.visit_id = i.{visit_id}
            JOIN p1_{site}_returns r ON r.visit_id = i.{visit_id}
            {inpatient_join}
        """)


def validate_cohorts(con):
    for site, visit_id, _, _, _ in SITE_FIELDS:
        n, unique, invalid = con.execute(f"""
            SELECT count(*), count(DISTINCT visit_id), count(*) FILTER (WHERE
                visit_id IS NULL OR patient_id IS NULL OR arrival_ts IS NULL
                OR departure_ts <= arrival_ts OR age_years < $adult
                OR is_c1 IS DISTINCT FROM coalesce(
                    bp_readings >= $min_readings AND bp_span_minutes >= $span, false
                )
                OR (is_c2 AND (NOT is_c1 OR incomplete_visit))
                OR is_c3 IS DISTINCT FROM (is_c1 AND NOT prior_htn)
                OR (return_7d_event AND return_30d_event IS DISTINCT FROM true)
                OR return_outcome_eligible IS DISTINCT FROM (NOT incomplete_visit)
                OR (NOT return_outcome_eligible AND (
                    return_7d_event IS NOT NULL OR return_30d_event IS NOT NULL
                ))
                OR (return_30d_event AND (return_minutes_within_30d <= 0
                    OR return_minutes_within_30d > $max_return_hours * 60))
            )
            FROM p1_{site}_cohort
        """, {"adult": cfg.MIN_ADULT_AGE, "min_readings": cfg.C1_MIN_BP_READINGS,
              "span": cfg.C1_MIN_SPAN_MINUTES, "max_return_hours": cfg.RETURN_30D_HOURS}).fetchone()
        initial = con.execute(f"SELECT count({visit_id}) FROM p1_{site}_index").fetchone()[0]
        retained_readings, recorded_readings, orphans = con.execute(f"""
            SELECT (SELECT count(*) FROM p1_bp_{site}),
                (SELECT coalesce(sum(bp_readings), 0) FROM p1_{site}_cohort),
                count(*) FROM p1_bp_{site} b
                LEFT JOIN p1_{site}_cohort c USING (visit_id)
                WHERE c.visit_id IS NULL OR b.minutes_from_arrival < 0
                    OR b.minutes_from_arrival > ?
        """, [cfg.ED_MAX_HOURS * 60]).fetchone()
        if n != unique or n != initial or invalid or retained_readings != recorded_readings or orphans:
            raise ValueError("P1 visit, cohort, outcome or BP referential checks failed.")


def write_parquets(con, literal):
    for site, _, _, _, _ in SITE_FIELDS:
        queries = (
            (f"cohort_{site}", f"""
                SELECT * FROM p1_{site}_cohort ORDER BY patient_id, arrival_ts, visit_id
            """),
            (f"bp_long_{site}", f"""
                SELECT visit_id, source, minutes_from_arrival, sbp, dbp FROM p1_bp_{site}
                ORDER BY visit_id, minutes_from_arrival, source, sbp, dbp
            """),
        )
        for name, query in queries:
            destination = DERIVED / f"{name}.parquet"
            pending = destination.with_suffix(".parquet.tmp")
            try:
                con.execute(f"COPY ({query}) TO {literal(pending)} (FORMAT PARQUET, COMPRESSION ZSTD)")
                pending.replace(destination)
            finally:
                pending.unlink(missing_ok=True)


def write_procedure_code_audit(con, write_csv):
    rows = con.execute("""
        SELECT icd_version, icd_code, long_title, count(*) AS dictionary_records
        FROM p1_mimic_procedure_labels
        GROUP BY icd_version, icd_code, long_title
        ORDER BY icd_version, icd_code
    """).fetchall()
    write_csv(TABLES / "p1_inpatient_ent_candidate_codes.csv",
              ("icd_version", "icd_code", "procedure_title", "dictionary_records"), rows)


def main():
    p0 = import_module("01_p0_audit")
    setup = import_module("02_p1_setup")
    random.seed(cfg.RANDOM_SEED)
    step = "P1 preparation check"
    baseline = OUTPUTS / "raw_fingerprint_p1_before.csv"
    try:
        setup.verify_project_runtime()
        if not baseline.is_file() or not (TABLES / "p1_hosp_checksum.csv").is_file():
            raise RuntimeError("Run the P1 setup and verification first.")
        with duckdb.connect(str(DERIVED / "epistaxis.duckdb"), read_only=True) as con:
            con.execute(f"SET temp_directory = {p0.literal(DERIVED / 'duckdb_tmp')}")
            con.execute("SET preserve_insertion_order = false")
            setup.register_hosp(con)
            for step, module in (
                ("1.2 finalized label mappings", "02_p1_dictionaries"),
                ("P1 cohort selection and BP", "02_p1_cohort_core"),
                ("P1 covariates", "02_p1_covariates"),
                ("P1 return and inpatient procedure outcomes", "02_p1_outcomes"),
            ):
                print(step, flush=True)
                import_module(module).run(con)
            step = "P1 cohort assembly and checks"
            combine_cohorts(con)
            validate_cohorts(con)
            step = "P1 aggregate tables and local derived Parquets"
            import_module("02_p1_tables").run(con, p0.write_csv)
            write_procedure_code_audit(con, p0.write_csv)
            write_parquets(con, p0.literal)
        step = "P1 raw fingerprint comparison"
        p0.fingerprint("p1_after")
        if baseline.read_bytes() != (OUTPUTS / "raw_fingerprint_p1_after.csv").read_bytes():
            raise RuntimeError("Original source fingerprints have changed.")
    except Exception as error:
        # CSV parser failures can disclose source values when formatted.
        print(f"P1 stopped at {step} ({type(error).__name__}).", file=sys.stderr)
        return 1
    print("P1 derived relations and aggregate tables passed checks.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
