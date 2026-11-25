"""Build the patient-level Aim 2 follow-up data and aggregate outputs."""

from importlib import import_module
from pathlib import Path
import sys

import duckdb
import numpy as np


DERIVED = Path("derived")
TABLES = Path("outputs/tables")
SITES = (("mcmed", "MC-MED"), ("mimic", "MIMIC-IV-ED"))
VARIANTS = {
    "mcmed": (
        ("primary", "index", "is_c3"),
        ("same_day_htn", "index_same_day_htn", "NOT prior_htn_same_day_sensitivity"),
    ),
    "mimic": (
        ("primary", "index", "is_c3"),
        ("no_prior_ed", "index_no_prior_ed", "NOT (prior_htn_medication OR prior_htn_hosp_diagnosis)"),
    ),
}
CONTINUOUS_FIELDS = (
    "age_years", "sbp_first", "dbp_first", "triage_heart_rate", "triage_pain_score",
)
CATEGORICAL_FIELDS = (
    "bp_phenotype", "sex", "triage_acuity_group", "race", "ethnicity", "arrival_mode",
)
BINARY_FIELDS = (
    "has_medication_list", "anticoagulant_use", "antiplatelet_use",
    "any_antithrombotic_use", "ckd", "liver_disease", "prior_epistaxis_365d",
)
BASELINE_FIELDS = CONTINUOUS_FIELDS + CATEGORICAL_FIELDS + BINARY_FIELDS


def create_indices(con):
    rows = []
    for short, site in SITES:
        con.execute(f"""
            CREATE OR REPLACE TEMP TABLE p3_{short}_c1 AS
            SELECT c.visit_id, c.patient_id, c.arrival_ts, c.departure_ts, c.incomplete_visit,
                c.is_c3, c.prior_htn, c.prior_htn_medication, c.prior_htn_hosp_diagnosis,
                c.prior_htn_ed_diagnosis, c.prior_htn_same_day_sensitivity,
                c.age_years, c.sex, c.triage_acuity_group, c.triage_heart_rate, c.triage_pain_score,
                c.race, c.ethnicity, c.arrival_mode, c.has_medication_list,
                c.anticoagulant_use, c.antiplatelet_use, c.any_antithrombotic_use,
                c.ckd, c.liver_disease, c.prior_epistaxis_365d,
                p.bp_phenotype, p.sbp_first, p.dbp_first
            FROM read_parquet('derived/cohort_{short}.parquet') c
            JOIN read_parquet('derived/phenotypes_{short}.parquet') p USING (visit_id)
            WHERE c.is_c1
        """)
        c1_count, c1_patients, phenotype_count, invalid = con.execute(f"""
            SELECT count(*), count(DISTINCT patient_id),
                (SELECT count(*) FROM read_parquet('derived/phenotypes_{short}.parquet')),
                count(*) FILTER (WHERE is_c3 IS DISTINCT FROM (NOT prior_htn)
                    OR incomplete_visit IS NULL OR patient_id IS NULL OR visit_id IS NULL
                    OR bp_phenotype NOT IN ('N', 'T', 'P') OR bp_phenotype IS NULL
                    OR NOT coalesce(isfinite(age_years), false)
                    OR NOT coalesce(isfinite(sbp_first), false)
                    OR arrival_ts IS NULL OR departure_ts IS NULL OR departure_ts <= arrival_ts)
            FROM p3_{short}_c1
        """).fetchone()
        if c1_count != phenotype_count or invalid:
            raise ValueError("C1 or pretreatment phenotype alignment failed.")
        for name, suffix, predicate in VARIANTS[short]:
            view = f"p3_{short}_{suffix}"
            counts = con.execute(f"""
                SELECT count(*), count(DISTINCT patient_id),
                    count(*) FILTER (WHERE incomplete_visit),
                    count(*) FILTER (WHERE NOT incomplete_visit),
                    count(DISTINCT patient_id) FILTER (WHERE NOT incomplete_visit)
                FROM p3_{short}_c1 WHERE {predicate}
            """).fetchone()
            # Choose the index without querying follow-up, after S3.2 eligibility.
            con.execute(f"""
                CREATE OR REPLACE TEMP TABLE {view} AS
                SELECT *
                FROM p3_{short}_c1
                WHERE ({predicate}) AND NOT incomplete_visit
                QUALIFY row_number() OVER (PARTITION BY patient_id ORDER BY arrival_ts, visit_id) = 1
            """)
            selected, patients = con.execute(
                f"SELECT count(*), count(DISTINCT patient_id) FROM {view}"
            ).fetchone()
            if selected != patients or patients != counts[-1]:
                raise ValueError("The earliest eligible index is not one visit per patient.")
            rows.append((site, name, c1_count, c1_patients, *counts, selected))
    return rows


def combine_followup(con):
    expressions = []
    for short, site in SITES:
        for variant, suffix, _ in VARIANTS[short]:
            expressions.append(f"""
                SELECT '{site}' AS site, '{variant}' AS index_variant, i.visit_id, i.patient_id,
                    {", ".join("i." + field for field in BASELINE_FIELDS)},
                    f.has_followup, f.htn_evidence, f.has_followup_strict, f.htn_evidence_strict
                FROM p3_{short}_{suffix} i
                JOIN p3_{short}_followup f
                  ON f.index_variant = '{variant}' AND f.visit_id = i.visit_id
                  AND f.patient_id = i.patient_id
            """)
    statement = " UNION ALL ".join(expressions)
    con.execute(f"CREATE OR REPLACE TEMP TABLE p3_patient_followup AS {statement}")


def validate_followup(con, flow):
    for site, variant, *_, expected_size in flow:
        rows, patients, invalid = con.execute("""
            SELECT count(*), count(DISTINCT patient_id),
                count(*) FILTER (WHERE has_followup IS NULL OR has_followup_strict IS NULL
                    OR has_followup IS DISTINCT FROM (htn_evidence IS NOT NULL)
                    OR has_followup_strict IS DISTINCT FROM (htn_evidence_strict IS NOT NULL)
                    OR (has_followup_strict AND NOT has_followup))
            FROM p3_patient_followup WHERE site = ? AND index_variant = ?
        """, [site, variant]).fetchone()
        if rows != expected_size or patients != rows or invalid:
            raise ValueError("Follow-up data do not match the selected index cohort.")


def smd(prop1, prop0, var1, var0):
    if any(value is None or not np.isfinite(value) for value in (prop1, prop0, var1, var0)):
        return None
    pooled = np.sqrt((var1 + var0) / 2)
    return float(abs(prop1 - prop0) / pooled) if pooled > 0 else None


def describe_followup(con, p0):
    p1 = import_module("02_p1_tables")
    totals = con.execute("""
        SELECT site, index_variant, count(*) AS indexed,
            count(*) FILTER (WHERE has_followup) AS observed,
            count(*) FILTER (WHERE has_followup_strict) AS strict,
            count(*) FILTER (WHERE has_medication_list) AS listed,
            count(*) FILTER (WHERE has_medication_list AND has_followup) AS listed_with_followup
        FROM p3_patient_followup GROUP BY site, index_variant ORDER BY site, index_variant
    """).fetchall()
    p0.write_csv(
        TABLES / "p3_followup_availability.csv",
        ("site", "index_variant", "patients", "with_followup", "followup_pct",
         "followup_wilson95_low_pct", "followup_wilson95_high_pct",
         "with_strict_followup", "strict_followup_pct", "medication_list_patients",
         "medication_list_with_followup", "medication_list_followup_pct"),
        [(site, variant, n, follow, p1.percent(follow, n), *p1.wilson(follow, n),
          strict, p1.percent(strict, n), listed, listed_followed, p1.percent(listed_followed, listed))
         for site, variant, n, follow, strict, listed, listed_followed in totals],
    )
    comparisons = []
    for _, site in SITES:
        n_follow, n_unfollow = con.execute("""
            SELECT count(*) FILTER (WHERE has_followup), count(*) FILTER (WHERE NOT has_followup)
            FROM p3_patient_followup WHERE site = ? AND index_variant = 'primary'
        """, [site]).fetchone()
        for field in CONTINUOUS_FIELDS:
            data = con.execute(f"""
                SELECT has_followup, count({field}), median({field}),
                    quantile_cont({field}, 0.25), quantile_cont({field}, 0.75),
                    avg({field}), var_samp({field})
                FROM p3_patient_followup WHERE site = ? AND index_variant = 'primary'
                GROUP BY has_followup
            """, [site]).fetchall()
            grouped = {bool(flag): (count, median, q1, q3, mean, variance)
                       for flag, count, median, q1, q3, mean, variance in data}
            yes, med_y, q1_y, q3_y, mean_y, var_y = grouped.get(True, (0, None, None, None, None, None))
            no, med_n, q1_n, q3_n, mean_n, var_n = grouped.get(False, (0, None, None, None, None, None))
            comparisons.append((site, field, "", n_follow, n_unfollow, yes, no,
                                None, None, None, None, med_y, q1_y, q3_y, med_n, q1_n, q3_n,
                                smd(mean_y, mean_n, var_y, var_n)))
        for field in (*CATEGORICAL_FIELDS, *BINARY_FIELDS):
            results = con.execute(f"""
                SELECT has_followup, cast({field} AS VARCHAR) AS level, count(*)
                FROM p3_patient_followup
                WHERE site = ? AND index_variant = 'primary' AND {field} IS NOT NULL
                GROUP BY has_followup, {field} ORDER BY level
            """, [site]).fetchall()
            yes = sum(count for flag, _, count in results if flag)
            no = sum(count for flag, _, count in results if not flag)
            levels = ("true",) if field in BINARY_FIELDS else sorted({level for _, level, _ in results})
            if not levels:
                levels = ("",)
            for level in levels:
                n_y = sum(count for flag, value, count in results if flag and value == level)
                n_n = sum(count for flag, value, count in results if not flag and value == level)
                p_y, p_n = (n_y / yes if yes else None), (n_n / no if no else None)
                diff = smd(p_y, p_n, p_y * (1 - p_y), p_n * (1 - p_n)) if yes and no else None
                comparisons.append((site, field, level, n_follow, n_unfollow, yes, no,
                                    n_y if yes else None, n_n if no else None,
                                    p1.percent(n_y, yes), p1.percent(n_n, no),
                                    None, None, None, None, None, None, diff))
    p0.write_csv(
        TABLES / "p3_followup_baseline.csv",
        ("site", "characteristic", "level", "with_followup_patients", "without_followup_patients",
         "with_followup_nonmissing", "without_followup_nonmissing",
         "with_followup_count", "without_followup_count", "with_followup_pct", "without_followup_pct",
         "with_followup_median", "with_followup_q1", "with_followup_q3",
         "without_followup_median", "without_followup_q1", "without_followup_q3", "absolute_standardized_difference"),
        comparisons,
    )


def export_patient_relation(con, p0):
    path = DERIVED / "p3_patient_followup.parquet"
    pending = path.with_suffix(".parquet.tmp")
    try:
        con.execute(f"""
            COPY (
                SELECT * FROM p3_patient_followup ORDER BY site, index_variant, patient_id
            ) TO {p0.literal(pending)} (FORMAT PARQUET, COMPRESSION ZSTD)
        """)
        pending.replace(path)
    finally:
        pending.unlink(missing_ok=True)


def main():
    p0 = import_module("01_p0_audit")
    setup = import_module("02_p1_setup")
    step = "P3 preflight"
    try:
        setup.verify_project_runtime()
        TABLES.mkdir(parents=True, exist_ok=True)
        before = p0.fingerprint("p3_before")
        with duckdb.connect(str(DERIVED / "epistaxis.duckdb"), read_only=True) as con:
            con.execute(f"SET temp_directory = {p0.literal(DERIVED / 'duckdb_tmp')}")
            setup.register_hosp(con)
            step = "First eligible C3 index"
            flow = create_indices(con)
            step = "MIMIC outpatient BP"
            import_module("04_p3_mimic_followup").run(con)
            step = "MC-MED recorded follow-up"
            import_module("04_p3_mcmed_followup").run(con)
            step = "Follow-up relation and descriptive statistics"
            combine_followup(con)
            validate_followup(con, flow)
            p0.write_csv(
                TABLES / "p3_patient_flow.csv",
                ("site", "index_variant", "c1_visits", "c1_patients", "c3_visits_before_incomplete",
                 "c3_patients_before_incomplete", "incomplete_c3_visits", "c3_eligible_visits",
                 "c3_eligible_patients", "first_c3_index_patients"), flow,
            )
            describe_followup(con, p0)
            export_patient_relation(con, p0)
            step = "Aim 2 modified Poisson and sensitivities"
            import_module("04_p3_statistics").run(con)
        step = "Raw fingerprints"
        if p0.fingerprint("p3_after") != before:
            raise RuntimeError("Source metadata changed.")
    except Exception as error:
        # Raw scanners may include source values in exception messages.
        print(f"P3 stopped at {step} ({type(error).__name__}).", file=sys.stderr)
        return 1
    print("P3 patient-level Aim 2 and aggregate outputs complete.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
