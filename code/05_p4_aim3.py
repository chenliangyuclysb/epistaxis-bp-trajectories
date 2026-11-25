"""Build the prespecified P4 analysis frames and coordinate the Aim 3 fits."""

import csv
from importlib import import_module
from pathlib import Path
import sys

import duckdb


DERIVED = Path("derived")
TABLES = Path("outputs/tables")
SITES = (("mcmed", "MC-MED"), ("mimic", "MIMIC-IV-ED"))
BASE_COLUMNS = (
    "visit_id", "patient_id", "arrival_ts", "departure_ts", "disposition",
    "is_c1", "is_c2", "dx_match", "cc_match", "is_first_eligible_visit", "incomplete_visit",
    "age_years", "sex", "triage_acuity_group", "any_antithrombotic_use",
    "prior_htn", "hht", "has_medication_list", "prior_htn_same_day_sensitivity",
    "prior_htn_hosp_diagnosis", "prior_htn_ed_diagnosis", "prior_htn_medication",
    "ed_acute_aht", "return_7d_event", "return_30d_event", "return_outcome_eligible",
    "ed_first_acute_aht_minutes",
    "hemostasis_level3", "hemostasis_level3_by_60min", "first_escalation_minutes",
    "landmark_eligible", "hemostasis_level3_after_60min",
    "hemostasis_level3_any_admission_diagnostic",
    "landmark_eligible_any_admission_diagnostic",
    "hemostasis_level3_by_60min_any_admission_diagnostic",
    "hemostasis_level3_after_60min_any_admission_diagnostic",
    "inpatient_ent_outcome_eligible", "inpatient_ent_hemostasis_event",
)
PHENOTYPE_COLUMNS = (
    "bp_phenotype", "sbp_first", "dbp_first", "sbp_last", "dbp_last",
    "time_first_minutes", "time_last_minutes",
)


def prepare_index(con):
    parts = []
    for short, site in SITES:
        dense = ("p.persistent_subtype, p.hypertension_load_sbp" if short == "mcmed" else
                 "NULL::VARCHAR AS persistent_subtype, NULL::DOUBLE AS hypertension_load_sbp")
        parts.append(f"""
            SELECT '{site}' AS site, {", ".join("c." + name for name in BASE_COLUMNS)},
                {", ".join("p." + name for name in PHENOTYPE_COLUMNS)}, {dense},
                epoch(c.departure_ts - c.arrival_ts) / 60.0 AS ed_duration_minutes
            FROM read_parquet('derived/cohort_{short}.parquet') c
            LEFT JOIN read_parquet('derived/phenotypes_{short}.parquet') p USING (visit_id)
        """)
    combined = " UNION ALL ".join(parts)
    con.execute(f"CREATE OR REPLACE TEMP TABLE p4_index AS {combined}")
    duplicates, invalid = con.execute("""
        SELECT count(*) - count(DISTINCT (site, visit_id)),
            count(*) FILTER (
                WHERE patient_id IS NULL OR visit_id IS NULL
                OR is_c1 IS DISTINCT FROM (bp_phenotype IS NOT NULL)
                OR (is_c2 AND (
                    NOT is_c1 OR incomplete_visit OR NOT return_outcome_eligible
                    OR return_7d_event IS NULL OR return_30d_event IS NULL
                ))
                OR (return_7d_event AND NOT return_30d_event)
            )
        FROM p4_index
    """).fetchone()
    if duplicates or invalid:
        raise ValueError("P4 cohort, phenotype or frozen endpoint flags are inconsistent.")
    return con.execute("SELECT * FROM p4_index ORDER BY site, patient_id, arrival_ts, visit_id").df()


def source_fingerprint(p0, label):
    return p0.fingerprint("p4_" + label)


def model_rows(filename):
    with (TABLES / filename).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_table4(p0):
    main_counts = {}
    for row in model_rows("p4_main_group_counts.csv"):
        key = row["endpoint"], row["site"]
        main_counts.setdefault(key, {})[row["phenotype"]] = (row["model_events"], row["model_visits"])
    rows = []
    for row in model_rows("p4_main_estimates.csv"):
        if row["analysis"] != "Main GEE":
            continue
        groups = main_counts[(row["endpoint"], row["site"])]
        covariates = ("Age and sex" if row["site"] == "MC-MED" else
                      "Full covariates, with center" if row["site"] == "Pooled" else
                      "Full covariates, within center")
        rows.append((
            row["endpoint"], row["site"], row["analysis"], "N/T/P phenotype",
            covariates, row["contrast"], row["visits"], row["patients"], row["events"],
            " | ".join(f"{level}: {groups[level][0]}/{groups[level][1]}" for level in ("N", "T", "P")),
            row["estimate"], row["ci95_low"], row["ci95_high"], row["converged"],
            row["fitted_mean_gt_one_n"], row["status"],
        ))
    secondary = model_rows("p4_secondary_estimates.csv")
    group_counts = {}
    for row in secondary:
        key = (row["site"], row["analysis"], row["outcome"], row["exposure"], row["adjustment"])
        groups = group_counts.setdefault(key, {})
        if row["reference"]:
            groups[row["contrast"]] = (row["fitted_contrast_events"], row["fitted_contrast_visits"])
            groups[row["reference"]] = (row["fitted_reference_events"], row["fitted_reference_visits"])
    exposure_labels = {
        "landmark_phenotype": "N60/T60/P60 phenotype",
        "arrival_sbp_per_10": "Arrival SBP per 10 mmHg",
        "bp_phenotype": "N/T/P phenotype",
    }
    covariate_labels = {
        "landmark": "Age, sex, antithrombotic use and prior hypertension",
        "age": "Age only",
        "none": "Crude",
    }
    for row in secondary:
        key = (row["site"], row["analysis"], row["outcome"], row["exposure"], row["adjustment"])
        groups = group_counts[key]
        contrast = (f"{row['contrast']} vs {row['reference']}" if row["reference"] else row["contrast"])
        rows.append((
            row["outcome"], row["site"], row["analysis"], exposure_labels[row["exposure"]],
            covariate_labels[row["adjustment"]], contrast,
            row["modeled_visits"], row["modeled_patients"], row["modeled_events"],
            " | ".join(f"{level}: {groups[level][0]}/{groups[level][1]}" for level in sorted(groups)),
            row["estimate"], row["ci95_low"], row["ci95_high"], row["converged"],
            row["fitted_means_above_one"], row["status"],
        ))
    p0.write_csv(
        TABLES / "p4_table4.csv",
        ("endpoint", "site", "analysis", "exposure", "adjustment", "contrast", "visits", "patients",
         "events", "modeled_events_per_visits_by_group", "rr", "ci95_low", "ci95_high", "converged",
         "fitted_mean_gt_one_n", "status"),
        rows,
    )


def main():
    p0 = import_module("01_p0_audit")
    setup = import_module("02_p1_setup")
    step = "P4 environment and source metadata"
    try:
        setup.verify_project_runtime()
        TABLES.mkdir(parents=True, exist_ok=True)
        before = source_fingerprint(p0, "before")
        with duckdb.connect(str(DERIVED / "epistaxis.duckdb"), read_only=True) as con:
            con.execute(f"SET temp_directory = {p0.literal(DERIVED / 'duckdb_tmp')}")
            setup.register_hosp(con)
            step = "P4 frozen primary cohorts"
            index = prepare_index(con)
            print(step, flush=True)
            step = "P4 main, 30-day, interaction and late-SBP models"
            import_module("05_p4_models").run(index.loc[index.is_c2])
            print(step, flush=True)
            step = "P4 landmark, inpatient procedure and exploratory analyses"
            import_module("05_p4_secondary").run(con, index)
            print(step, flush=True)
            step = "P4 preplanned sensitivity cohort rebuilds"
            variants = import_module("05_p4_sensitivity").build(con, index)
            step = "P4 sensitivity and C1-availability weighting"
            import_module("05_p4_sensitivity_models").run(index, variants)
            print(step, flush=True)
            step = "P4 aggregate Table 4"
            write_table4(p0)
        step = "P4 source metadata verification"
        if source_fingerprint(p0, "after") != before:
            raise RuntimeError("An original source fingerprint changed during P4.")
    except Exception as error:
        # CSV scanner and model diagnostics can contain individual records.
        print(f"P4 stopped at {step} ({type(error).__name__}).", file=sys.stderr)
        return 1
    print("P4 Aim 3 aggregate outputs complete.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
