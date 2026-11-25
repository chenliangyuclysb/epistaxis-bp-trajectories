"""Run the G4-approved S3.3 sensitivity without changing P4 artifacts."""

import csv
from importlib import import_module
from pathlib import Path
import sys

import duckdb


TABLES = Path("outputs/tables")
COHORT_FIELDS = (
    "visit_id", "patient_id", "age_years", "sex", "any_antithrombotic_use", "prior_htn", "triage_acuity_group",
    "prior_epistaxis_365d", "return_7d_event",
)


def selected_c2(con, models):
    parts = [f"""
        SELECT '{label}' AS site, {', '.join('c.' + name for name in COHORT_FIELDS)}, p.bp_phenotype
        FROM read_parquet('derived/cohort_{short}.parquet') c
        LEFT JOIN read_parquet('derived/phenotypes_{short}.parquet') p USING (visit_id)
        WHERE c.is_c2
    """ for short, label in zip(("mcmed", "mimic"), models.SITES)]
    cohort = con.execute(" UNION ALL ".join(parts) + " ORDER BY site, patient_id, visit_id").df()
    if cohort.empty or cohort.duplicated(["site", "visit_id"]).any() or not cohort.bp_phenotype.isin(models.LEVELS).all():
        raise ValueError("The frozen C2 phenotype cohort is inconsistent.")
    return cohort


def p4_references():
    with (TABLES / "p4_main_estimates.csv").open(encoding="utf-8-sig", newline="") as stream:
        return {row["site"]: row for row in csv.DictReader(stream)
                if row["endpoint"] == "7-day return" and row["analysis"] == "Main GEE" and row["contrast"] == "P vs T"}


def run(con):
    models = import_module("05_p4_models")
    frame = selected_c2(con, models)
    references = p4_references()
    rows = []
    mimic = models.SITES[1]
    for site, subset, adjustment in (("Pooled", frame, "pooled"), (mimic, frame.loc[frame.site.eq(mimic)], "mimic")):
        fit = models.fit_model(subset, "return_7d_event", continuous=("prior_epistaxis_365d",), adjustment=adjustment)
        effect = fit.estimate("bp_phenotype=P")
        reference = references[site]
        if effect["status"] != "Estimated":
            raise RuntimeError("The G4-approved S3.3 consistency estimate is unavailable.")
        if fit.n != int(reference["visits"]) or fit.events != int(reference["events"]) or fit.patients != int(reference["patients"]):
            raise RuntimeError("The supplemental model did not preserve P4 model denominators.")
        observed_reference = int(reference["exposed_visits"]), int(reference["exposed_events"])
        observed_contrast = int(reference["reference_visits"]), int(reference["reference_events"])
        if fit.group_counts["P"] != observed_reference or fit.group_counts["T"] != observed_contrast:
            raise RuntimeError("The supplemental model changed P/T counts from the P4 fit.")
        if fit.n_gt_one or not fit.converged:
            if fit.cross_check is None or fit.cross_check.estimate("bp_phenotype=P")["status"] != "Estimated":
                raise RuntimeError("The supplemental fit requires a successful Poisson cross-check.")
            check_effect = fit.cross_check.estimate("bp_phenotype=P")
        else:
            check_effect = dict.fromkeys(("estimate", "ci95_low", "ci95_high"))
        rows.append((site, "P vs T", fit.method, int(subset.prior_epistaxis_365d.sum()), fit.missing_counts["prior_epistaxis_365d"],
                     *fit.group_counts["P"], *fit.group_counts["T"], reference["estimate"], reference["ci95_low"], reference["ci95_high"],
                     *(effect[key] for key in ("estimate", "ci95_low", "ci95_high", "log_estimate", "se", "status")),
                     *(check_effect[key] for key in ("estimate", "ci95_low", "ci95_high")), *fit.diagnostics()))
    models.write_csv("p5_s33_consistency.csv", (
        "site", "contrast", "method", "known_prior_epistaxis_visits", "missing_prior_epistaxis_visits",
        "p_visits", "p_events", "t_visits", "t_events", "p4_main_estimate", "p4_main_ci95_low", "p4_main_ci95_high",
        "estimate", "ci95_low", "ci95_high", "log_estimate", "se", "status",
        "poisson_crosscheck_estimate", "poisson_crosscheck_ci95_low", "poisson_crosscheck_ci95_high", *models.MODEL_FIELDS,
    ), rows)


def main():
    try:
        import_module("02_p1_setup").verify_project_runtime()
        with duckdb.connect() as con:
            run(con)
    except Exception as error:
        # Omit exception messages, which may contain protected source values.
        print(f"P5 5.0 stopped ({type(error).__name__}).", file=sys.stderr)
        return 1
    print("P5 5.0 S3.3 consistency models completed and preserved P4 counts.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
