"""Compare the audited Aim 1 LMM with the G5 Satterthwaite model."""

import csv
from importlib import import_module
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

import pyarrow as pa
import pyarrow.csv as arrow_csv
import pyarrow.parquet as parquet

import config as cfg


TERMS = ("estimate", "std_error", "df", "ci_low", "ci_high", "p_value")
DIAG_FIELDS = (
    "modeled_reading_n", "modeled_visit_n", "modeled_patient_n", "spline_knot_1_hours", "spline_knot_2_hours",
    "converged", "singular", "optimizer", "warning_count", "reml_loglik", "aic", "residual_sd",
    "patient_intercept_sd", "visit_intercept_sd", "visit_slope_sd_mmhg_per_hour", "visit_intercept_slope_corr",
    "joint_interaction_f", "joint_interaction_numerator_df", "joint_interaction_denominator_df", "joint_interaction_p_value",
)


def rows(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def fit_models(runner):
    source = parquet.read_table(runner.INPUT, columns=list(runner.MODEL_COLUMNS))
    buffer = pa.BufferOutputStream()
    arrow_csv.write_csv(source, buffer)
    with TemporaryDirectory(dir=runner.DERIVED, prefix="p6_lmm_aggregates_") as temp:
        stage = Path(temp)
        executable = runner.r_command(stage)[0]
        command = [executable, "--vanilla", str(runner.ROOT / "code" / "07_p6_lmm.R"), str(runner.R_LIBRARY),
                   str(runner.R_CODE), str(stage), *map(str, cfg.P2_MODEL_HOURS), str(cfg.P2_MODEL_SPLINE_DF),
                   str(cfg.P2_MODEL_CONFIDENCE_LEVEL), str(cfg.RANDOM_SEED)]
        child = subprocess.run(command, env=runner.r_environment(), cwd=runner.ROOT, input=buffer.getvalue().to_pybytes(),
                               capture_output=True, timeout=cfg.P2_R_TIMEOUT_SECONDS, check=False)
        if child.returncode or b"P6_LMM_COMPLETED" not in child.stdout.splitlines():
            raise RuntimeError("The project-local lme4/lmerTest model was unavailable.")
        return rows(stage / "coefficients.csv"), rows(stage / "diagnostics.csv")[0]


def compare(runner, new_terms, new_diag):
    old_terms = {row["term"]: row for row in rows(runner.OUTPUT / "p2_mcmed_lmm_coefficients.csv")
                 if row["specification"] == "pre_acute_dose"}
    old_diag = next(row for row in rows(runner.OUTPUT / "p2_mcmed_lmm_diagnostics.csv")
                    if row["specification"] == "pre_acute_dose")
    old_diag["converged"] = str(old_diag["status"] == "converged").upper()
    new_by_term = {row["term"]: row for row in new_terms}
    if set(old_terms) != set(new_by_term) or any(old_diag[name] != new_diag[name] for name in DIAG_FIELDS[:5]):
        raise ValueError("The Satterthwaite fit must use the original cohort, spline knots and fixed effects.")
    writer = import_module("05_p4_models").write_csv
    writer("p6_lmm_comparison.csv", ("term", *("original_" + field for field in TERMS),
                                    *("satterthwaite_" + field for field in TERMS)),
           [(name, *(old_terms[name][field] for field in TERMS), *(new_by_term[name][field] for field in TERMS))
            for name in old_terms])
    method_data = (("P2 nlme containment", old_diag["random_effects"], old_diag),
                   ("P6 lmerTest Satterthwaite", "visit_intercept_slope", new_diag))
    writer("p6_lmm_diagnostics.csv", ("model", "random_effects", *DIAG_FIELDS),
           [(model, random, *(record.get(field) for field in DIAG_FIELDS)) for model, random, record in method_data])


def main():
    step = "project environment"
    try:
        import_module("02_p1_setup").verify_project_runtime()
        runner = import_module("03_p2_r_models")
        step = "Satterthwaite model"
        new_terms, new_diag = fit_models(runner)
        step = "P2 coefficient and denominator comparison"
        compare(runner, new_terms, new_diag)
    except Exception as error:
        print(f"P6 LMM stopped at {step} ({type(error).__name__}).", file=sys.stderr)
        return 1
    print("P6 Aim 1 LMM robustness comparison completed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
