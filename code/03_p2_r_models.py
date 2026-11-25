"""Run the prespecified P2 R models with an outcome-free, in-memory CSV pipe."""

import argparse
import csv
import os
from pathlib import Path
import shutil
import subprocess
from tempfile import TemporaryDirectory

import pyarrow as pa
import pyarrow.csv as arrow_csv
import pyarrow.parquet as parquet

import config as cfg


ROOT = Path(__file__).resolve().parents[1]
DERIVED = ROOT / "derived"
OUTPUT = ROOT / "outputs" / "tables"
INPUT = DERIVED / "p2_mcmed_model_input.parquet"
R_LIBRARY = ROOT / ".venv" / "R" / "library"
R_CODE = ROOT / "code" / "03_p2_models.R"
MODEL_COLUMNS = (
    "visit_id", "patient_id", "minutes_from_arrival", "sbp", "prior_htn", "age_years",
    "sex", "has_medication_list", "ed_acute_aht", "ed_first_acute_aht_minutes",
)
OUTPUT_NAMES = (
    "p2_mcmed_lmm_coefficients.csv", "p2_mcmed_lmm_diagnostics.csv",
    "p2_mcmed_lmm_adjusted_curves.csv", "p2_mcmed_lcmm_models.csv",
    "p2_mcmed_lcmm_classes.csv", "p2_mcmed_lcmm_trajectories.csv",
)


def r_command(output_dir):
    executable = shutil.which("Rscript")
    if executable is None or not R_LIBRARY.is_dir():
        raise RuntimeError("Rscript or the local R package library is unavailable.")
    return [
        executable, "--vanilla", str(R_CODE), str(R_LIBRARY), str(output_dir), str(cfg.RANDOM_SEED),
        *(str(value) for value in cfg.P2_MODEL_HOURS), str(cfg.P2_MODEL_SPLINE_DF),
        str(cfg.P2_MODEL_CURVE_GRID_POINTS), str(cfg.P2_MODEL_CONFIDENCE_LEVEL),
        ",".join(map(str, cfg.P2_LCMM_CLASS_COUNTS)), str(cfg.P2_LCMM_MIN_CLASS_PROPORTION),
        str(cfg.P2_LCMM_MIN_MEAN_POSTERIOR), str(cfg.P2_LCMM_GRID_STARTS),
        str(cfg.P2_LCMM_GRID_ITERATIONS), str(cfg.P2_LCMM_FINAL_ITERATIONS),
    ]


def r_environment():
    environment = os.environ.copy()
    environment.pop("R_LIBS", None)
    environment.pop("R_LIBS_SITE", None)
    environment["R_LIBS_USER"] = str(R_LIBRARY)
    if os.name == "nt":
        environment["LC_ALL"] = "Chinese (Simplified)_China.utf8"
        environment["LC_CTYPE"] = environment["LC_ALL"]
    for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        environment[variable] = "1"
    return environment


def run(input_path=INPUT, output_dir=OUTPUT):
    input_path = Path(input_path).resolve()
    output_dir = Path(output_dir).resolve()
    if not input_path.is_relative_to(DERIVED.resolve()) or not output_dir.is_relative_to(ROOT):
        raise RuntimeError("Keep the P2 modeling inputs and outputs inside the project.")
    try:
        data = parquet.read_table(input_path, columns=list(MODEL_COLUMNS))
        buffer = pa.BufferOutputStream()
        arrow_csv.write_csv(data, buffer)
        csv_bytes = buffer.getvalue().to_pybytes()
    except Exception:
        raise RuntimeError("Could not prepare the restricted P2 Parquet input.") from None

    output_dir.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=DERIVED, prefix="p2_model_aggregates_") as temporary:
        stage = Path(temporary)
        try:
            completed = subprocess.run(
                r_command(stage), cwd=ROOT, env=r_environment(), input=csv_bytes,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=cfg.P2_R_TIMEOUT_SECONDS, check=False,
            )
        except subprocess.TimeoutExpired:
            raise RuntimeError("P2 R fitting exceeded its runtime limit.") from None
        if completed.returncode or b"P2_MODELS_COMPLETED" not in completed.stdout.splitlines():
            phases = ("setup", "input", "lmm", "lcmm")
            phase = next((name for name in phases if f"P2_MODELS_FAILED_{name}".encode() in completed.stdout.splitlines()), "unknown")
            raise RuntimeError(f"P2 R fitting failed during {phase}. R output was withheld to protect individual records.")
        for filename in OUTPUT_NAMES:
            result = stage / filename
            target = output_dir / filename
            if result.is_file():
                result.replace(target)
            elif target.is_file():
                target.unlink()
    with (output_dir / "p2_mcmed_lmm_diagnostics.csv").open(newline="", encoding="utf-8") as source:
        lmm = list(csv.DictReader(source))
    with (output_dir / "p2_mcmed_lcmm_models.csv").open(newline="", encoding="utf-8") as source:
        lcmm = list(csv.DictReader(source))
    return lmm, lcmm


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=INPUT)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT)
    arguments = parser.parse_args()
    try:
        diagnostics_lmm, diagnostics_lcmm = run(arguments.input, arguments.output_dir)
    except RuntimeError as error:
        raise SystemExit(str(error)) from None
    except OSError:
        raise SystemExit("Could not execute the project R modeling launcher.") from None
    fit_count = sum(row["status"] == "converged" for row in diagnostics_lmm)
    chosen = next((row["k"] for row in diagnostics_lcmm if row["selected"] == "TRUE"), "none")
    print(f"P2 R models completed. LMM converged {fit_count}/{len(diagnostics_lmm)}. Selected LCMM k={chosen}.")
