"""Replay P1–P6 and compare aggregate table hashes and physical CSV cells."""

import argparse
import csv
from decimal import Decimal, InvalidOperation
import hashlib
from importlib import import_module
import io
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
from time import monotonic


ROOT = Path(__file__).resolve().parents[1]
OUTPUTS = ROOT / "outputs"
DERIVED = ROOT / "derived"
AUDIT = OUTPUTS / "reproducibility"
PRIOR = (
    ("P1", "02_p1_setup.py", "02_p1_cohort.py", "02_p1_flow_figure.py"),
    ("P2", "03_p2_apply_g1.py", "03_p2_phenotypes.py", "03_p2_aim1.py", "03_p2_r_models.py",
     "03_p2_figure2.py", "03_p2_figure3.py", "03_p2_lcmm_figure.py"),
    ("P3", "04_p3_aim2.py"),
    ("P4", "05_p4_aim3.py", "05_p4_figures.py"),
    ("P5", "06_p5_s33.py", "06_p5_aim4.py", "06_p5_figure.py"),
)
FINAL = (("P6", "07_p6_lmm.py", "07_p6_overlap.py", "07_p6_tables.py", "07_p6_figures.py", "07_p6_ledger.py",
          "07_p6_release.py"),)


def covered(stage):
    tables = OUTPUTS / "tables"
    if stage == "prior":
        return sorted((*tables.glob("p[1-5]_*.csv"), *(OUTPUTS / "dictionary_final").glob("*.csv")))
    artifacts = [*tables.glob("p6_*.csv"), *(tables / "final").glob("*.csv"), OUTPUTS / "key_numbers.csv"]
    return sorted(path for path in artifacts if path.is_file())


def files_metadata():
    roots = (ROOT / "mc-med", ROOT / "data_raw", ROOT / "mimic-iv-3.1")
    path_set = set(ROOT.glob("mimic*.zip"))
    for directory in roots:
        path_set.update(path for path in directory.rglob("*") if path.is_file())
    return {path: (path.stat().st_size, path.stat().st_mtime_ns) for path in path_set}


def verify_hosp_inputs():
    # The P1 replay verifies existing members. Initial extraction is a separate setup task.
    setup = import_module("02_p1_setup")
    directory = ROOT / setup.ARCHIVE_PREFIX
    if not all((directory / name).is_file() for name in (*setup.HOSP_PAYLOAD, setup.MANIFEST)):
        raise RuntimeError("P1 replay requires the original allowlisted hosp files to be present.")


def run_commands(plan):
    records = []
    for phase, *scripts in plan:
        for script in scripts:
            start = monotonic()
            result = subprocess.run([sys.executable, "-B", str(ROOT / "code" / script)], cwd=ROOT,
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
            seconds = round(monotonic() - start, 2)
            records.append((phase, f"code/{script}", result.returncode, seconds))
            # Suppress child output, including parser and optimizer diagnostics.
            print(f"{phase} {script}: {'passed' if result.returncode == 0 else 'failed'} ({seconds:.1f} s)", flush=True)
            if result.returncode:
                return records, False
    return records, True


def sha256(content):
    return hashlib.sha256(content).hexdigest() if content is not None else None


def csv_records(content):
    reader = csv.reader(io.StringIO(content.decode("utf-8-sig"), newline=""))
    records = []
    line = 1
    for row in reader:
        records.append((line, row))
        line = reader.line_num + 1
    return records


def classify(before, after):
    try:
        left, right = Decimal(before), Decimal(after)
    except InvalidOperation:
        return "text or structure", None
    if not left.is_finite() or not right.is_finite():
        return "nonfinite representation", None
    difference = right - left
    return "numeric representation" if difference == 0 else "numeric value", difference


def differences(path, prior, current):
    left, right = csv_records(prior), csv_records(current)
    if left == right:
        yield (path, None, None, None, None, "CSV byte formatting only", None)
        return
    headings = (right[0] if right else left[0])[1]
    for index in range(max(len(left), len(right))):
        old_line, row_l = left[index] if index < len(left) else (None, [])
        new_line, row_r = right[index] if index < len(right) else (None, [])
        for col in range(max(len(row_l), len(row_r))):
            a = row_l[col] if col < len(row_l) else ""
            b = row_r[col] if col < len(row_r) else ""
            if a != b:
                name = headings[col] if col < len(headings) else f"column_{col + 1}"
                classification, delta = classify(a, b)
                yield (path, new_line or old_line, name, a, b, classification, delta)


def compare(stage, before):
    after = {path.relative_to(ROOT).as_posix(): path for path in covered(stage)}
    hash_rows, change_rows = [], []
    for key in sorted(before.keys() | after.keys()):
        previous = before[key].read_bytes() if key in before else None
        current = after[key].read_bytes() if key in after else None
        old_hash, new_hash = sha256(previous), sha256(current)
        status = ("Identical" if old_hash == new_hash else "Added" if old_hash is None else
                  "Missing" if new_hash is None else "Changed")
        hash_rows.append((key, old_hash, new_hash, status))
        if status == "Changed":
            change_rows.extend(differences(key, previous, current))
    return hash_rows, change_rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("prior", "final"), required=True)
    args = parser.parse_args()
    step = "environment and source metadata"
    try:
        import_module("02_p1_setup").verify_project_runtime()
        p0 = import_module("01_p0_audit")
        AUDIT.mkdir(parents=True, exist_ok=True)
        DERIVED.mkdir(parents=True, exist_ok=True)
        raw_before = files_metadata()
        if args.stage == "prior":
            verify_hosp_inputs()
        with TemporaryDirectory(dir=DERIVED, prefix="p6_aggregate_reference_") as directory:
            before = {}
            for position, path in enumerate(covered(args.stage)):
                snapshot = Path(directory) / f"{position}.csv"
                snapshot.write_bytes(path.read_bytes())
                before[path.relative_to(ROOT).as_posix()] = snapshot
            step = "stage replay"
            logs, passed = run_commands(PRIOR if args.stage == "prior" else FINAL)
            p0.write_csv(AUDIT / f"p6_{args.stage}_steps.csv", ("phase", "script", "exit_code", "elapsed_s"), logs)
            step = "read-only and aggregate consistency verification"
            raw_after = files_metadata()
            p0.write_csv(AUDIT / f"p6_{args.stage}_raw.csv", ("files", "size_and_mtime_unchanged"),
                         [(len(raw_before), raw_before == raw_after)])
            hashes, cells = compare(args.stage, before)
            p0.write_csv(AUDIT / f"p6_{args.stage}_sha256.csv", ("file", "sha256_before", "sha256_after", "status"), hashes)
            p0.write_csv(AUDIT / f"p6_{args.stage}_changed_cells.csv",
                         ("file", "line", "column", "before", "after", "change_type", "numeric_delta"), cells)
            unchanged = sum(row[-1] == "Identical" for row in hashes)
            print(f"{args.stage}: {unchanged}/{len(hashes)} aggregate files byte-identical; {len(cells)} changed cells.")
            if raw_before != raw_after or not passed:
                return 1
    except Exception as error:
        print(f"P6 replay stopped at {step} ({type(error).__name__}).", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
