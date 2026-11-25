"""Prepare the approved project environment and verify the P1 hosp allowlist."""

import argparse
import csv
import gzip
from importlib import import_module
import os
from pathlib import Path, PurePosixPath
import sys
import tempfile
from zipfile import ZipFile


ROOT = Path(__file__).resolve().parents[1]
DERIVED = ROOT / "derived"
TABLES = ROOT / "outputs" / "tables"
ARCHIVE = ROOT / "mimic-iv_3.1.zip"
ARCHIVE_PREFIX = "mimic-iv-3.1/"
MANIFEST = "SHA256SUMS.txt"
HOSP_STEMS = (
    "patients", "admissions", "diagnoses_icd", "d_icd_diagnoses",
    "procedures_icd", "d_icd_procedures", "omr",
)
HOSP_PAYLOAD = tuple(f"hosp/{name}.csv.gz" for name in HOSP_STEMS)
P8_LAB_STEMS = ("d_labitems", "labevents")
P8_LAB_PAYLOAD = tuple(f"hosp/{name}.csv.gz" for name in P8_LAB_STEMS)


def verify_project_runtime():
    if Path.cwd() != ROOT or Path(sys.prefix).resolve() != (ROOT / ".venv").resolve():
        raise RuntimeError("P1 requires the project root and the project .venv.")
    local_r = ROOT / ".venv" / "R" / "library"
    if local_r.is_dir():
        os.environ["R_LIBS_USER"] = str(local_r)
    if os.name == "nt":
        os.environ["LC_ALL"] = "Chinese (Simplified)_China.utf8"
        os.environ["LC_CTYPE"] = os.environ["LC_ALL"]
    os.environ["MPLCONFIGDIR"] = str(ROOT / ".venv" / "matplotlib")


def parse_checksums(content, required=HOSP_PAYLOAD):
    checksums = {}
    for line in content.decode("utf-8-sig").splitlines():
        if not line.strip():
            continue
        checksum, name = line.strip().split(maxsplit=1)
        relative = PurePosixPath(name.lstrip("*").strip()).as_posix()
        relative = relative.removeprefix(ARCHIVE_PREFIX)
        if len(checksum) != 64 or int(checksum, 16) < 0 or relative in checksums:
            raise ValueError("Invalid or repeated checksum entry.")
        checksums[relative] = checksum.lower()
    if not set(required) <= checksums.keys():
        raise ValueError("The checksum manifest omits an allowlisted file.")
    return checksums


def extract_hosp(include_p8_labs=False):
    p0 = import_module("01_p0_audit")
    raw_root = ROOT / ARCHIVE_PREFIX
    raw_hosp = raw_root / "hosp"
    for directory in (raw_root, raw_hosp, DERIVED):
        if directory.resolve() != directory:
            raise RuntimeError("The raw destination must remain inside the project.")
    payload = P8_LAB_PAYLOAD if include_p8_labs else HOSP_PAYLOAD
    allowed_files = {raw_root / path for path in (*payload, MANIFEST)}
    protected_sources = sorted(
        set((ROOT / "mc-med").glob("*.csv"))
        | set(ROOT.glob("mimic*.zip"))
        | set(ROOT.glob("mimic*字段.xlsx"))
        | {path for path in (ROOT / "data_raw").rglob("*") if path.is_file()}
        | {path for path in raw_root.rglob("*") if path.is_file()}
    )
    protected_metadata = {path: (path.stat().st_size, path.stat().st_mtime_ns)
                          for path in protected_sources}
    rows = []
    with ZipFile(ARCHIVE) as source, tempfile.TemporaryDirectory(dir=DERIVED) as temporary:
        stage = Path(temporary)
        manifest_bytes = source.read(ARCHIVE_PREFIX + MANIFEST)
        expected = parse_checksums(manifest_bytes, payload)
        pending = []
        manifest = raw_root / MANIFEST
        if include_p8_labs and not manifest.is_file():
            raise RuntimeError("Complete the original P1 hosp extraction before opting in to P8 labs.")
        for relative in payload:
            destination = raw_root / relative
            if destination not in allowed_files or destination.is_symlink():
                raise RuntimeError("Invalid raw destination.")
            if destination.exists():
                with destination.open("rb") as handle:
                    existing = p0.stream_sha256(handle)
                if existing != expected[relative]:
                    raise ValueError("An existing raw file failed checksum verification.")
                actual = existing
            else:
                staged = stage / destination.name
                with source.open(ARCHIVE_PREFIX + relative) as zipped, staged.open("wb") as target:
                    actual = p0.stream_sha256(zipped, target)
                if actual != expected[relative]:
                    raise ValueError("An archive member failed checksum verification.")
                pending.append((staged, destination))
            rows.append((relative, expected[relative], actual, "pass"))
        if manifest.is_symlink():
            raise RuntimeError("Invalid checksum manifest destination.")
        if manifest.exists():
            if manifest.read_bytes() != manifest_bytes:
                raise ValueError("The existing raw manifest differs from the archived member.")
        else:
            staged_manifest = stage / MANIFEST
            staged_manifest.write_bytes(manifest_bytes)
            pending.append((staged_manifest, manifest))
        raw_hosp.mkdir(parents=True, exist_ok=True)
        for staged, destination in pending:
            # Windows rename refuses to replace an already existing raw file.
            staged.rename(destination)
    if protected_metadata != {path: (path.stat().st_size, path.stat().st_mtime_ns)
                              for path in protected_sources}:
        raise RuntimeError("Pre-existing raw files changed during extraction.")
    checksum_file = "p8_hosp_labs_checksum.csv" if include_p8_labs else "p1_hosp_checksum.csv"
    p0.write_csv(TABLES / checksum_file,
                 ("member", "expected_sha256", "actual_sha256", "status"), rows)
    if not include_p8_labs:
        p0.fingerprint("p1_before")
    return len(rows)


def register_hosp(con, include_p8_labs=False):
    p0 = import_module("01_p0_audit")
    for name in (*HOSP_STEMS, *P8_LAB_STEMS) if include_p8_labs else HOSP_STEMS:
        path = ROOT / ARCHIVE_PREFIX / "hosp" / f"{name}.csv.gz"
        with gzip.open(path, mode="rt", encoding="utf-8-sig", newline="") as source:
            header = next(csv.reader(source))
        columns = ", ".join(f"{p0.literal(column)}: 'VARCHAR'" for column in header)
        con.execute(f"""
            CREATE OR REPLACE TEMP VIEW mimic_hosp_{name} AS
            SELECT * FROM read_csv(
                {p0.literal(path.relative_to(ROOT).as_posix())},
                columns = {{{columns}}}, auto_detect = false, header = true,
                delim = ',', strict_mode = true
            )
        """)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--include-p8-labs", action="store_true", help="Verify and extract the two P8 lab files only.")
    args = parser.parse_args()
    step = "project environment"
    try:
        verify_project_runtime()
        if not args.include_p8_labs:
            p0 = import_module("01_p0_audit")
            import_module("statsmodels.api")
            p0.record_environment(None)
        stage = "P8 8.1" if args.include_p8_labs else "P1 1.1"
        step = f"{stage} hosp extraction and hash verification"
        count = extract_hosp(include_p8_labs=args.include_p8_labs)
        print(f"{stage} verified {count} allowlisted hosp payloads.")
    except Exception as error:
        # CSV, SQL and archive parser error details are kept out of logs.
        print(f"Stopped at {step} ({type(error).__name__}).", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
