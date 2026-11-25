"""Run the P0 audit from the project root and export aggregate results."""

import argparse
import csv
from datetime import datetime, timezone
import gzip
import hashlib
from importlib import import_module
import io
from pathlib import Path
import random
import shutil
import subprocess
import sys
from zipfile import ZipFile

import config as cfg


ROOT = Path.cwd()
OUTPUTS = ROOT / "outputs"
TABLES = OUTPUTS / "tables"
DERIVED = ROOT / "derived"
ED_ARCHIVE = ROOT / "mimic-iv-ed_2.2.zip"
HOSP_ARCHIVE = ROOT / "mimic-iv_3.1.zip"
ED_ROOT = ROOT / "data_raw" / "mimic-iv-ed"
HOSP_ROOT = ROOT / "mimic-iv-3.1" / "hosp"
MC_TIMES = {
    "visits": ("Arrival_time", "Roomed_time", "Dispo_time", "Admit_time", "Departure_time"),
    "numerics": ("Time",),
    "orders": ("Order_time", "First_admin_time", "Result_time"),
    "meds": ("Entry_date", "Start_date", "End_date"),
    "pmh": ("Noted_date",),
    "labs": ("Order_time", "Result_time"),
    "rads": ("Order_time", "Result_time"),
}
MIMIC_TIMES = {
    "edstays": ("intime", "outtime"),
    "diagnosis": (),
    "triage": (),
    "vitalsign": ("charttime",),
    "pyxis": ("charttime",),
    "medrecon": ("charttime",),
}
ED_PAYLOAD = {f"ed/{name}.csv.gz" for name in MIMIC_TIMES} | {"README.txt"}
ED_MANIFEST = "SHA256SUMS.txt"
HOSP_REQUIRED = (
    "omr", "patients", "admissions", "diagnoses_icd", "d_icd_diagnoses",
    "procedures_icd", "d_icd_procedures",
)
HOSP_OPTIONAL = ("labevents", "d_labitems")
PACKAGES = ("pandas", "numpy", "duckdb", "pyarrow", "scipy", "statsmodels", "lifelines", "matplotlib")


def write_csv(path, header, rows):
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)
    temporary.replace(path)


def python_package_version(name, executable=None):
    if executable is None:
        try:
            return str(import_module(name).__version__)
        except ImportError:
            return "MISSING"

    process = subprocess.run(
        [executable, "-B", "-c", f"import {name}; print({name}.__version__)"],
        capture_output=True, text=True, check=False,
    )
    return process.stdout.strip() if process.returncode == 0 else "UNAVAILABLE"


def record_environment(statsmodels_python):
    versions = {name: python_package_version(name) for name in PACKAGES}
    lines = [
        f"Python {sys.version.split()[0]}",
        f"Python executable: {sys.executable}",
        f"Random seed: {cfg.RANDOM_SEED}",
        *(f"{name}: {version}" for name, version in versions.items()),
    ]
    if statsmodels_python and versions["statsmodels"] == "MISSING":
        version = python_package_version("statsmodels", statsmodels_python)
        lines.extend([
            f"statsmodels (existing local environment, import checked): {version}",
            f"statsmodels Python executable: {statsmodels_python}",
        ])
    rscript = shutil.which("Rscript")
    if rscript is None:
        lines.extend(["R: UNAVAILABLE", "lcmm: NOT CHECKED"])
    else:
        result = subprocess.run(
            [rscript, "--vanilla", "-e", "writeLines(R.version.string)", "-e",
             'writeLines(paste0("lcmm: ", if (requireNamespace("lcmm", quietly=TRUE)) '
             'as.character(packageVersion("lcmm")) else "MISSING"))'],
            capture_output=True, text=True, check=False,
        )
        lines.extend(result.stdout.strip().splitlines() if result.returncode == 0 else [
            "R: UNAVAILABLE", "lcmm: NOT CHECKED",
        ])
    (OUTPUTS / "env.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def fingerprint(label):
    source_files = set((ROOT / "mc-med").glob("*.csv"))
    source_files.update(ROOT.glob("mimic*.zip"))
    source_files.update(ROOT.glob("mimic*字段.xlsx"))
    if (ROOT / "mimic-iv-3.1").is_dir():
        source_files.update(path for path in (ROOT / "mimic-iv-3.1").rglob("*") if path.is_file())
    rows = []
    for path in sorted(source_files):
        info = path.stat()
        mtime_utc = datetime.fromtimestamp(info.st_mtime, tz=timezone.utc).isoformat()
        rows.append((path.relative_to(ROOT).as_posix(), info.st_size, info.st_mtime_ns, mtime_utc))
    write_csv(OUTPUTS / f"raw_fingerprint_{label}.csv",
              ("relative_path", "bytes", "mtime_ns", "mtime_utc"), rows)
    return rows


def stream_sha256(handle, destination=None):
    digest = hashlib.sha256()
    for block in iter(lambda: handle.read(1024 * 1024), b""):
        digest.update(block)
        if destination is not None:
            destination.write(block)
    return digest.hexdigest()


def extract_ed():
    with ZipFile(ED_ARCHIVE) as archive:
        prefix = "mimic-iv-ed/"
        manifest_bytes = archive.read(prefix + ED_MANIFEST)
        checksums = {}
        for line in manifest_bytes.decode("utf-8-sig").splitlines():
            if not line.strip():
                continue
            expected, relative = line.split(maxsplit=1)
            checksums[relative.lstrip("*")] = expected.lower()
        if not ED_PAYLOAD <= checksums.keys():
            raise ValueError("An allowlisted file is missing from the checksum manifest.")

        rows = []
        for relative in sorted(ED_PAYLOAD):
            path = ED_ROOT / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            temporary = path.with_suffix(path.suffix + ".tmp")
            try:
                with archive.open(prefix + relative) as source, temporary.open("wb") as destination:
                    checksum = stream_sha256(source, destination)
                status = "pass" if checksum == checksums[relative] else "fail"
                rows.append((relative, checksums[relative], checksum, "extracted", status))
                if status == "pass":
                    temporary.replace(path)
            finally:
                temporary.unlink(missing_ok=True)

        # LICENSE is in the manifest but outside the extraction allowlist.
        if "LICENSE.txt" in checksums:
            with archive.open(prefix + "LICENSE.txt") as source:
                checksum = stream_sha256(source)
            status = "pass" if checksum == checksums["LICENSE.txt"] else "fail"
            rows.append(("LICENSE.txt", checksums["LICENSE.txt"], checksum, "archive_only", status))
        unknown_entries = checksums.keys() - ED_PAYLOAD - {"LICENSE.txt"}
        if unknown_entries:
            raise ValueError("Unexpected checksum manifest entry.")
        write_csv(TABLES / "p0_ed_checksum.csv",
                  ("member", "expected_sha256", "actual_sha256", "location", "status"), rows)
        if any(row[-1] == "fail" for row in rows):
            raise ValueError("SHA-256 mismatch.")
        (ED_ROOT / ED_MANIFEST).write_bytes(manifest_bytes)


def audit_hosp():
    names = set()
    if HOSP_ARCHIVE.is_file():
        with ZipFile(HOSP_ARCHIVE) as archive:
            names = {"/".join(Path(entry).parts[-2:]) for entry in archive.namelist()}
    rows = []
    for requirement, stems in (("required", HOSP_REQUIRED), ("optional", HOSP_OPTIONAL)):
        for stem in stems:
            name = f"{stem}.csv.gz"
            rows.append((name, requirement, (HOSP_ROOT / name).is_file(), f"hosp/{name}" in names))
    write_csv(TABLES / "p0_hosp_inventory.csv",
              ("file", "requirement", "present_in_hosp_directory", "listed_in_zip"), rows)


def ident(value):
    return '"' + value.replace('"', '""') + '"'


def literal(value):
    return "'" + str(value).replace("'", "''") + "'"


def csv_columns(path):
    open_csv = gzip.open if path.suffix == ".gz" else Path.open
    with open_csv(path, mode="rt" if path.suffix == ".gz" else "r",
                  encoding="utf-8-sig", newline="") as source:
        return next(csv.reader(source))


def register_and_audit(con):
    table_rows, time_rows = [], []
    sources = [
        *(("mcmed", name, ROOT / "mc-med" / f"{name}.csv", fields)
          for name, fields in MC_TIMES.items()),
        *(("mimic", name, ED_ROOT / "ed" / f"{name}.csv.gz", fields)
          for name, fields in MIMIC_TIMES.items()),
    ]
    ids = ("CSN", "MRN", "stay_id", "subject_id")
    for site, name, path, timestamp_fields in sources:
        print(f"0.4 {site}_{name}", flush=True)
        header = csv_columns(path)
        view = ident(f"{site}_{name}")
        raw_fields = ", ".join(f"{literal(column)}: 'VARCHAR'" for column in header)
        con.execute(f"""
            CREATE OR REPLACE VIEW {view} AS
            SELECT * FROM read_csv(
                {literal(path.relative_to(ROOT).as_posix())},
                columns = {{{raw_fields}}}, auto_detect = false, header = true,
                delim = ',', strict_mode = true
            )
        """)
        clean_ids = [
            f"count(DISTINCT nullif(trim({ident(column)}), ''))" if column in header else "NULL"
            for column in ids
        ]
        base_fields = ", ".join(["count(*)", *clean_ids])
        raw_times = ", ".join(f"nullif(trim({ident(field)}), '') AS t{i}"
                              for i, field in enumerate(timestamp_fields))
        typed_times = ", ".join(f"try_cast(t{i} AS TIMESTAMP) AS ts{i}"
                                for i, _ in enumerate(timestamp_fields))
        time_counts = ", ".join(
            f"count(t{i}), count(ts{i}), year(min(ts{i})), year(max(ts{i}))"
            for i, _ in enumerate(timestamp_fields)
        )
        additions = f", {time_counts}" if timestamp_fields else ""
        reader = (
            f"(SELECT *, {typed_times} FROM (SELECT *, {raw_times} FROM {view}))"
            if timestamp_fields else view
        )
        counts = con.execute(f"SELECT {base_fields}{additions} FROM {reader}").fetchone()
        row_count = counts[0]
        table_rows.append((site, name, *counts[:5]))
        for index, field in enumerate(timestamp_fields):
            offset = 5 + 4 * index
            populated, valid, year_start, year_end = counts[offset:offset + 4]
            percent = valid * 100 / populated if populated else None
            status = ("no_populated_values" if percent is None else
                      "pass" if percent >= cfg.TIMESTAMP_MIN_PARSE_PERCENT else "fail")
            time_rows.append((
                site, name, field, row_count, populated, row_count - populated,
                round((row_count - populated) * 100 / row_count, 6) if row_count else None,
                valid, populated - valid, round(percent, 6) if percent is not None else None,
                year_start, year_end, status,
            ))
    write_csv(TABLES / "p0_table_inventory.csv",
              ("site", "table", "rows", *(f"unique_{name.lower()}" for name in ids)), table_rows)
    write_csv(TABLES / "p0_timestamp_audit.csv",
              ("site", "table", "field", "total_rows", "populated", "missing", "missing_pct",
               "parsed", "failed", "parsed_pct_of_populated", "min_year", "max_year",
               "g0_status"), time_rows)


def prepare_candidates(con):
    con.execute("""
        CREATE OR REPLACE TEMP TABLE p0_mcmed_candidates AS
        WITH flagged AS (
            SELECT *,
                try_cast(Arrival_time AS TIMESTAMP) AS arrival_ts,
                try_cast(Departure_time AS TIMESTAMP) AS departure_ts,
                try_cast(Age AS DOUBLE) AS age_years,
                coalesce(regexp_matches(trim(Dx_ICD10), $dx), false) AS dx_match,
                coalesce(regexp_matches(CC, $cc), false) AS cc_match
            FROM mcmed_visits
        )
        SELECT * FROM flagged WHERE dx_match OR cc_match
    """, {"dx": cfg.MCMED_DX_PATTERN, "cc": cfg.EPIS_CC_PATTERN})
    con.execute("""
        CREATE OR REPLACE TEMP TABLE p0_mcmed_index AS
        SELECT *
        FROM p0_mcmed_candidates
        WHERE age_years >= $adult AND arrival_ts IS NOT NULL AND MRN IS NOT NULL
            AND NOT coalesce(regexp_matches(trim(Dx_ICD10), $trauma), false)
            AND NOT coalesce(regexp_matches(CC, $postop_cc), false)
            AND NOT coalesce(regexp_matches(Dx_name, $postop_dx), false)
        QUALIFY row_number() OVER (PARTITION BY MRN ORDER BY arrival_ts, CSN) = 1
    """, {
        "adult": cfg.MIN_ADULT_AGE, "trauma": cfg.TRAUMA_DX_PATTERN,
        "postop_cc": cfg.POST_OP_CC_PATTERN, "postop_dx": cfg.POST_OP_DX_PATTERN,
    })
    con.execute("""
        CREATE OR REPLACE TEMP TABLE p0_mimic_dx_flags AS
        WITH coded AS (
            SELECT stay_id, try_cast(seq_num AS INTEGER) AS sequence,
                coalesce((icd_version = '10' AND regexp_matches(trim(icd_code), $dx10))
                    OR (icd_version = '9' AND trim(icd_code) = $dx9), false) AS epistaxis,
                coalesce(icd_version = '10' AND regexp_matches(trim(icd_code), $trauma), false)
                    AS trauma,
                coalesce(regexp_matches(icd_title, $postop_dx), false) AS postoperative
            FROM mimic_diagnosis
        )
        SELECT stay_id, bool_or(epistaxis) AS dx_any_match,
            bool_or(epistaxis AND sequence = $primary) AS dx_match,
            bool_or(trauma) AS trauma_dx_match,
            bool_or(postoperative) AS post_op_dx_match
        FROM coded GROUP BY stay_id
    """, {
        "dx10": cfg.MIMIC_DX10_PATTERN, "dx9": cfg.MIMIC_DX9_CODE,
        "trauma": cfg.TRAUMA_DX_PATTERN, "postop_dx": cfg.POST_OP_DX_PATTERN,
        "primary": cfg.MIMIC_PRIMARY_DX_SEQ,
    })
    con.execute("""
        CREATE OR REPLACE TEMP TABLE p0_mimic_candidates AS
        WITH complaints AS (
            SELECT stay_id,
                bool_or(coalesce(regexp_matches(chiefcomplaint, $cc), false)) AS cc_match,
                bool_or(coalesce(regexp_matches(chiefcomplaint, $postop_cc), false)) AS postop_cc
            FROM mimic_triage GROUP BY stay_id
        )
        SELECT ed.*,
            try_cast(ed.intime AS TIMESTAMP) AS arrival_ts,
            try_cast(ed.outtime AS TIMESTAMP) AS departure_ts,
            coalesce(dx.dx_match, false) AS dx_match,
            coalesce(c.cc_match, false) AS cc_match,
            coalesce(c.postop_cc, false) AS postop_cc
        FROM mimic_edstays ed
        LEFT JOIN p0_mimic_dx_flags dx USING (stay_id)
        LEFT JOIN complaints c USING (stay_id)
        WHERE dx.dx_match OR c.cc_match
    """, {"cc": cfg.EPIS_CC_PATTERN, "postop_cc": cfg.POST_OP_CC_PATTERN})
    con.execute("""
        CREATE OR REPLACE TEMP TABLE p0_mimic_index AS
        SELECT c.*
        FROM p0_mimic_candidates c
        LEFT JOIN p0_mimic_dx_flags d USING (stay_id)
        WHERE c.arrival_ts IS NOT NULL AND c.subject_id IS NOT NULL AND NOT c.postop_cc
            AND NOT coalesce(d.trauma_dx_match, false)
            AND NOT coalesce(d.post_op_dx_match, false)
        QUALIFY row_number() OVER (PARTITION BY c.subject_id ORDER BY c.arrival_ts, c.stay_id) = 1
    """)


def open_database():
    duckdb = import_module("duckdb")
    con = duckdb.connect(str(DERIVED / "epistaxis.duckdb"))
    con.execute(f"SET temp_directory = {literal(DERIVED / 'duckdb_tmp')}")
    con.execute("SET preserve_insertion_order = false")
    return con


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prepare-only", action="store_true",
                        help="Record the environment, extract/check ED and inventory hosp.")
    parser.add_argument("--statsmodels-python", help="Existing local interpreter for the version check.")
    args = parser.parse_args()
    random.seed(cfg.RANDOM_SEED)
    for directory in (OUTPUTS, TABLES, DERIVED):
        directory.mkdir(parents=True, exist_ok=True)
    step = "0.2 initial fingerprints"
    before = fingerprint("before")
    try:
        step = "0.1 environment"
        record_environment(args.statsmodels_python)
        step = "0.3 allowlisted ED extraction and SHA-256 verification"
        extract_ed()
        step = "0.8 hosp inventory"
        audit_hosp()
        if not args.prepare_only:
            step = "0.4 DuckDB view, identifier and timestamp audit"
            with open_database() as con:
                register_and_audit(con)
                step = "0.5/0.6/0.7 temporary candidate tables"
                prepare_candidates(con)
                for step, module in (
                    ("0.6 blood-pressure audit", "01_p0_bp"),
                    ("0.5 reconnaissance validation", "01_p0_recon_check"),
                    ("0.7 candidate dictionary tables", "01_p0_dictionaries"),
                ):
                    print(step, flush=True)
                    import_module(module).run(con)
    except Exception as error:
        # Scanner errors can include values from the source rows.
        print(f"P0 stopped at {step} ({type(error).__name__}).", file=sys.stderr)
        return 1
    finally:
        after = fingerprint("after")
        if before != after:
            print("P0 raw file fingerprint comparison failed.", file=sys.stderr)
    if before != after:
        return 1
    print("P0 raw file fingerprints match.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
