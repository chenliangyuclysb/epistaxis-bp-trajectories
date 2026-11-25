"""Append P8 aggregate provenance and check every numerical value in the supplement."""

import argparse
import ast
import csv
from decimal import Decimal, InvalidOperation
from functools import lru_cache
from importlib import import_module
from io import StringIO
from pathlib import Path
import re
import sys


ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "outputs/key_numbers.csv"
SUPPLEMENT = ROOT / "manuscript/supplement/supplement_v1.md"
CHECK = SUPPLEMENT.with_name("number_check_supplement.csv")
ANALYSIS_SOURCES = (
    "p8_posthoc_severity.csv", "p8_posthoc_severity_descriptive.csv",
    "p8_pyxis_local_hemostatic_audit.csv", "p8_hemoglobin_audit.csv", "p8_hemoglobin_coverage.csv",
)
CONFIG_NUMBERS = {
    "WIDE_DUP_MIN": ("MC_DENSE_ALT_DUP_MINUTES", ()),
    "LCMM_K_MIN": ("P2_LCMM_CLASS_COUNTS", (0,)),
    "LCMM_K_MAX": ("P2_LCMM_CLASS_COUNTS", (3,)),
    "LCMM_MIN_PCT": ("P2_LCMM_MIN_CLASS_PROPORTION", ()),
    "LCMM_POSTERIOR": ("P2_LCMM_MIN_MEAN_POSTERIOR", ()),
    "RANDOM_SEED": ("RANDOM_SEED", ()),
    "HB_INCREMENT": ("P8_HEMOGLOBIN_MODEL_INCREMENT_G_DL", ()),
    "HB_MIN_AVAIL": ("P8_MIMIC_HEMOGLOBIN_MIN_AVAILABLE_PERCENT", ()),
    "STRATUM_EVENT_MIN": ("P8_MIN_STRATUM_EVENTS", ()),
    "EPV_MIN": ("P8_MIN_EVENTS_PER_PARAMETER", ()),
}
ENV_NUMBERS = {
    "PYTHON_VERSION": "Python ", "R_VERSION": "R version ", "DUCKDB_VERSION": "duckdb: ",
    "PANDAS_VERSION": "pandas: ", "NUMPY_VERSION": "numpy: ", "PYARROW_VERSION": "pyarrow: ",
    "SCIPY_VERSION": "scipy: ", "STATSMODELS_VERSION": "statsmodels: ", "LIFELINES_VERSION": "lifelines: ",
    "NLME_VERSION": "nlme: ", "LCMM_VERSION": "lcmm: ", "LME4_VERSION": "lme4: ", "LMERTEST_VERSION": "lmerTest: ",
}
OLD_ALIASES = {"KNOT_LOW": "K07454", "KNOT_MEDIAN": "K07455", "KNOT_HIGH": "K07456"}
KEY = re.compile(r"\{(K\d+)\}")
PLACEHOLDER = re.compile(r"\[\[([^\[\]\r\n]+)\]\]")
NUMBER = re.compile(r"(?<![\w.])(?:R\d+\.\d+|v\d+(?:\.\d+)+|\d+(?:\.\d+){2,}|\d+\.\d+-\d+|"
                    r"[+−-]?(?:(?:\d{1,3}(?:,\d{3})+|\d+)(?:\.\d+)?|\.\d+)(?:[eE][+-]?\d+)?%?)")


@lru_cache(maxsize=None)
def csv_rows(relative):
    path = (ROOT / relative).resolve()
    if not path.is_relative_to(ROOT / "outputs") or path.suffix != ".csv":
        raise ValueError("Numeric CSV sources must be aggregate project outputs.")
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        header = tuple(reader.fieldnames)
        line = reader.line_num + 1
        rows = {}
        for row in reader:
            rows[line] = row
            line = reader.line_num + 1
    return header, rows


@lru_cache(maxsize=None)
def text_lines(relative):
    path = (ROOT / relative).resolve()
    if not path.is_relative_to(ROOT) or path.suffix not in (".md", ".py", ".txt"):
        raise ValueError("Text locators must refer to in-project protocol or environment text.")
    return path.read_text(encoding="utf-8-sig").splitlines()


def decimal(value):
    try:
        converted = Decimal(value.strip().rstrip("%").replace(",", "").replace("−", "-"))
    except (InvalidOperation, AttributeError):
        return None
    return converted if converted.is_finite() else None


def locator(relative, fragment):
    line, text = next((i, text) for i, text in enumerate(text_lines(relative), 1) if fragment in text)
    return line, text.strip()


def method_rows():
    cfg = import_module("config")
    config = "code/config.py"
    names = {}
    for node in ast.parse("\n".join(text_lines(config))).body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                targets = target.elts if isinstance(target, ast.Tuple) else (target,)
                names.update((target.id, node.lineno) for target in targets if isinstance(target, ast.Name))
    for symbol, (name, indices) in CONFIG_NUMBERS.items():
        raw = getattr(cfg, name)
        for index in indices:
            raw = raw[index]
        shown = f"{100 * raw:.1f}%" if symbol == "LCMM_MIN_PCT" else str(raw)
        if symbol == "HB_MIN_AVAIL":
            shown = f"{raw:.1f}%"
        line = names[name]
        yield symbol, shown, config, line, name, text_lines(config)[line - 1].strip()

    sap = "docs/01_研究方案_SAP.md" if (ROOT / "docs/01_研究方案_SAP.md").exists() else "docs/SAP_v1.0_zh.md"
    documented = (
        ("MC_YEAR_START", "2020", "MC-MED（2020–2022）"), ("MC_YEAR_END", "2022", "MC-MED（2020–2022）"),
        ("MIMIC_YEAR_START", "2011", "MIMIC-IV-ED 2.2（2011–2019）"),
        ("MIMIC_YEAR_END", "2019", "MIMIC-IV-ED 2.2（2011–2019）"),
        ("EPISTAXIS_ICD10", "R04.0", "`Dx_ICD10` 以 R04.0"),
        ("EPISTAXIS_ICD9", "784.7", "或 ICD-9 7847"),
    )
    for symbol, value, fragment in documented:
        line, text = locator(sap, fragment)
        yield symbol, value, sap, line, "frozen scope or diagnosis", text

    env = "outputs/env.txt"
    for symbol, prefix in ENV_NUMBERS.items():
        line, raw = next((i, text.strip()) for i, text in enumerate(text_lines(env), 1) if text.startswith(prefix))
        yield symbol, raw.removeprefix(prefix).split()[0], env, line, "software version", raw


def make_entry(fields, key, description, value, path, line, column, cell, extra=None):
    entry = dict.fromkeys(fields, "")
    entry.update(key_number_id=key, description=description, value=value, source_file=path,
                 source_line=str(line), source_column=column, source_cell_value=cell)
    if extra:
        entry.update(extra)
    return entry


class Ledger:
    def __init__(self, fields, rows):
        self.fields = fields
        self.items = {r["key_number_id"]: r for r in rows}
        self.next = max(int(key[1:]) for key in self.items) + 1
        self.pending = []
        self.locators = {}
        self.aliases = OLD_ALIASES.copy()
        for key, row in self.items.items():
            position = (row["source_file"], row["source_line"], row["source_column"], row["value"])
            self.locators.setdefault(position, key)
            for prefix in ("Manuscript method ", "P8 method "):
                if row["description"].startswith(prefix):
                    self.aliases[row["description"].removeprefix(prefix).split()[0]] = key

    def add(self, description, value, path, line, column, source_cell, extra=None):
        position = (path, str(line), column, value)
        if position in self.locators:
            return self.locators[position]
        key = f"K{self.next:05d}"
        self.next += 1
        item = make_entry(self.fields, key, description, value, path, line, column, source_cell, extra)
        self.items[key] = item
        self.pending.append(item)
        self.locators[position] = key
        return key

    def ingest(self):
        for alias, value, path, line, column, cell in method_rows():
            self.aliases[alias] = self.add(f"P8 method {alias}", value, path, line, column, cell)
        for name in ANALYSIS_SOURCES:
            path = f"outputs/tables/{name}"
            _, rows = csv_rows(path)
            for line, row in rows.items():
                context = ", ".join(row.get(label, "") for label in ("analysis", "center", "site", "model", "contrast", "phenotype")
                                    if row.get(label))
                for column, value in row.items():
                    if decimal(value) is None:
                        continue
                    label = f"P8 post hoc {Path(name).stem} {context} {column}".strip()
                    self.add(label, value, path, line, column, value, ci_reference(row, path, line, column))

    def render(self, expression):
        if expression.startswith("CSV|"):
            _, name, line, column, shown = expression.split("|", 4)
            if name not in ANALYSIS_SOURCES:
                raise ValueError("Post hoc prose may cite only the aggregate results allowlist.")
            path = f"outputs/tables/{name}"
            value = csv_rows(path)[1][int(line)][column]
            key = self.locators[path, line, column, value]
        else:
            key = self.aliases[expression]
            shown = self.items[key]["value"]
            if expression in OLD_ALIASES:
                shown = f"{decimal(shown):g}".removesuffix(".0")
            if expression == "MATCH_BOOT_REPS":
                shown = f"{int(shown):,}"
        return f"{shown} {{{key}}}"


def ci_reference(row, path, line, column):
    if column not in ("rr", "estimate"):
        return None
    for lo, hi in (("ci95_low", "ci95_high"), ("ci_low", "ci_high")):
        if decimal(row.get(lo)) is not None and decimal(row.get(hi)) is not None:
            return {"ci95_low": row[lo], "ci95_high": row[hi],
                    "ci95_low_reference": f"{path}:{line}:{lo}", "ci95_high_reference": f"{path}:{line}:{hi}"}
    return None


def matches_value(shown, ledger_value):
    if shown == ledger_value:
        return True
    left, right = decimal(shown), decimal(ledger_value)
    if left is None or right is None:
        return False
    return abs(left - right) <= Decimal(5).scaleb(left.as_tuple().exponent - 1)


def verify_citation(key, shown, ledger):
    row = ledger.items[key]
    source, line = row["source_file"], int(row["source_line"])
    if source.endswith(".csv"):
        valid = csv_rows(source)[1][line][row["source_column"]] == row["source_cell_value"]
    else:
        valid = row["source_cell_value"] in text_lines(source)[line - 1]
    if not valid or not matches_value(shown, row["value"]):
        raise ValueError(f"Supplementary numeric check failed for {key}.")
    return (key, shown, row["value"], "Yes", source, line, row["source_column"])


def check_numbers(rendered, ledger):
    checks = []
    for line, text in enumerate(rendered.splitlines(), 1):
        expected = {m.start(): m for m in KEY.finditer(text)}
        for token in NUMBER.finditer(text):
            shown = token[0]
            prefix = text[:token.start()]
            if prefix.endswith("Supplementary Data ") and shown == "1":
                continue
            if shown.startswith("-") and prefix.endswith("}"):
                shown = shown[1:]
            key_start = token.end() + len(text[token.end():]) - len(text[token.end():].lstrip(" "))
            citation = expected.pop(key_start, None)
            if citation is None:
                raise ValueError(f"Unlabeled numerical token in supplementary prose line {line}.")
            checks.append(verify_citation(citation[1], shown, ledger))
        if expected:
            raise ValueError(f"A supplementary ledger label is detached in line {line}.")
    return checks


def write_ledger(prefix, ledger):
    if not ledger.pending:
        return
    if not prefix.endswith(b"\n") or LEDGER.read_bytes() != prefix:
        raise RuntimeError("Refusing to replace or reformat the existing ledger prefix.")
    stream = StringIO(newline="")
    csv.DictWriter(stream, ledger.fields).writerows(ledger.pending)
    stage = LEDGER.with_suffix(".csv.tmp")
    try:
        stage.write_bytes(prefix + stream.getvalue().encode("utf-8"))
        stage.replace(LEDGER)
    finally:
        stage.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-only", action="store_true", help="Check and resolve in memory without changing any file.")
    args = parser.parse_args()
    step = "project environment"
    try:
        import_module("02_p1_setup").verify_project_runtime()
        step = "append-only table and methods provenance"
        before = LEDGER.read_bytes()
        fields, rows = csv_rows(LEDGER.relative_to(ROOT).as_posix())
        ledger = Ledger(fields, rows.values())
        ledger.ingest()
        step = "supplement number check"
        current = SUPPLEMENT.read_text(encoding="utf-8-sig")
        prose = PLACEHOLDER.sub(lambda match: ledger.render(match[1]), current)
        checks = check_numbers(prose, ledger)
        if not args.check_only:
            step = "atomic ledger append and supplementary output"
            write_ledger(before, ledger)
            assert LEDGER.read_bytes().startswith(before)
            if prose != current:
                SUPPLEMENT.write_text(prose, encoding="utf-8", newline="\n")
            import_module("01_p0_audit").write_csv(CHECK, ("Ledger tag", "In supplement", "Ledger value", "Consistent",
                                                            "Source file", "Source line", "Source column"), checks)
    except Exception as error:
        if isinstance(error, ValueError) and str(error).startswith(("Unlabeled numerical token", "Supplementary numeric check failed",
                                                                       "A supplementary ledger label")):
            print(str(error), file=sys.stderr)
        print(f"P8 ledger stopped at {step} ({type(error).__name__}).", file=sys.stderr)
        return 1
    print(f"P8 checked {len(checks)} numerical occurrences; {len(ledger.pending)} append-only ledger entries.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
