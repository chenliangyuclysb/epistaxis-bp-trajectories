"""Build the manuscript number ledger from checked aggregate source locators."""

from collections import defaultdict
import ast
import csv
from decimal import Decimal, InvalidOperation
from functools import lru_cache
from importlib import import_module
from pathlib import Path
import re
import sys


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "outputs" / "tables"
PROVENANCE = TABLES / "final" / "p6_table_sources.csv"
EXTRA_SOURCES = (
    "p2_aim1_km_summary.csv", "p2_aim1_km_timepoints.csv", "p2_aim1_phenotype_htn.csv",
    "p2_aim1_mimic_paired_sbp.csv", "p2_last_bp_to_departure.csv",
    "p3_mimic_omr_flow.csv", "p3_mcmed_followup_audit.csv", "p4_main_estimates.csv",
    "p4_main_interaction.csv", "p4_main_evalue.csv", "p4_sbp_tests.csv", "p4_sbp_knots.csv", "p4_sbp_estimates.csv",
    "p5_matching_flow.csv", "p5_matching_group_means.csv",
    "p2_aim1_prevalence_ratio.csv",
)
DESCRIPTORS = {
    "center", "site", "site_or_scope", "analysis", "measure", "metric", "characteristic", "contrast", "outcome",
    "group", "pheno", "phenotype", "prior_htn", "population", "comparison", "scenario", "exposure", "specification",
    "endpoint", "model", "section", "treatment_definition", "minutes", "time_min", "time_minutes",
    "minutes_from_arrival", "level", "cohort", "term_or_contrast", "test",
}
NUMERIC_TOKEN = re.compile(r"(?<![\w.,])[-+−]?(?:\d{1,3}(?:,\d{3})+|\d+)(?:\.\d+)?(?:[eE][+-]?\d+)?")
GROUP_TOKEN = re.compile(r"([A-Za-z0-9_-]+):\s*(\d[\d,]*)/(\d[\d,]*)")
INTERVAL_COLUMNS = {
    "estimate": (("ci95_low", "ci95_high"), ("ci_low", "ci_high")),
    "prevalence_ratio": (("ci95_low", "ci95_high"),),
    "estimate_mmhg": (("ci95_low", "ci95_high"), ("ci_low", "ci_high")),
    "rr": (("ci95_low", "ci95_high"),),
    "crude_rr": (("katz_ci95_low", "katz_ci95_high"),),
    "mean_sbp": (("ci95_low", "ci95_high"),),
    "recovered_pct": (("pointwise95_low_pct", "pointwise95_high_pct"),),
    "mean_change_mmhg": (("cluster_robust_mean95_low_mmhg", "cluster_robust_mean95_high_mmhg"),),
    "original_estimate": (("original_ci_low", "original_ci_high"),),
    "satterthwaite_estimate": (("satterthwaite_ci_low", "satterthwaite_ci_high"),),
    "treated_minus_control_change": (("ci95_low", "ci95_high"),),
    "set_weighted_percent_normal": (("ci95_low_pct", "ci95_high_pct"),),
}

# Manuscript methods numbers retain exact source locations in the frozen configuration.
METHOD_NUMBERS = {
    "ADULT_AGE": "MIN_ADULT_AGE", "C1_BP_MIN": "C1_MIN_BP_READINGS", "C1_SPAN": "C1_MIN_SPAN_MINUTES",
    "EPISODE_WINDOW": "RETURN_7D_HOURS", "RETURN7": "RETURN_7D_HOURS", "RETURN30": "RETURN_30D_HOURS",
    "BP_LIMIT": "ED_MAX_HOURS", "CHART_DUP_WINDOW": "BP_CHART_TRIAGE_DUP_MINUTES",
    "NEAR_DUP_TIME": "G0_BP_DUP_MINUTES", "NEAR_DUP_SBP": "BP_CROSS_SOURCE_MAX_SBP_DIFF",
    "SBP_MIN": "BP_SBP_MIN", "SBP_MAX": "BP_SBP_MAX", "DBP_MIN": "BP_DBP_MIN", "DBP_MAX": "BP_DBP_MAX",
    "PULSE_PRESSURE": "BP_MIN_PULSE_PRESSURE", "ELEVATED_SBP": "BP_ELEVATED_SBP", "ELEVATED_DBP": "BP_ELEVATED_DBP",
    "SEVERE_SBP": "BP_SEVERE_SBP", "SEVERE_DBP": "BP_SEVERE_DBP", "ARV_MIN_READINGS": "BP_MIN_ARV_READINGS",
    "MIMIC_FOLLOWUP_FIRST": "P3_MIMIC_FOLLOWUP_FIRST_DAY", "MC_FOLLOWUP_FIRST": "P3_MCMED_FOLLOWUP_FIRST_DAY",
    "FOLLOWUP_LAST": "P3_FOLLOWUP_LAST_DAY", "OMR_MIN_READINGS": "P3_OMR_MIN_READINGS",
    "LANDMARK": "ESC_LANDMARK_MINUTES", "CI95": "P2_MODEL_CONFIDENCE_LEVEL", "MODEL_MAX_HOURS": "P2_MODEL_HOURS",
    "LMM_SPLINE_DF": "P2_MODEL_SPLINE_DF", "THINNING_REPS": "SPARSE_BP_REPETITIONS",
    "MATCH_LOOKBACK": "P5_MATCH_ANCHOR_LOOKBACK_MINUTES", "MATCH_TIME_CALIPER": "P5_MATCH_TIME_CALIPER_MINUTES",
    "MATCH_SBP_CALIPER": "P5_MATCH_SBP_CALIPER_MMHG", "MATCH_MAX_CONTROLS": "P5_MATCH_MAX_CONTROLS",
    "MATCH_FOLLOWUP": "P5_MATCH_FOLLOWUP_MINUTES", "MATCH_TIME_FIRST": "P5_MATCH_RESPONSE_WINDOWS",
    "MIN120": "P5_MATCH_RESPONSE_WINDOWS", "MATCH_BOOT_REPS": "P5_MATCH_BOOTSTRAPS",
    "OVERLAP_EVENT_MIN": "P5_RARE_TREATED_EVENT_COUNT", "FIG2_MAX_HOURS": "P6_FIG2_MAIN_MAX_HOURS",
    "SENS_INTENSITY_SBP": "P4_S11_HIGH_SBP", "SENS_INTENSITY_DBP": "P4_S11_HIGH_DBP",
    "OMR_STRICT_SBP": "P3_OMR_STRICT_SBP", "OMR_STRICT_DBP": "P3_OMR_STRICT_DBP",
    "OMR_STRICT_MIN": "P3_OMR_STRICT_MIN_READINGS", "PRIOR_EPISODES_DAYS": "PRIOR_EPISODES_LOOKBACK_DAYS",
    "SBP_INCREMENT": "P4_SBP_UNIT",
}


def manuscript_method_value(symbol, configured):
    if symbol in ("RETURN7", "RETURN30"):
        return str(configured // 24), "hours divided by 24"
    if symbol == "CI95":
        return f"{configured * 100:g}%", "proportion expressed as percent"
    if symbol == "MODEL_MAX_HOURS":
        return str(configured[1]), "upper bound of the window"
    if symbol in ("MATCH_TIME_FIRST", "MIN120"):
        target = configured[0] if symbol == "MATCH_TIME_FIRST" else configured[1]
        return str(target[0]), "target time in minutes"
    return str(configured), "configured value"


def methods_numbers():
    cfg = import_module("config")
    file = "code/config.py"
    code = (ROOT / file).read_text(encoding="utf-8-sig")
    lines = code.splitlines()
    locations = {}
    for node in ast.parse(code).body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                names = target.elts if isinstance(target, ast.Tuple) else (target,)
                locations.update((item.id, node.lineno) for item in names if isinstance(item, ast.Name))
    for symbol, variable in METHOD_NUMBERS.items():
        line = locations[variable]
        value, units = manuscript_method_value(symbol, getattr(cfg, variable))
        yield (f"Manuscript method {symbol} ({units})", value, None, None, None, None,
               file, line, variable, lines[line - 1].strip(), None, None, None, None)

    readme = "README_code.md" if (ROOT / "README_code.md").is_file() else "README.md"
    versions = {"MCMED_VERSION": "v1.0.1", "MIMIC_ED_VERSION": "v2.2", "MIMIC_HOSP_VERSION": "v3.1"}
    line, source = next((line, text) for line, text in enumerate((ROOT / readme).read_text(encoding="utf-8").splitlines(), 1)
                        if text.startswith("Obtain access to the confirmed study releases"))
    for symbol, value in versions.items():
        if value not in source:
            raise ValueError("A manuscript dataset version differs from the confirmed release.")
        yield (f"Manuscript method {symbol} (confirmed release)", value, None, None, None, None,
               readme, line, "release declaration", value, None, None, None, None)

    writing = "docs/05_写作要求.md"
    if (ROOT / writing).is_file():
        line, ethics = next((line, text) for line, text in enumerate((ROOT / writing).read_text(encoding="utf-8").splitlines(), 1)
                            if text.lstrip().startswith("- **Ethics**"))
        for symbol, value, anchor in (("ETHICS_ARTICLE", "32", "Article 32"), ("ETHICS_YEAR", "2023", "(2023)")):
            if anchor not in ethics:
                raise ValueError("The supplied ethics statement differs from its manuscript reference.")
            yield (f"Manuscript method {symbol} (supplied ethics wording)", value, None, None, None, None,
                   writing, line, "user-confirmed ethics text", anchor, None, None, None, None)

    reviewed_values = (
        ("PRIMARY_EXCLUDED_INCREASE", "50%", "review/G4_审核意见.md", "7 天止血失败风险增加 50%"),
        ("MATCH_BOUNDARY_MODELS", "7", "review/G5_审核意见.md", "7 个随机截距方差接近 0"),
    )
    for symbol, value, path, anchor in reviewed_values:
        if (ROOT / path).is_file():
            line = next(line for line, text in enumerate((ROOT / path).read_text(encoding="utf-8-sig").splitlines(), 1)
                        if anchor in text)
            yield (f"Manuscript method {symbol} (reviewed reporting value)", value, None, None, None, None,
                   path, line, "reviewed interpretation", anchor, None, None, None, None)


def as_number(value):
    try:
        number = Decimal(value.replace(",", "").replace("−", "-"))
    except InvalidOperation:
        return None
    return number if number.is_finite() else None


@lru_cache(maxsize=None)
def csv_rows(relative):
    path = (ROOT / relative).resolve()
    if not path.is_relative_to(ROOT / "outputs") or path.suffix.lower() != ".csv":
        raise ValueError("The ledger only accepts aggregate output CSV paths.")
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        reader.fieldnames
        rows = {}
        start = reader.line_num + 1
        for row in reader:
            rows[start] = row
            start = reader.line_num + 1
    return rows


def describe(record, table):
    labels = [f"{field}: {value}" for field, value in record.items()
              if field.lower().replace(" ", "_") in DESCRIPTORS and value]
    return f"{table} — {'; '.join(labels)}" if labels else table


def source_interval(source_file, source_line, column, source_row):
    for lo, hi in INTERVAL_COLUMNS.get(column, ()):
        if as_number(source_row.get(lo, "")) is not None and as_number(source_row.get(hi, "")) is not None:
            return source_row[lo], source_row[hi], f"{source_file}:{source_line}:{lo}", f"{source_file}:{source_line}:{hi}"
    return (None,) * 4


def bounds(record, cell_group, source_row):
    if record["component"].lower() in ("estimate", "percent", "rr", "estimate_mmhg"):
        low = next((item for item in cell_group if item["component"].lower() == "95% ci lower"), None)
        high = next((item for item in cell_group if item["component"].lower() == "95% ci upper"), None)
        if low and high:
            references = [f"{part['source_file']}:{part['source_line']}:{part['source_column']}" for part in (low, high)]
            return low["source_value"], high["source_value"], *references
    return source_interval(record["source_file"], record["source_line"], record["source_column"], source_row)


def compatible_display(number, display):
    for token in NUMERIC_TOKEN.findall(display):
        shown = as_number(token)
        if shown is None:
            continue
        tolerance = Decimal(5).scaleb(shown.as_tuple().exponent - 1)
        if abs(shown - number) <= tolerance or ("%" in display and abs(shown - 100 * number) <= tolerance):
            return True
    return False


def cite_final_tables(provenance):
    groups = defaultdict(list)
    for row in provenance.values():
        if Path(row["final_file"]).name.startswith("supp_dictionary_") or not row["source_file"].endswith(".csv"):
            continue
        groups[row["final_file"], int(row["final_line"]), row["final_column"]].append(row)
    seen = set()
    records = []
    for (final_file, final_line, final_column), locators in groups.items():
        final_row = csv_rows(final_file)[final_line]
        table = Path(final_file).stem
        display = final_row[final_column]
        description = describe(final_row, table)
        for link in locators:
            source_file = link["source_file"]
            source_line = int(link["source_line"])
            source_column = link["source_column"]
            source_row = csv_rows(source_file)[source_line]
            value = source_row[source_column]
            if value != link["source_value"] or display != link["display_value"]:
                raise ValueError("Final-cell provenance differs from a referenced CSV cell.")
            if as_number(value) is not None:
                if table.startswith("table") and not compatible_display(as_number(value), display):
                    raise ValueError(f"Numeric display mismatch at {table} line {final_line}, column {final_column}.")
                values = [(link["component"], value, bounds(link, locators, source_row))]
            elif table == "table4_hemostatic_outcomes" and link["component"] == "Modeled events/visits by group":
                values = [(f"{group} {field}", number, (None,) * 4)
                          for group, events, visits in GROUP_TOKEN.findall(value)
                          for field, number in (("events", events), ("visits", visits))]
            else:
                continue
            for component, estimate, (lo, hi, lo_ref, hi_ref) in values:
                records.append((f"{description} — {final_column} / {component}", estimate, lo, hi, lo_ref, hi_ref,
                                source_file, source_line, source_column, value, final_file, final_line, final_column, display))
            seen.add((source_file, source_line, source_column))
    return records, seen


def additional_sources(seen):
    for filename in EXTRA_SOURCES:
        path = TABLES / filename
        if not path.is_file():
            continue
        relative = path.relative_to(ROOT).as_posix()
        for line, row in csv_rows(relative).items():
            description = describe(row, path.stem)
            for column, value in row.items():
                if as_number(value) is None or (relative, line, column) in seen:
                    continue
                low, high, low_ref, high_ref = source_interval(relative, line, column, row)
                yield (f"{description} — {column}", value, low, high, low_ref, high_ref,
                       relative, line, column, value, None, None, None, None)


def main():
    stage = "environment"
    try:
        import_module("02_p1_setup").verify_project_runtime()
        writer = import_module("01_p0_audit").write_csv
        stage = "main and supplemental cell cross-check"
        provenance = csv_rows(PROVENANCE.relative_to(ROOT).as_posix())
        records, seen = cite_final_tables(provenance)
        stage = "figure and text-support inputs"
        records.extend(additional_sources(seen))
        stage = "manuscript method source locators"
        records.extend(methods_numbers())
        stage = "aggregate ledger export"
        header = ("key_number_id", "description", "value", "ci95_low", "ci95_high", "ci95_low_reference",
                  "ci95_high_reference", "source_file", "source_line", "source_column", "source_cell_value",
                  "manuscript_file", "manuscript_line", "manuscript_column", "display_cell")
        writer(ROOT / "outputs" / "key_numbers.csv", header,
               [(f"K{index:05d}", *record) for index, record in enumerate(records, start=1)])
    except Exception as error:
        # All mismatches use only aggregate table names and column labels.
        if isinstance(error, ValueError) and str(error).startswith("Numeric display mismatch at "):
            print(str(error), file=sys.stderr)
        print(f"P6 ledger stopped at {stage} ({type(error).__name__}).", file=sys.stderr)
        return 1
    print(f"P6 ledger verified {len(records)} numeric entries.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
