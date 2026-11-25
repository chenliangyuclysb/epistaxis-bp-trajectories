"""Format the audited P1–P6 aggregate results and retain cell-level provenance."""

from collections import Counter
import csv
from dataclasses import dataclass
from importlib import import_module
from io import StringIO
from math import isfinite
from pathlib import Path
import sys

import duckdb


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "outputs" / "tables"
FINAL = TABLES / "final"
DICTIONARY = ROOT / "outputs" / "dictionary_final"
SITES = (("mcmed", "MC-MED"), ("mimic", "MIMIC-IV-ED"))
NA = "—"

# SAP S9, G5 6.1c: use literal upper-case prefixes for reported race.
RACE_PREFIXES = (
    ("White", ("WHITE",)), ("Black", ("BLACK",)),
    ("Hispanic or Latino", ("HISPANIC", "LATINO")), ("Asian", ("ASIAN",)),
)
OTHER_RACE = "Other or Unknown"
RACE_LEVELS = tuple(label for label, _ in RACE_PREFIXES) + (OTHER_RACE,)
HISPANIC = "Hispanic or Latino"
COMORBIDITIES = ("ckd", "liver_disease", "hht")
COMORBIDITY_SOURCE = TABLES / "p7_table1_comorbidity.csv"


@dataclass(frozen=True)
class SourceRow:
    path: str
    line: int
    values: dict

    def __getitem__(self, field):
        return self.values[field]


@dataclass(frozen=True)
class Cell:
    value: str
    references: tuple


def number(value):
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if isfinite(result) else None


def numeric_cell(row, field, text=None, component=None):
    raw = row[field]
    if not raw or str(raw).lower() in ("nan", "na", "none"):
        return NA
    return Cell(str(raw) if text is None else text, ((component or field, row, field),))


def formatted(row, field, digits=0, multiplier=1, component=None):
    value = number(row[field])
    if value is None:
        return NA
    return numeric_cell(row, field, f"{value * multiplier:,.{digits}f}", component)


def compound(parts, render):
    if any(not isinstance(part, Cell) for part in parts):
        return NA
    refs = tuple(ref for part in parts for ref in part.references)
    return Cell(render(*(part.value for part in parts)), refs)


def interval(row, estimate, lo, hi, digits=2, multiplier=1, suffix=""):
    parts = (formatted(row, estimate, digits, multiplier, "Estimate"),
             formatted(row, lo, digits, multiplier, "95% CI lower"),
             formatted(row, hi, digits, multiplier, "95% CI upper"))
    return compound(parts, lambda point, low, high: f"{point}{suffix} ({low}–{high}{suffix})")


def iq_range(row, median, low, high, digits=1):
    parts = (formatted(row, median, digits, component="Median"),
             formatted(row, low, digits, component="First quartile"),
             formatted(row, high, digits, component="Third quartile"))
    return compound(parts, lambda point, q1, q3: f"{point} ({q1}–{q3})")


def ratio_counts(row, numerator, denominator):
    parts = (formatted(row, numerator, component="Event or exposed count"),
             formatted(row, denominator, component="Denominator"))
    return compound(parts, lambda count, total: f"{count}/{total}")


def count_percent(row, count, percent):
    return compound((formatted(row, count, component="Count"), formatted(row, percent, 1, component="Percent")),
                    lambda n, pct: f"{n} ({pct}%)")


def input_rows(path):
    path = Path(path)
    relative = path.relative_to(ROOT).as_posix()
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = tuple(reader.fieldnames or ())
        previous = reader.line_num
        rows = []
        for values in reader:
            if None in values or any(value is None for value in values.values()):
                raise ValueError("An audited aggregate CSV has a row with an unexpected width.")
            rows.append(SourceRow(relative, previous + 1, values))
            previous = reader.line_num
    if not fields:
        raise ValueError("An audited aggregate CSV has no header.")
    return fields, rows


def csv_line(values):
    buffer = StringIO(newline="")
    csv.writer(buffer, lineterminator="\n").writerow(values)
    return buffer.getvalue()


class TableWriter:
    def __init__(self):
        self.sources = []
        self.index = []

    def emit(self, filename, fields, rows, title, note=""):
        if not rows:
            raise ValueError(f"P6 final table {filename} has no aggregate records.")
        destination = FINAL / filename
        final_path = destination.relative_to(ROOT).as_posix()
        pending = destination.with_suffix(".csv.tmp")
        line = 2
        try:
            with pending.open("w", encoding="utf-8-sig", newline="") as handle:
                handle.write(csv_line(fields))
                for row in rows:
                    if len(row) != len(fields):
                        raise ValueError(f"P6 final table {filename} has a row with an unexpected width.")
                    values = []
                    for column, item in zip(fields, row):
                        if isinstance(item, Cell):
                            values.append(item.value)
                            for part, source, source_column in item.references:
                                self.sources.append((final_path, line, column, part, source.path, source.line,
                                                     source_column, source[source_column], item.value))
                        else:
                            values.append(item)
                    record = csv_line(values)
                    handle.write(record)
                    line += record.count("\n")
            pending.replace(destination)
        finally:
            pending.unlink(missing_ok=True)
        self.index.append((final_path, title, len(rows), note))

    def copy_table(self, path, filename, title, note="", select=None):
        fields, rows = input_rows(path)
        if select is not None:
            rows = [row for row in rows if select(row)]
        self.emit(filename, fields, [tuple(numeric_cell(row, key) if number(row[key]) is not None else row[key]
                                           for key in fields) for row in rows], title, note)

    def finish(self):
        write_csv = import_module("01_p0_audit").write_csv
        write_csv(FINAL / "p6_table_sources.csv", ("final_file", "final_line", "final_column", "component",
            "source_file", "source_line", "source_column", "source_value", "display_value"), self.sources)
        write_csv(FINAL / "p6_table_index.csv", ("final_file", "table_title", "data_rows", "table_note"), self.index)


def grouped_race(center, race, ethnicity):
    race_label = (race or "").strip().upper()
    ethnicity_label = (ethnicity or "").strip().upper()
    if center == "MC-MED" and ethnicity_label.startswith(("HISPANIC", "LATINO")):
        return HISPANIC
    return next((label for label, prefixes in RACE_PREFIXES if race_label.startswith(prefixes)), OTHER_RACE)


def write_grouped_race():
    totals, detail = [], []
    with duckdb.connect() as con:
        for short, site in SITES:
            grouped = con.execute("""
                SELECT race, ethnicity, count(*) AS visits FROM read_parquet(?) WHERE is_c1
                GROUP BY race, ethnicity ORDER BY race, ethnicity
            """, [str(ROOT / "derived" / f"cohort_{short}.parquet")]).fetchall()
            denominator = sum(count for _, _, count in grouped)
            counts = Counter()
            for race, ethnicity, visits in grouped:
                label = grouped_race(site, race, ethnicity)
                counts[label] += visits
                detail.append((site, race or "Missing", ethnicity or "Missing or unavailable", label, visits))
            if denominator <= 0 or sum(counts.values()) != denominator:
                raise ValueError("The P6 grouped race denominator did not reconcile.")
            totals.extend((site, label, denominator, counts[label], 100 * counts[label] / denominator) for label in RACE_LEVELS)
    write_csv = import_module("01_p0_audit").write_csv
    write_csv(TABLES / "p6_table1_grouped_race.csv", ("site", "level", "denominator", "count", "pct"), totals)
    write_csv(TABLES / "p6_table1_detailed_race_mapping.csv", ("site", "source_race", "source_ethnicity", "level", "visits"), detail)


def write_g6_comorbidities():
    # G6 3.1 and SAP S9: use the G1-updated analysis flags on the fixed C1 population.
    rows = []
    fields = ", ".join(f"count({name}), count(*) FILTER (WHERE {name})" for name in COMORBIDITIES)
    with duckdb.connect() as con:
        for short, site in SITES:
            counts = con.execute(f"""
                SELECT count(*), {fields}
                FROM read_parquet(?) WHERE is_c1
            """, [str(ROOT / "derived" / f"cohort_{short}.parquet")]).fetchone()
            total, *subgroups = counts
            for name, (observed, yes) in zip(COMORBIDITIES, zip(subgroups[::2], subgroups[1::2])):
                if not total or observed != total:
                    raise ValueError("The G1-updated comorbidity flags must be known for all C1 visits.")
                rows.append((site, name, "Yes", observed, yes, 100 * yes / observed))
    write_csv = import_module("01_p0_audit").write_csv
    write_csv(COMORBIDITY_SOURCE, ("site", "characteristic", "level", "denominator", "count", "pct"), rows)


BASELINE_LABELS = {
    "visits": "Index visits", "patients": "Unique patients", "age_years": "Age, years", "sex": "Sex",
    "race_group": "Race or ethnicity", "arrival_mode": "Arrival mode", "triage_acuity_group": "Triage acuity",
    "triage_heart_rate": "Triage heart rate, beats/min", "triage_pain_score": "Pain score, 0–10",
    "prior_htn": "Known hypertension", "has_medication_list": "Medication list documented",
    "anticoagulant_use": "Anticoagulant use", "antiplatelet_use": "Antiplatelet use",
    "any_antithrombotic_use": "Any antithrombotic use", "ckd": "Chronic kidney disease",
    "liver_disease": "Liver disease", "hht": "Hereditary hemorrhagic telangiectasia",
    "prior_epistaxis_365d": "Prior epistaxis ED visit within 365 days", "hemoglobin": "Hemoglobin, g/dL",
    "platelet_count": "Platelet count, ×10³/μL", "inr": "International normalized ratio",
    "hemostasis_level2": "Level 2 hemostatic treatment",
}


def table1(writer):
    _, historical = input_rows(TABLES / "p1_table1.csv")
    _, grouped = input_rows(TABLES / "p6_table1_grouped_race.csv")
    _, revised = input_rows(COMORBIDITY_SOURCE)
    by_site = {(row["site"], row["characteristic"], row["level"]): row for row in historical + revised}
    by_site.update({(row["site"], "race_group", row["level"]): row for row in grouped})
    labels = [site for _, site in SITES]
    for site in labels:
        eligible = by_site[(site, "visits", "")]["count"]
        if any(row["denominator"] != eligible for row in grouped if row["site"] == site):
            raise ValueError("P6 grouped race and audited Table 1 C1 denominators differ.")
        if any(row["denominator"] != eligible for row in revised if row["site"] == site):
            raise ValueError("The G6 comorbidity and audited Table 1 C1 denominators differ.")

    rows = []
    for variable, label in BASELINE_LABELS.items():
        levels = RACE_LEVELS if variable == "race_group" else sorted({key[2] for key in by_site if key[1] == variable})
        if len(levels) > 1 and "" in levels:
            levels.remove("")
        for level in levels:
            values, denominators = [], []
            for site in labels:
                row = by_site.get((site, variable, level))
                if row is None:
                    values.append(NA)
                    denominators.append(NA)
                    continue
                denominators.append(formatted(row, "denominator"))
                if number(row.values.get("median")) is not None:
                    values.append(iq_range(row, "median", "q1", "q3", 2 if variable == "inr" else 1))
                elif number(row["pct"]) is not None:
                    values.append(count_percent(row, "count", "pct"))
                else:
                    values.append(formatted(row, "count"))
            rows.append((label, level, *values, *denominators))
    writer.emit("table1_baseline.csv", ("Characteristic", "Level", *labels,
        "MC-MED observed visits", "MIMIC-IV-ED observed visits"), rows,
        "Table 1. Baseline characteristics of C1 by center",
        "Values are median (IQR) or count (%). Grouped race or ethnicity includes unknown values in the full C1 denominator. Chronic kidney and liver disease can include documentation through the index visit. Hereditary hemorrhagic telangiectasia has no time restriction. An em dash indicates unavailable data or a level absent from the site's source table.")


PHENOTYPE_LABELS = {
    "arrival_elevated": "Elevated arrival BP", "arrival_severe": "Severely elevated arrival BP",
    "last_normal": "Normal final BP", "N": "N, normal arrival", "T": "T, transient elevation",
    "P": "P, persistent elevation", "P-never": "P-never, no normal reading", "P-rebound": "P-rebound, normal then elevated",
    "n_bp": "BP readings per visit", "twa_sbp": "Time-weighted SBP, mmHg", "hypertension_load_sbp": "SBP burden, mmHg",
    "sbp_arv": "Average real variability in SBP, mmHg", "peak_sbp": "Peak SBP, mmHg", "ever_severe": "Any severely elevated BP",
}


def table2(writer):
    _, source = input_rows(TABLES / "p2_table2.csv")
    rows = []
    for row in source:
        continuous = number(row["median"]) is not None
        values = iq_range(row, "median", "q1", "q3") if continuous else count_percent(row, "n", "pct")
        limits = NA if continuous else interval(row, "pct", "wilson95_lo_pct", "wilson95_hi_pct", 1, suffix="%")
        sample = formatted(row, "n") if continuous else formatted(row, "denominator")
        rows.append((row["site"], row["phenotype"], PHENOTYPE_LABELS[row["metric"]], sample, values, limits))
    writer.emit("table2_phenotypes.csv", ("Center", "Population", "Measure", "Observed visits", "Summary", "Percent (Wilson 95% CI)"),
        rows, "Table 2. BP phenotypes and MC-MED dense-reading metrics",
        "Continuous measures are median (IQR). Proportions are count (%) with Wilson 95% CI. The P subtype denominator is the number of P visits.")


def table3(writer):
    _, source = input_rows(TABLES / "p3_table3.csv")
    rows = []
    for row in sorted(source, key=lambda item: item["site"] != "MIMIC-IV-ED"):
        measure = row["estimand"]
        if measure == "Negative predictive value audit":
            continue
        risk = measure != "Adjusted risk ratio"
        rows.append((row["site"], "Main evidence" if row["site"] == "MIMIC-IV-ED" else "Supporting evidence",
            measure, row["group_or_contrast"], formatted(row, "eligible_patients"), formatted(row, "followup_patients"),
            ratio_counts(row, "numerator", "denominator"), ratio_counts(row, "reference_numerator", "reference_denominator"),
            interval(row, "estimate", "ci95_low", "ci95_high", 1 if risk else 2, 100 if risk else 1, "%" if risk else ""), row["status"]))
    writer.emit("table3_followup_hypertension.csv", ("Center", "Evidence role", "Measure", "Group or contrast", "C3 patients",
        "Patients with follow-up", "Numerator / denominator", "Reference numerator / denominator", "Estimate (95% CI)", "Status"), rows,
        "Table 3. Subsequent hypertension evidence in patients without known hypertension",
        "Each patient contributes their first eligible C3 visit. MIMIC provides the primary evidence from outpatient BP. MC-MED provides supporting evidence from PMH coding or new antihypertensive treatment. The negative predictive value uses T as the negative group.")


def table4(writer):
    _, source = input_rows(TABLES / "p4_table4.csv")
    _, contextual = input_rows(TABLES / "p4_secondary_estimates.csv")
    exposure_types = {"Arrival SBP per 10 mmHg": "arrival_sbp_per_10", "N60/T60/P60 phenotype": "landmark_phenotype",
                      "N/T/P phenotype": "bp_phenotype"}
    rows = []
    for row in source:
        if row["analysis"] == "Any admission landmark sensitivity":
            continue
        context = NA
        clinical = [item for item in contextual if item["site"] == row["site"] and item["analysis"] == row["analysis"]
                    and item["exposure"] == exposure_types[row["exposure"]]]
        if clinical:
            context = ratio_counts(clinical[0], "raw_risk_set_events", "raw_risk_set_visits")
        rows.append((row["endpoint"], row["site"], row["analysis"], row["exposure"], row["adjustment"], row["contrast"],
            context, formatted(row, "visits"), formatted(row, "patients"), formatted(row, "events"),
            numeric_cell(row, "modeled_events_per_visits_by_group", component="Modeled events/visits by group"),
            interval(row, "rr", "ci95_low", "ci95_high"), row["status"]))
    writer.emit("table4_hemostatic_outcomes.csv", ("Outcome", "Center", "Analysis", "Exposure", "Adjustment", "Contrast",
        "Clinical risk set events / visits", "Modeled visits", "Modeled patients", "Modeled events", "Modeled group counts",
        "Risk ratio (95% CI)", "Status"), rows, "Table 4. Early hemostatic failure and secondary outcomes",
        "The primary comparison is pooled P vs T for 7-day return. Landmark phenotype requires a reading within (0, 60] minutes. The inpatient model uses completed C1 visits within the linked-admission risk set.")


def table5(writer):
    _, frequencies = input_rows(TABLES / "p5_exposure_frequency.csv")
    _, timing = input_rows(TABLES / "p5_exposure_timing.csv")
    _, no_admin = input_rows(TABLES / "p5_exposure_mcmed_no_administration.csv")
    _, response = input_rows(TABLES / "p5_matched_response.csv")
    _, returns = input_rows(TABLES / "p5_treatment_return_7d.csv")
    rows = []
    for row in frequencies:
        if row["exposure_definition"] == "Class A":
            rows.append(("A-class treatment", row["site"], row["population"], "ED administration", ratio_counts(row, "treated_visits", "eligible_visits"),
                         interval(row, "treated_percent", "wilson_ci95_low", "wilson_ci95_high", 1, suffix="%"), "Wilson 95% CI", "", ""))
    for row in timing:
        if row["exposure_definition"] != "Class A":
            continue
        site, population = row["site"], row["population"]
        rows.extend((
            ("Timing", site, population, "Time to first dose, min", formatted(row, "treated_visits"),
             iq_range(row, "first_dose_minutes_median", "first_dose_minutes_q25", "first_dose_minutes_q75"), "Median (IQR)", "", ""),
            ("Timing", site, population, "Last SBP before first dose, mmHg", formatted(row, "strict_pre_dose_bp_visits"),
             iq_range(row, "last_pre_dose_sbp_median", "last_pre_dose_sbp_q25", "last_pre_dose_sbp_q75"), "Median (IQR)", "", ""),
            ("Timing", site, population, "Severe BP elevation before first dose", ratio_counts(row, "severe_before_dose_visits", "strict_pre_dose_bp_visits"),
             interval(row, "severe_before_dose_percent_observed", "severe_wilson_ci95_low", "severe_wilson_ci95_high", 1, suffix="%"), "Wilson 95% CI", "", ""),
        ))
    for row in no_admin:
        if row["exposure_definition"] == "Class A":
            rows.append(("Order audit", row["site"], row["population"], "A-class ED orders with no documented administration",
                formatted(row, "orders_without_documented_administration"),
                numeric_cell(row, "visits_with_orders_without_documented_administration", f'{row["visits_with_orders_without_documented_administration"]} visits'), "", "", ""))
    for row in response:
        if row["scenario"] == "primary":
            row_time = numeric_cell(row, "minutes", f'+{row["minutes"]} min')
            side = compound((formatted(row, "treated_mean_sbp_change", 1, component="Treated mean change"),
                             formatted(row, "control_mean_sbp_change", 1, component="Control mean change")),
                            lambda treated, control: f"{treated} vs {control} mmHg")
            rows.append(("Matched response", row["site"], row_time, "SBP change difference, mmHg", formatted(row, "matched_sets"),
                interval(row, "treated_minus_control_change", "ci95_low", "ci95_high", 1), "Matched-set bootstrap 95% CI", side, ""))
    for row in returns:
        if row["treatment_definition"] == "A":
            rows.append(("Crude 7-day return", row["site"], "C2", "A-class treatment vs no A-class treatment",
                ratio_counts(row, "treated_events", "treated_visits"), interval(row, "crude_rr", "katz_ci95_low", "katz_ci95_high"),
                "Katz approximate 95% CI", ratio_counts(row, "untreated_events", "untreated_visits"), formatted(row, "fisher_two_sided_p", 3)))
    writer.emit("table5_treatment.csv", ("Panel", "Center", "Population or time", "Measure", "Count or denominator",
        "Estimate with interval (Katz approximate 95% CI for crude RR)", "Interval method", "Comparator or group means", "Fisher exact two-sided p"), rows,
        "Table 5. ED antihypertensive treatment and matched BP response",
        "A-class drug analyses are descriptive. Timing and pre-dose SBP use median (IQR). Treatment and severe BP proportions use Wilson 95% CI. Matched differences use matched-set bootstrap 95% CI. Crude return RRs use Katz approximate 95% CI and a separate two-sided Fisher exact p value.")


def table_s11(writer):
    rows = []
    _, main = input_rows(TABLES / "p4_sensitivity_estimates.csv")
    for row in main:
        rows.append((row["scenario"], "Aim 3", row["site"], row["analysis"], row["outcome"], row["contrast"],
            formatted(row, "model_visits"), formatted(row, "model_events"), interval(row, "rr", "ci95_low", "ci95_high"), row["status"]))

    _, outcomes = input_rows(TABLES / "p4_table4.csv")
    for row in outcomes:
        if row["endpoint"] == "30-day return":
            rows.append(("S11-8", "Aim 3", row["site"], "30-day outcome", row["endpoint"], row["contrast"],
                formatted(row, "visits"), formatted(row, "events"), interval(row, "rr", "ci95_low", "ci95_high"), row["status"]))

    _, latent = input_rows(TABLES / "p2_mcmed_lcmm_models.csv")
    for row in latent:
        category = numeric_cell(row, "k", f'{row["k"]} classes')
        status = "Selected" if row["selected"] == "TRUE" else "Eligible" if row["eligible"] == "TRUE" else "Not eligible"
        for label, field, multiplier, digits in (
            ("BIC", "bic", 1, 2), ("Minimum assigned class, %", "min_class_proportion", 100, 2),
            ("Minimum class mean assigned posterior", "min_mean_assigned_posterior", 1, 3),
        ):
            rows.append(("S11-10", "Aim 1", "MC-MED", "Latent class mixed model", label, category,
                formatted(row, "visit_n"), NA, formatted(row, field, digits, multiplier), status))

    _, evalue = input_rows(TABLES / "p4_main_evalue.csv")
    for row in evalue:
        for label, field in (("Point estimate E-value", "point_e_value"), ("Closest CI limit E-value", "ci_limit_e_value")):
            rows.append(("S11-11", "Aim 3", row["site"], "E-value", label, "P vs T", NA, NA,
                         formatted(row, field, 2), row["status"]))

    _, dedup = input_rows(TABLES / "p2_dense_dedup_sensitivity.csv")
    for row in dedup:
        window = numeric_cell(row, "window_minutes", f'±{row["window_minutes"]} min')
        summary = iq_range(row, "median", "q1", "q3") if number(row["median"]) is not None else formatted(row, "true_count")
        label = PHENOTYPE_LABELS.get(row["metric"], {
            "pretreatment_first_normal_minutes": "Pre-treatment first normal BP, min",
            "recorded_first_normal_minutes": "Recorded first normal BP, min", "would_meet_c1": "Meets C1 under the window",
        }.get(row["metric"], row["metric"]))
        rows.append(("S11-14", "Aim 1", "MC-MED", "Dense-reading deduplication", label, window,
            formatted(row, "n"), NA, summary, "Descriptive"))
    _, recovery = input_rows(TABLES / "p2_dense_km_sensitivity.csv")
    for row in recovery:
        group = "Known hypertension" if row["prior_htn"].lower() == "true" else "No known hypertension"
        rows.append(("S11-14", "Aim 1", row["site"], row["analysis"], "Median recovery, min", group,
            formatted(row, "km_visits"), formatted(row, "recovered_visits"), formatted(row, "median_recovery_minutes", 1),
            "Descriptive"))

    _, first_aim = input_rows(TABLES / "p2_aim1_prevalence_ratio.csv")
    aim1_s11 = {"Exclude acute class A treatment (S10)": "S11-6", "Medication list available (S11-16)": "S11-16",
                "Ignore prior ED hypertension diagnosis (S11-17)": "S11-17"}
    for row in first_aim:
        if row["analysis"] in aim1_s11:
            rows.append((aim1_s11[row["analysis"]], "Aim 1", row["site"], row["analysis"], "Prior hypertension prevalence ratio",
                row["contrast"], formatted(row, "patient_clusters"), NA, interval(row, "prevalence_ratio", "ci95_low", "ci95_high"),
                row["status"]))
    _, first_proportions = input_rows(TABLES / "p2_aim1_proportions.csv")
    for row in first_proportions:
        if row["analysis"] in aim1_s11:
            rows.append((aim1_s11[row["analysis"]], "Aim 1", row["site"], row["analysis"], row["measure"], "",
                formatted(row, "denominator"), formatted(row, "numerator"),
                interval(row, "percent", "wilson95_low_pct", "wilson95_high_pct", 1, suffix="%"), "Descriptive"))

    _, hypertension = input_rows(TABLES / "p3_sensitivity.csv")
    aim2_s11 = {"Count same-day HTN diagnosis": "S11-13", "Medication list documented": "S11-16",
                "Exclude prior ED HTN diagnosis": "S11-17", "At least 2 readings, threshold 130/80": "S6.2",
                "Follow-up starts day 7": "S6.2", "Follow-up IPW": "S10-IPW"}
    for row in hypertension:
        risk = row["estimand"] != "Adjusted risk ratio"
        rows.append((aim2_s11[row["analysis"]], "Aim 2", row["site"], row["analysis"], row["estimand"], row["group_or_contrast"],
            formatted(row, "denominator"), formatted(row, "numerator"),
            interval(row, "estimate", "ci95_low", "ci95_high", 1 if risk else 2, 100 if risk else 1, "%" if risk else ""),
            row["status"]))

    _, frequency = input_rows(TABLES / "p5_exposure_frequency.csv")
    for row in frequency:
        if row["exposure_definition"] in ("Class B", "Class B + alpha-1"):
            rows.append(("S11-15", "Aim 4", row["site"], row["exposure_definition"], "ED administration",
                row["population"], formatted(row, "eligible_visits"), formatted(row, "treated_visits"),
                interval(row, "treated_percent", "wilson_ci95_low", "wilson_ci95_high", 1, suffix="%"), "Descriptive"))
    _, timing = input_rows(TABLES / "p5_exposure_timing.csv")
    for row in timing:
        if row["exposure_definition"] in ("Class B", "Class B + alpha-1"):
            for label, count_field, columns in (
                ("Time to first dose, min", "treated_visits", ("first_dose_minutes_median", "first_dose_minutes_q25", "first_dose_minutes_q75")),
                ("Last pre-dose SBP, mmHg", "strict_pre_dose_bp_visits", ("last_pre_dose_sbp_median", "last_pre_dose_sbp_q25", "last_pre_dose_sbp_q75")),
            ):
                rows.append(("S11-15", "Aim 4", row["site"], row["exposure_definition"], label, row["population"],
                    formatted(row, count_field), NA, iq_range(row, *columns), "Descriptive"))
            rows.append(("S11-15", "Aim 4", row["site"], row["exposure_definition"], "Severe BP before first dose", row["population"],
                formatted(row, "strict_pre_dose_bp_visits"), formatted(row, "severe_before_dose_visits"),
                interval(row, "severe_before_dose_percent_observed", "severe_wilson_ci95_low", "severe_wilson_ci95_high", 1, suffix="%"),
                "Descriptive"))
    _, treatment = input_rows(TABLES / "p5_treatment_return_7d.csv")
    for row in treatment:
        if row["treatment_definition"] in ("B", "B + alpha-1"):
            rows.append(("S11-15", "Aim 4", row["site"], f'Class {row["treatment_definition"]}', "Crude 7-day return RR, Katz CI",
                "Treated vs untreated", formatted(row, "visits"), formatted(row, "treated_events"),
                interval(row, "crude_rr", "katz_ci95_low", "katz_ci95_high"), row["status"]))
    _, overlap = input_rows(TABLES / "p6_bclass_overlap_estimate.csv")
    for row in overlap:
        rows.append(("S11-15", "Aim 4", "Pooled", "B-class overlap weighting", "Exploratory 7-day return RR", "Treated vs untreated",
            formatted(row, "visits"), formatted(row, "treated_events"), interval(row, "estimate", "ci95_low", "ci95_high"), row["status"]))

    assigned = {row[0] for row in rows if str(row[0]).startswith("S11-")}
    if assigned != {f"S11-{index}" for index in range(1, 18)}:
        raise ValueError("The combined P6 sensitivity table does not cover all S11 items.")
    rows.sort(key=lambda row: int(row[0].split("-")[1]) if row[0].startswith("S11-") else 99)
    writer.emit("supp_s11_sensitivity_all.csv", ("SAP item", "Aim", "Center", "Analysis", "Measure", "Group or contrast",
        "Sample", "Events or exposed count", "Estimate or summary", "Status"), rows,
        "Supplementary table. Prespecified sensitivity results",
        "All 17 S11 items are included. Other prespecified S6.2 and S10 analyses are appended. Depending on the row, the summary is an estimate (95% CI), median (IQR), or a scalar statistic. Each exact value is traced to the source row and column.")


def deviation_table(writer):
    source = ROOT / "review" / "deviations.md"
    if not source.is_file():
        # The code-only release can run without its private project review folder.
        return
    rows = []
    section, header = "General", ()
    relative = source.relative_to(ROOT).as_posix()
    for line_no, line in enumerate(source.read_text(encoding="utf-8-sig").splitlines(), start=1):
        if line.startswith("## "):
            section = line.removeprefix("## ")
        if not line.startswith("|"):
            continue
        tokens = tuple(part.strip() for part in line.strip().strip("|").split("|"))
        if tokens[0] == "日期":
            header = tokens
        elif tokens[0].replace("-", "").strip() == "":
            continue
        elif len(tokens) == len(header):
            original = SourceRow(relative, line_no, dict(zip(header, tokens)))
            rows.append((section, *(numeric_cell(original, field, component=field) for field in header), *([""] if len(header) == 3 else [])))
        else:
            raise ValueError("An entry in the supplied deviation log has an unexpected column count.")
    writer.emit("supp_protocol_deviations.csv", ("Section", "Date", "Basis", "Protocol detail", "Rationale"), rows,
        "Supplementary table. Deviations and implementation clarifications",
        "The approved project log is transcribed by original line. Entries with one combined detail and rationale field appear in Protocol detail. The protocol log's original language and wording are retained.")


def supplementary_tables(writer):
    dictionary_titles = {
        "mcmed_orders_mapping": "MC-MED medication orders and hemostatic escalation",
        "mimic_pyxis_mapping": "MIMIC ED medication dispensing",
        "mcmed_meds_mapping": "MC-MED outpatient treatment, by class and subclass",
        "mimic_medrecon_mapping": "MIMIC reconciled medications, by class",
        "mcmed_labs_mapping": "MC-MED laboratory components",
        "mcmed_disposition_mapping": "MC-MED disposition",
        "mcmed_reversal_g1_mapping": "MC-MED anticoagulant reversal review at G1",
    }
    for stem, title in dictionary_titles.items():
        note = "G2 subsequently approved andexanet. No matching order occurred in P4." if stem == "mcmed_reversal_g1_mapping" else ""
        writer.copy_table(DICTIONARY / f"{stem}.csv", f"supp_dictionary_{stem}.csv", f"Dictionary. {title}", note)

    supplements = (
        ("p1_inpatient_ent_candidate_codes", "dictionary_mimic_inpatient_procedures", "Accepted inpatient nasal hemostasis codes"),
        ("p1_table1", "baseline_detailed", "Original detailed P1 C1 baseline characteristics"),
        ("p6_table1_detailed_race_mapping", "baseline_race_mapping", "Source race and ethnicity categories and final groups"),
        ("p1_flow_counts", "flow", "Index cohort flow by stage and center"),
        ("p1_subcohort_exclusions", "flow_subcohorts", "C1 and C2 subcohort exclusions"),
        ("p3_patient_flow", "flow_c3", "Aim 2 C3 selection and first qualifying index"),
        ("p2_htn_sources", "hypertension_sources", "Recorded prior hypertension by source"),
        ("p2_sparse_summary", "sparse_sampling", "Monte Carlo sparse sampling and classification agreement"),
        ("p2_sparse_confusion", "sparse_sampling_confusion", "Sparse sampling confusion matrix"),
        ("p2_mcmed_lcmm_models", "lcmm_models", "Latent class model selection"),
        ("p2_mcmed_lcmm_classes", "lcmm_classes", "Latent class size and posterior probability"),
        ("p6_lmm_comparison", "lmm_satterthwaite", "Original and Satterthwaite linear mixed-model estimates"),
        ("p6_lmm_diagnostics", "lmm_diagnostics", "Linear mixed-model sample and diagnostic comparison"),
        ("p3_followup_availability", "aim2_followup_availability", "Observed follow-up by center and phenotype"),
        ("p3_followup_baseline", "aim2_followup_baseline", "Baseline characteristics with and without follow-up"),
        ("p3_table3", "aim2_primary_detail", "Full primary Aim 2 results and broader NPV reference coding"),
        ("p4_secondary_estimates", "landmark_inpatient", "Landmark and inpatient modeled outcomes"),
        ("p4_landmark_audit", "landmark_flow", "Landmark risk set, BP availability and G2 reversal update"),
        ("p4_exploratory_estimates", "aim3_exploratory", "MC-MED exploratory P subtype and SBP burden estimates"),
        ("p4_exploratory_burden_cutpoints", "aim3_burden_cutpoints", "SBP burden cutpoints and group counts"),
        ("p5_s33_consistency", "aim3_prior_epistaxis_history", "S3.3 consistency sensitivity for prior epistaxis history"),
        ("p5_treatment_features", "aim4_treatment_characteristics", "Characteristics of treated and untreated visits"),
        ("p5_treatment_selection_models", "aim4_treatment_selection", "Descriptive probability of antihypertensive treatment"),
        ("p5_exposure_mcmed_no_administration", "aim4_no_documented_administration", "MC-MED medication orders without documented administration"),
        ("p5_matched_response", "aim4_matched_response", "Matched BP response and B-class control sensitivity"),
        ("p5_matching_recovery", "aim4_matched_recovery", "Normal BP and recovery at the 120-minute target"),
        ("p5_matching_mixedlm", "aim4_matched_mixed_model", "Matched BP response mixed model cross-check"),
        ("p5_treatment_return_7d", "aim4_crude_return", "Descriptive treatment RRs, Katz approximate CI and Fisher exact test"),
        ("p6_bclass_overlap_estimate", "aim4_bclass_overlap", "Exploratory overlap-weighted B-class RR"),
        ("p6_bclass_overlap_balance", "aim4_bclass_overlap_balance", "Covariate balance before and after overlap weighting"),
        ("p6_bclass_overlap_weights", "aim4_bclass_overlap_diagnostics", "Overlap weights and effective sample sizes"),
        ("p6_bclass_overlap_missingness", "aim4_bclass_overlap_missingness", "Overlap-weighting covariate availability"),
    )
    for stem, output, title in supplements:
        note = ("Crude RR intervals are Katz approximate 95% CIs. Two-sided Fisher exact p values are reported separately."
                if stem == "p5_treatment_return_7d" else "")
        writer.copy_table(TABLES / f"{stem}.csv", f"supp_{output}.csv", f"Supplementary table. {title}", note)
    table_s11(writer)
    deviation_table(writer)


def run():
    FINAL.mkdir(parents=True, exist_ok=True)
    write_grouped_race()
    write_g6_comorbidities()
    writer = TableWriter()
    for table in (table1, table2, table3, table4, table5):
        table(writer)
    supplementary_tables(writer)
    writer.finish()
    return len(writer.index), len(writer.sources)


if __name__ == "__main__":
    try:
        tables, references = run()
    except Exception as error:
        print(f"P6 table assembly failed ({type(error).__name__}).", file=sys.stderr)
        raise SystemExit(1) from None
    print(f"P6 tables completed: {tables} aggregate tables and {references} source-cell references.")
