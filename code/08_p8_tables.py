"""Assemble the OHNS supplementary tables from approved aggregate files."""

import argparse
import csv
from math import isfinite
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "outputs" / "tables"
FINAL = TABLES / "final"
SUPPLEMENT = ROOT / "manuscript" / "supplement"
MANUSCRIPT_V2 = ROOT / "manuscript" / "draft_v2_OHNS_tagged.md"
NA = "—"
MCMED = "MC-MED"
MIMIC = "MIMIC-IV-ED"
PHENOTYPES = {
    "N": "normal on arrival", "T": "transiently elevated", "P": "persistently elevated",
    "N60": "normal on arrival at the landmark", "T60": "transiently elevated at the landmark",
    "P60": "persistently elevated at the landmark", "P-rebound": "P-rebound",
    "P-never": "persistently elevated without an intervening normal reading",
}
LABELS = {
    "C1": "BP cohort", "C2": "return cohort", "C3": "follow-up cohort", "Both centers": "Pooled",
    "known HTN": "documented prior hypertension", "prior_htn": "documented prior hypertension",
    "hemostasis_level2": "local hemostatic treatment", "A": "class A", "B": "class B",
    "B + alpha-1": "class B plus alpha-1", "exclude_b_controls": "Exclude class B-treated controls",
    "primary": "Primary", "Same-day HTN": "Same-day documentation of hypertension",
    "Prior HTN": "Documented prior hypertension", "Prior hypertension": "Documented prior hypertension",
    "Exploratory P subtype": "Exploratory persistently elevated subtypes",
}
ACRONYMS = {"bp", "sbp", "dbp", "ci", "ed", "rr", "pmh", "icd", "inr", "bic", "aic", "iqr", "htn", "txa", "npv", "ppv", "g1", "g2"}
FORBIDDEN_HEADERS = {"csn", "mrn", "visit_id", "patient_id", "stay_id", "subject_id", "hadm_id"}
GENERAL_FIELDS = (
    "Section", "Center", "Analysis or cohort", "Group or contrast", "Measure",
    "Sample or count", "Summary", "Context", "Status",
)


def aggregate_rows(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = tuple(reader.fieldnames or ())
        if not fields or any(column.lower() in FORBIDDEN_HEADERS for column in fields):
            raise ValueError("An expected aggregate input is missing a header or contains an identifier column.")
        rows = list(reader)
    if not rows or any(None in row or any(value is None for value in row.values()) for row in rows):
        raise ValueError("An expected aggregate input is empty or structurally inconsistent.")
    return fields, rows


def source(stem, final=True):
    return aggregate_rows((FINAL if final else TABLES) / f"{'supp_' if final else ''}{stem}.csv")[1]


def nval(value):
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if isfinite(result) else None


def text(value):
    return NA if value is None or str(value).strip().lower() in ("", "nan", "na", "none") else str(value)


def num(value, decimals=0):
    value = nval(value)
    return NA if value is None else f"{value:,.{decimals}f}".replace("-", "−")


def pct(value, fraction=False):
    value = nval(value)
    return NA if value is None else f"{value * (100 if fraction else 1):.1f}%"


def counts(numerator, denominator):
    return NA if nval(numerator) is None or nval(denominator) is None else f"{num(numerator)}/{num(denominator)}"


def count_pct(count, percent, fraction=False):
    number, percent = num(count), pct(percent, fraction)
    return number if percent == NA else f"{number} ({percent})"


def span(low, high, decimals=2, percent=False):
    a, b = nval(low), nval(high)
    if a is None or b is None:
        return NA
    join = " to " if min(a, b) < 0 else "-"
    suffix = "%" if percent else ""
    return f"{num(a, decimals)}{join}{num(b, decimals)}{suffix}"


def estimate_ci(estimate, low, high, decimals=2, fraction=False, percent=False):
    value = nval(estimate)
    if value is None:
        return NA
    factor = 100 if fraction else 1
    point = f"{num(value * factor, decimals)}{'%' if percent else ''}"
    a, b = nval(low), nval(high)
    interval = span(a * factor, b * factor, decimals, percent) if a is not None and b is not None else NA
    return point if interval == NA else f"{point} ({interval})"


def median_iqr(median, q1, q3, decimals=1):
    return estimate_ci(median, q1, q3, decimals)


def p_value(value):
    value = nval(value)
    if value is None:
        return NA
    if value < 0.001:
        return "P < .001"
    digits = 2 if value >= 0.01 else 3
    return f"P = {f'{value:.{digits}f}'.removeprefix('0')}"


def readable(field):
    terms = field.replace("_", " ").split()
    if not terms:
        return NA
    words = [word.upper() if word.lower() in ACRONYMS else word.lower() for word in terms]
    return " ".join((words[0][0].upper() + words[0][1:], *words[1:]))


def terminology(value):
    value = text(value)
    if value in PHENOTYPES:
        return PHENOTYPES[value]
    if value in LABELS:
        return LABELS[value]
    parts = value.split(" vs ")
    if len(parts) == 2 and all(part in PHENOTYPES for part in parts):
        return " vs ".join(PHENOTYPES[part] for part in parts)
    for label, name in (("C1", LABELS["C1"]), ("C2", LABELS["C2"]), ("C3", LABELS["C3"])):
        value = " ".join(name if word == label else word for word in value.split(" "))
    return (value.replace("Recorded prior hypertension", "Documented prior hypertension")
                 .replace("Known hypertension", "Documented prior hypertension").replace(" HTN", " hypertension"))


def write_csv(path, fields, rows):
    rows = list(rows)
    if not rows or any(len(row) != len(fields) for row in rows):
        raise ValueError(f"Unexpected shape for {path.name}.")
    path.parent.mkdir(parents=True, exist_ok=True)
    pending = path.with_suffix(".csv.tmp")
    try:
        with pending.open("w", encoding="utf-8-sig", newline="") as handle:
            output = csv.writer(handle, lineterminator="\n")
            output.writerow(fields)
            output.writerows(tuple(text(value) for value in row) for row in rows)
        pending.replace(path)
    finally:
        pending.unlink(missing_ok=True)


def general(section, center=None, analysis=None, group=None, measure=None, sample=None, summary=None, context=None, status=None):
    return (section, terminology(center), terminology(analysis), terminology(group), terminology(measure),
            sample, summary, context, status)


def write_general(directory, table, rows):
    write_csv(directory / f"TableS{table}.csv", GENERAL_FIELDS, rows)


def sample(visits, patients=None, events=None):
    parts = [f"{num(value)} {label}" for label, value in (("visits", visits), ("patients", patients), ("events", events))
             if nval(value) is not None]
    return ", ".join(parts) if parts else NA


def format_range_summary(value):
    # P6 combines verified estimates into a single display cell with an en dash.
    value = text(value)
    if "–" not in value:
        return value
    left, separator, interval = value.partition(" (")
    if not separator:
        return value.replace("–", "-")
    a, b = interval.rstrip(")").split("–", maxsplit=1)
    join = " to " if a.startswith(("-", "−")) or b.startswith(("-", "−")) else "-"
    return f"{left.replace('-', '−')} ({a.replace('-', '−')}{join}{b.replace('-', '−')})"


def table1(directory):
    rows = []
    for row in source("flow"):
        rows.append(("Index flow", row["site"], num(row["step"]), terminology(row["stage"]), num(row["visits"]),
                     num(row["patients"]), num(row["visits_removed"]), num(row["patients_removed"])))
    for row in source("flow_subcohorts"):
        reason = f'{terminology(row["transition"])} / {terminology(row["reason"])}'
        rows.append(("Subcohort exclusions", row["site"], NA, reason, NA, NA, num(row["excluded_visits"]), num(row["affected_patients"])))
    for row in source("flow_c3"):
        for name, visit, person, removed in (
            ("BP cohort", "c1_visits", "c1_patients", None),
            ("Potential follow-up cohort before incomplete visit exclusion", "c3_visits_before_incomplete", "c3_patients_before_incomplete", None),
            ("Incomplete visits excluded", None, None, "incomplete_c3_visits"),
            ("Eligible visits before single index selection", "c3_eligible_visits", "c3_eligible_patients", None),
            ("First eligible follow-up cohort index", None, "first_c3_index_patients", None),
        ):
            rows.append(("Follow-up cohort flow", row["site"], readable(row["index_variant"]), name, num(row.get(visit)),
                         num(row.get(person)), num(row.get(removed)), NA))
    write_csv(directory / "TableS1.csv", ("Section", "Center", "Step or scenario", "Stage or reason", "Visits", "Patients", "Excluded visits", "Affected patients"), rows)


def table2(directory):
    rows = []
    for row in source("hypertension_sources"):
        rows.append((row["site"], terminology(row["cohort"]), terminology(row["source"]), num(row["visits"]),
            num(row["prior_htn_visits"]), count_pct(row["source_visits"], row["source_pct"]),
            span(row["source_wilson95_lo_pct"], row["source_wilson95_hi_pct"], 1, True),
            count_pct(row["source_only_visits"], row["source_only_pct_of_all"]), pct(row["source_only_pct_of_htn"])))
    write_csv(directory / "TableS2.csv", ("Center", "Cohort", "Documented hypertension source", "All visits",
        "Visits with documented hypertension", "Source visits, n (%)", "Wilson 95% CI", "Source only, n (%) of all visits",
        "Source only, % of visits with documented hypertension"), rows)


def table3(directory):
    rows = []
    for row in source("sparse_sampling"):
        percentage = row["metric"].startswith("Sensitivity")
        rows.append(general("Agreement and detection", measure=row["metric"], sample=f'{num(row["replications"])} repetitions',
            summary=estimate_ci(row["mean"], row["p025"], row["p975"], 1 if percentage else 3, percentage, percentage),
            context="Monte Carlo mean and central 95% interval"))
    for row in source("sparse_sampling_confusion"):
        rows.append(general("Classification matrix", group=f'{terminology(row["full_data"])} to {terminology(row["sparse"])}',
            measure="Average classification count and row percentage", sample=num(row["pooled_count"]),
            summary=f'{num(row["mean_per_replication"], 2)} ({pct(row["row_pct"])})', context="Sample column is pooled across repetitions"))
    write_general(directory, 3, rows)


def table4(directory):
    rows = []
    for row in source("lmm_satterthwaite"):
        results = [estimate_ci(row[f"{prefix}_estimate"], row[f"{prefix}_ci_low"], row[f"{prefix}_ci_high"], 3)
                   for prefix in ("original", "satterthwaite")]
        rows.append(("Coefficients", terminology(row["term"]), *results,
                     num(row["original_df"], 1), num(row["satterthwaite_df"], 1),
                     p_value(row["original_p_value"]), p_value(row["satterthwaite_p_value"]), NA))
    diagnostics = source("lmm_diagnostics")
    if len(diagnostics) != 2:
        raise ValueError("The audited linear mixed-model comparison requires the two fit diagnostics.")
    for field in diagnostics[0]:
        if field == "model":
            continue
        if field == "joint_interaction_p_value":
            values = [p_value(row[field]) for row in diagnostics]
        else:
            values = [num(row[field], 0 if field.endswith("_n") or field == "warning_count" else 3)
                      if nval(row[field]) is not None else text(row[field]) for row in diagnostics]
        rows.append(("Diagnostics", readable(field), *values, NA, NA, NA, NA, NA))
    write_csv(directory / "TableS4.csv", ("Section", "Term or diagnostic", "Nested model estimate (95% CI)",
        "Satterthwaite model estimate (95% CI)", "Nested model df", "Satterthwaite model df", "Nested model P value",
        "Satterthwaite model P value", "Notes"), rows)


def table5(directory):
    rows = []
    for row in source("lcmm_models"):
        state = "Selected" if row["selected"].lower() == "true" else "Eligible" if row["eligible"].lower() == "true" else "Not eligible"
        rows.append(("Model selection", num(row["k"]), NA, num(row["visit_n"]), num(row["patient_n"]), num(row["reading_n"]),
            num(row["bic"], 2), NA, pct(row["min_class_proportion"], True), num(row["min_mean_assigned_posterior"], 3),
            num(row["iterations"]), row["converged"], state))
    for row in source("lcmm_classes"):
        state = "Selected model" if row["selected"].lower() == "true" else "Candidate model"
        rows.append(("Assigned classes", num(row["k"]), num(row["class"]), num(row["visits"]), NA, NA, NA,
            count_pct(row["visits"], row["proportion"], True), NA, num(row["mean_assigned_posterior"], 3), NA, NA, state))
    write_csv(directory / "TableS5.csv", ("Section", "Classes fitted", "Assigned class", "Visits", "Patients", "Readings",
        "BIC", "Assigned visits, n (%)", "Minimum assigned class, %", "Mean assigned posterior or model minimum",
        "Iterations", "Converged", "Selection status"), rows)


def table6(directory):
    rows = []
    for row in source("aim2_followup_availability"):
        for label, numerator, denominator, ratio, confidence in (
            ("Primary follow-up", "with_followup", "patients", "followup_pct", ("followup_wilson95_low_pct", "followup_wilson95_high_pct")),
            ("Strict follow-up", "with_strict_followup", "patients", "strict_followup_pct", None),
            ("Follow-up among documented medication lists", "medication_list_with_followup", "medication_list_patients", "medication_list_followup_pct", None),
        ):
            rows.append(("Availability", row["site"], readable(row["index_variant"]), label, NA, num(row[denominator]),
                count_pct(row[numerator], row[ratio]), NA, span(row[confidence[0]], row[confidence[1]], 1, True) if confidence else NA))
    for row in source("aim2_followup_baseline"):
        groups = []
        for prefix in ("with_followup", "without_followup"):
            summary = (median_iqr(row[f"{prefix}_median"], row[f"{prefix}_q1"], row[f"{prefix}_q3"])
                       if nval(row[f"{prefix}_median"]) is not None else count_pct(row[f"{prefix}_count"], row[f"{prefix}_pct"]))
            groups.append(f'{summary} [n = {num(row[f"{prefix}_nonmissing"])}]' if summary != NA else NA)
        rows.append(("Characteristics by follow-up status", row["site"], NA, readable(row["characteristic"]),
                     terminology(row["level"]), NA, *groups, num(row["absolute_standardized_difference"], 3)))
    write_csv(directory / "TableS6.csv", ("Section", "Center", "Index scenario", "Characteristic or criterion", "Level",
        "Eligible for availability estimate", "With follow-up", "Without follow-up", "Wilson 95% CI or absolute standardized difference"), rows)


def table7(directory):
    rows = []
    primary = [("Primary", row) for row in source("aim2_primary_detail")]
    key_fields = ("site", "analysis", "estimand", "group_or_contrast")
    sensitivity = {tuple(row[field] for field in key_fields): row for row in source("p3_sensitivity", False)}
    selected = []
    for row in source("s11_sensitivity_all"):
        if row["Aim"] == "Aim 2":
            key = (row["Center"], row["Analysis"], row["Measure"], row["Group or contrast"])
            selected.append((row["SAP item"], sensitivity[key]))
    if len(selected) != len(sensitivity):
        raise ValueError("The primary and sensitivity inputs for subsequent hypertension do not reconcile.")
    for label, row in (*primary, *selected):
        risk = row["estimand"] != "Adjusted risk ratio"
        rows.append((label, row["site"], terminology(row["analysis"]), terminology(row["estimand"]), terminology(row["group_or_contrast"]),
            num(row["eligible_patients"]), num(row["followup_patients"]), counts(row["numerator"], row["denominator"]),
            counts(row["reference_numerator"], row["reference_denominator"]),
            estimate_ci(row["estimate"], row["ci95_low"], row["ci95_high"], 1 if risk else 2, risk, risk), row["status"]))
    write_csv(directory / "TableS7.csv", ("Protocol item", "Center", "Analysis", "Measure", "Group or contrast", "Eligible patients",
        "Followed patients", "Events or positive results / denominator", "Reference events or positive results / denominator",
        "Estimate (95% CI)", "Status"), rows)


def table8(directory):
    candidates = [row for row in source("p4_sensitivity_estimates", False)
                  if row["site"] == "Both centers" and row["contrast"] == "P vs T"
                  and row["outcome"] == "7-day epistaxis return" and nval(row["rr"]) is not None]
    low = min(candidates, key=lambda row: nval(row["rr"]))
    high = max(candidates, key=lambda row: nval(row["rr"]))
    if (low["scenario"], high["scenario"]) != ("S11-4", "S11-1") or (num(low["rr"], 2), num(high["rr"], 2)) != ("0.86", "1.22"):
        raise ValueError("The Table S8 sensitivity range disagrees with the manuscript-audited end points.")
    rows = []
    marked = []
    for row in source("s11_sensitivity_all"):
        if row["Aim"] != "Aim 3":
            continue
        note = NA
        if row["Center"] == "Both centers" and row["Group or contrast"] == "P vs T" and row["Measure"] == "7-day epistaxis return":
            if row["SAP item"] in (low["scenario"], high["scenario"]):
                chosen = low if row["SAP item"] == low["scenario"] else high
                direction = "Lower" if chosen is low else "Upper"
                note = f"{direction} end point cited in main text"
                if row["Estimate or summary"].split(" ")[0] != num(chosen["rr"], 2):
                    raise ValueError("A marked Table S8 estimate disagrees with its unrounded source.")
                marked.append(direction)
        rows.append((row["SAP item"], terminology(row["Center"]), terminology(row["Analysis"]), row["Measure"],
            terminology(row["Group or contrast"]), row["Sample"], row["Events or exposed count"],
            format_range_summary(row["Estimate or summary"]), row["Status"], note))
    if sorted(marked) != ["Lower", "Upper"]:
        raise ValueError("Both range end points must appear exactly once in Table S8.")
    write_csv(directory / "TableS8.csv", ("Protocol item", "Center", "Analysis", "Outcome or statistic", "Contrast", "Sample",
        "Events or count", "Estimate (95% CI) or statistic", "Status", "Main-text range"), rows)


def table9(directory):
    rows = [(row["site"], terminology(row["contrast"]), sample(row["visits"], row["patients"], row["events"]),
             counts(row["p_events"], row["p_visits"]), counts(row["t_events"], row["t_visits"]),
             num(row["known_prior_epistaxis_visits"]), num(row["missing_prior_epistaxis_visits"]),
             estimate_ci(row["p4_main_estimate"], row["p4_main_ci95_low"], row["p4_main_ci95_high"]),
             estimate_ci(row["estimate"], row["ci95_low"], row["ci95_high"]), row["status"])
            for row in source("aim3_prior_epistaxis_history")]
    write_csv(directory / "TableS9.csv", ("Center", "Contrast", "Analytic sample", "Persistently elevated events / visits",
        "Transiently elevated events / visits", "Prior epistaxis visits recorded", "Missing history",
        "Locked estimate (95% CI)", "Estimate with prior history (95% CI)", "Status"), rows)


def table10(directory):
    rows = []
    for row in source("landmark_flow"):
        for label, visits, events in (
            ("Clinical risk set", "landmark_risk_set_visits", "landmark_events"),
            ("Early phenotype available", "phenotype_observed_visits", "phenotype_observed_events"),
            ("Early phenotype unavailable", "phenotype_missing_visits", "phenotype_missing_events"),
            ("Arrival SBP available", "arrival_sbp_observed_visits", "arrival_sbp_observed_events"),
        ):
            rows.append(general("Landmark flow", MCMED, row["analysis"], measure=label,
                sample=counts(row[events], row[visits]), context="Escalations / visits"))
        for field in ("g2_named_medication_order_records", "g2_medication_orders_in_index_ed",
                      "g2_additional_early_event_exclusions", "g2_additional_late_events_in_revised_risk_set"):
            rows.append(general("Andexanet rule audit", MCMED, row["analysis"], measure=readable(field), sample=num(row[field])))
    for row in source("landmark_inpatient"):
        is_inpatient = row["site"] == MIMIC
        group = (terminology(row["contrast"]) if not row["reference"] else
                 f'{terminology(row["contrast"])} vs {terminology(row["reference"])}')
        n = sample(row["modeled_visits"], row["modeled_patients"], row["modeled_events"])
        contrast_counts = counts(row["fitted_contrast_events"], row["fitted_contrast_visits"])
        reference_counts = counts(row["fitted_reference_events"], row["fitted_reference_visits"])
        context = f'Clinical risk set {counts(row["raw_risk_set_events"], row["raw_risk_set_visits"])} events/visits'
        if contrast_counts != NA:
            context += f". Modeled events/visits {contrast_counts} vs {reference_counts}"
        adjust = {"landmark": "Adjusted for age, sex, antithrombotic use and documented hypertension",
                  "age": "Age-adjusted", "none": "Crude"}[row["adjustment"]]
        rows.append(general("Inpatient procedures" if is_inpatient else "Escalation-of-care models", row["site"],
            row["analysis"], group, f"{adjust} risk ratio (95% CI)", n,
            estimate_ci(row["estimate"], row["ci95_low"], row["ci95_high"]), context, row["status"]))
    write_general(directory, 10, rows)


def table11(directory):
    rows = []
    for row in source("aim3_exploratory"):
        group = f'{terminology(row["contrast"])} vs {terminology(row["reference"])}'
        groups = (f'{counts(row["fitted_contrast_events"], row["fitted_contrast_visits"])} vs '
                  f'{counts(row["fitted_reference_events"], row["fitted_reference_visits"])} events/visits')
        rows.append(general("Exploratory return models", MCMED, row["analysis"], group, "Crude risk ratio (95% CI)",
            sample(row["modeled_visits"], row["modeled_patients"], row["modeled_events"]),
            estimate_ci(row["estimate"], row["ci95_low"], row["ci95_high"]), groups, row["status"]))
    for row in source("aim3_burden_cutpoints"):
        lower, upper = num(row["strict_lower_mmhg"], 5), num(row["inclusive_upper_mmhg"], 5)
        rows.append(general("SBP burden grouping", MCMED, "Return cohort", row["tertile"], "Group boundary, mmHg",
            counts(row["events"], row["visits"]), f"Lower exclusive {lower}. Upper inclusive {upper}", "Events / visits"))
    write_general(directory, 11, rows)


def table12(directory):
    rows = []
    for row in source("p8_posthoc_severity_descriptive", False):
        n = sample(row["visits"], row["patients"])
        group = terminology(row["phenotype"])
        rows.append(general("Post hoc description", row["site"], "Return cohort", group, "Local hemostatic treatment, n (%)", n,
            count_pct(row["local_hemostatic_visits"], row["local_hemostatic_percent"]),
            f'Wilson 95% CI {span(row["local_wilson_ci95_low"], row["local_wilson_ci95_high"], 1, True)}'))
        rows.append(general("Post hoc description", row["site"], "Return cohort", group, "First ED hemoglobin, g/dL, median (IQR)", n,
            median_iqr(row["hemoglobin_median_g_dl"], row["hemoglobin_q1_g_dl"], row["hemoglobin_q3_g_dl"]),
            f'Available {count_pct(row["hemoglobin_observed_visits"], row["hemoglobin_available_percent"])}'))
    for row in source("p8_posthoc_severity", False):
        n = sample(row["modeled_visits"], row["modeled_patients"], row["modeled_events"])
        if n == NA:
            n = sample(row["eligible_visits"])
        groups = ". ".join(f'{PHENOTYPES[k.upper()]} {counts(row[f"{k}_events"], row[f"{k}_visits"])}'
                           for k in ("n", "t", "p"))
        epv = nval(row["events_per_parameter"])
        if epv is not None:
            groups += f'. Events per parameter {num(epv, 1)}'
        if text(row["notes"]) != NA:
            groups += f'. {row["notes"]}'
        if text(row["stratum"]) != NA:
            groups = f'{row["stratum"]}. {groups}'
        state = f'{row["fit_status"]} / {row["stability"]}' if text(row["stability"]) != NA else row["fit_status"]
        rows.append(general("Post hoc outcome models", row["site"], row["analysis"], row["contrast"], row["model"], n,
            estimate_ci(row["rr"], row["ci95_low"], row["ci95_high"]), groups, state))
    for row in source("p8_hemoglobin_coverage", False):
        rows.append(general("Hemoglobin availability", row["site"], "Return cohort", measure="First ED hemoglobin available, n (%)",
            sample=sample(row["return_visits"]), summary=count_pct(row["hemoglobin_available_visits"], row["hemoglobin_available_percent"])))
    write_general(directory, 12, rows)


def table13(directory):
    rows = []
    for row in source("aim4_treatment_characteristics"):
        quant = nval(row["median"]) is not None
        summary = median_iqr(row["median"], row["q1"], row["q3"]) if quant else count_pct(row["count"], row["percent"])
        interval = NA if quant else span(row["ci95_percent_low"], row["ci95_percent_high"], 1, True)
        rows.append(general("Treatment characteristics", row["site"], row["treatment_definition"], row["group"],
            terminology(readable(row["characteristic"])), f'{num(row["observed"])} observed of {num(row["visits"])} visits', summary, interval))
    for row in source("aim4_treatment_selection"):
        details = f'{num(row["treated_visits"])} treated. {num(row["untreated_visits"])} untreated. {p_value(row["p_value"])}'
        rows.append(general("Descriptive treatment models", "Pooled", row["treatment_definition"], row["method"],
            terminology(readable(row["predictor"])), sample(row["visits"], row["patients"], row["events"]),
            estimate_ci(row["estimate"], row["ci95_low"], row["ci95_high"]), details, row["status"]))
    write_general(directory, 13, rows)


def table14(directory):
    rows = []
    for row in source("p5_matching_flow", False):
        for field, value in row.items():
            if field not in ("scenario", "site"):
                rows.append(general("Matching flow", row["site"], row["scenario"], measure=terminology(readable(field)),
                    sample=num(value, 1 if field.startswith("median_") else 0)))
    for row in source("aim4_matched_response"):
        group = f'{num(row["minutes"])} minutes'
        detail = (f'Mean change {num(row["treated_mean_sbp_change"], 1)} vs {num(row["control_mean_sbp_change"], 1)} mmHg. '
                  f'{num(row["control_slots"])} control slots and {num(row["unique_control_visits"])} distinct control visits')
        rows.append(general("Matched-set difference", row["site"], row["scenario"], group, "Treated minus control SBP change, mmHg",
            f'{num(row["matched_sets"])} sets', estimate_ci(row["treated_minus_control_change"], row["ci95_low"], row["ci95_high"], 1),
            detail, row["status"]))
    for row in source("p5_matching_group_means", False):
        rows.append(general("Matched group means", row["site"], row["scenario"], row["group"], f'Mean SBP at {num(row["minutes"])} minutes, mmHg',
            f'{num(row["matched_sets"])} sets, {num(row["readings"])} readings',
            estimate_ci(row["mean_sbp"], row["ci95_low"], row["ci95_high"], 1)))
    for row in source("aim4_matched_recovery"):
        rows.append(general("Documented normalization", row["site"], row["scenario"], row["group"], row["population"],
            counts(row["normal_observations"], row["observations"]),
            estimate_ci(row["set_weighted_percent_normal"], row["ci95_low_pct"], row["ci95_high_pct"], 1, percent=True),
            f'{num(row["matched_sets"])} matched sets. Sample is normal readings / observed readings'))
    for row in source("aim4_matched_mixed_model"):
        rows.append(general("Mixed-model check", row["site"], row["scenario"], f'{num(row["minutes"])} minutes',
            "Treated minus control SBP change, mmHg", f'{num(row["matched_sets"])} sets, {num(row["measurements"])} measurements',
            estimate_ci(row["treated_minus_control_change"], row["ci95_low"], row["ci95_high"], 1),
            f'Matched-set variance {num(row["estimated_matched_set_variance"], 3)}. Residual variance {num(row["residual_variance"], 3)}', row["status"]))
    write_general(directory, 14, rows)


def table15(directory):
    rows = []
    for row in source("aim4_crude_return"):
        group = (f'Treated {counts(row["treated_events"], row["treated_visits"])} ({pct(row["treated_risk_percent"])}) vs '
                 f'untreated {counts(row["untreated_events"], row["untreated_visits"])} ({pct(row["untreated_risk_percent"])})')
        rows.append(general("Crude return estimates", row["site"], row["treatment_definition"], group,
            "Crude risk ratio (Katz approximate 95% CI)", sample(row["visits"], row["patients"]),
            estimate_ci(row["crude_rr"], row["katz_ci95_low"], row["katz_ci95_high"]),
            f'Two-sided Fisher exact {p_value(row["fisher_two_sided_p"])}', row["status"]))
    for row in source("aim4_bclass_overlap"):
        summary = estimate_ci(row["estimate"], row["ci95_low"], row["ci95_high"])
        details = f'Weighted risk {pct(row["weighted_treated_risk_percent"])} vs {pct(row["weighted_untreated_risk_percent"])}'
        rows.append(general("Overlap weighting", "Pooled", "class B", "Treated vs untreated", "Exploratory weighted risk ratio (95% CI)",
                            sample(row["visits"], row["patients"]), summary, details, row["status"]))
    for row in source("aim4_bclass_overlap_balance"):
        labels = ("Unweighted standardized difference", "Weighted standardized difference", "Unweighted pooled SD")
        for label, field in zip(labels, ("unweighted_smd", "weighted_smd", "unweighted_pooled_sd")):
            value = nval(row[field])
            display = f"{value:.3g}".replace("-", "−") if value is not None else NA
            rows.append(general("Covariate balance", "Pooled", terminology(readable(row["covariate"])), row["level"], label, summary=display))
        for prefix in ("unweighted", "weighted"):
            interval = f'{num(row[prefix+"_treated_mean"], 3)} vs {num(row[prefix+"_untreated_mean"], 3)}'
            rows.append(general("Covariate balance", "Pooled", terminology(readable(row["covariate"])), row["level"],
                                f"{prefix.capitalize()} covariate mean, treated vs untreated", summary=interval))
    for row in source("aim4_bclass_overlap_diagnostics"):
        for field, value in row.items():
            if field not in ("site", "treatment_arm", "visits", "patients", "events"):
                rows.append(general("Weight diagnostics", row["site"], "class B", row["treatment_arm"], readable(field),
                                    sample(row["visits"], row["patients"], row["events"]), num(value, 3)))
    for row in source("aim4_bclass_overlap_missingness"):
        rows.append(general("Covariate availability", "Pooled", "class B", measure=terminology(readable(row["covariate"])),
            sample=f'{num(row["eligible_visits"])} eligible, {num(row["modeled_visits"])} modeled',
            summary=count_pct(row["missing_visits"], row["missing_pct"]), context="Missing visits, n (%)"))
    write_general(directory, 15, rows)


def table16(directory):
    rows = []
    for row in source("aim4_no_documented_administration"):
        if row["exposure_definition"] != "Class A":
            continue
        rows.append((row["site"], terminology(row["population"]), num(row["eligible_visits"]),
            num(row["qualifying_medication_orders_placed_in_ed"]), num(row["orders_without_documented_administration"]),
            num(row["visits_with_orders_without_documented_administration"]), num(row["patients_with_orders_without_documented_administration"]),
            num(row["visits_with_no_admin_order_and_no_in_class_ed_dose"]), num(row["visits_with_no_admin_order_and_in_class_ed_dose"]),
            num(row["orders_with_unparseable_administration_timestamp"]), num(row["orders_with_administration_outside_ed"])))
    fields = ("Center", "Cohort", "Eligible visits", "Class A orders issued in the ED", "Orders without documented administration",
              "Visits with such orders", "Patients with such orders", "Visits without a class A ED dose",
              "Visits with a class A ED dose", "Orders with invalid administration time", "Orders with administration outside the ED")
    write_csv(directory / "TableS16.csv", fields, rows)


def english_deviations(path):
    section = None
    rows = []
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        if line.startswith("## "):
            section = line.removeprefix("## ")
        if not line.startswith("|"):
            continue
        cells = [part.strip() for part in line.strip().strip("|").split("|")]
        if cells[0] == "Date" or set(cells[0]) == {"-"}:
            continue
        if not section or len(cells) not in (3, 4):
            raise ValueError("An English protocol log entry has an unsupported shape.")
        rows.append((section, *cells, *([NA] if len(cells) == 3 else [])))
    if not any("post hoc" in row[0].lower() for row in rows):
        raise ValueError("The English protocol log lacks the P8 post hoc entry.")
    return rows


def table17(directory, translation):
    write_csv(directory / "TableS17.csv", ("Section", "Date", "Basis", "Decision and rationale", "Additional rationale"),
              english_deviations(translation))


DICTIONARIES = (
    ("supp_dictionary_mcmed_orders_mapping.csv", "MC_MED_ED_orders.csv", True),
    ("supp_dictionary_mimic_pyxis_mapping.csv", "MIMIC_IV_ED_Pyxis.csv", True),
    ("supp_dictionary_mcmed_meds_mapping.csv", "MC_MED_home_medications.csv", True),
    ("supp_dictionary_mimic_medrecon_mapping.csv", "MIMIC_IV_ED_medication_reconciliation.csv", True),
    ("supp_dictionary_mcmed_labs_mapping.csv", "MC_MED_lab_components.csv", True),
    ("supp_dictionary_mcmed_disposition_mapping.csv", "MC_MED_ED_dispositions.csv", True),
    ("supp_dictionary_mcmed_reversal_g1_mapping.csv", "MC_MED_reversal_review_G1.csv", True),
    ("supp_dictionary_mimic_inpatient_procedures.csv", "MIMIC_IV_inpatient_nasal_hemostatic_codes.csv", True),
    ("supp_baseline_race_mapping.csv", "Race_ethnicity_mapping.csv", True),
    ("p8_pyxis_local_hemostatic_audit.csv", "MIMIC_IV_ED_local_hemostatic_Pyxis_audit.csv", False),
    ("p8_hemoglobin_audit.csv", "MIMIC_IV_ED_hemoglobin_item_audit.csv", False),
)
DICTIONARY_HEADERS = {
    "aht_acute": "Class A agent", "aht_any": "Class B agent", "aht_alpha1_sensitivity": "Alpha-1 blocker sensitivity definition",
    "aht_ed_category": "ED antihypertensive category", "aht_ed_exclusion_reason": "Reason for ED antihypertensive exclusion",
    "aht_home": "Home antihypertensive medication", "aht_home_sensitivity_only": "Home antihypertensive sensitivity definition only",
    "aht_home_exclusion_reason": "Reason for home antihypertensive exclusion", "acap_exclusion_reason": "Reason for antithrombotic exclusion",
    "hemostasis_level": "Hemostatic treatment tier", "hemostasis_exclusion_reason": "Reason for hemostatic classification exclusion",
    "hemostasis_category": "Hemostatic treatment category", "g1_hemostasis_level": "Hemostatic treatment tier at G1",
    "g1_reversal_exclusion_reason": "Reason for reversal exclusion at G1", "txa_route": "Tranexamic acid route",
    "ed_dispo": "Recorded ED disposition", "etcdescription": "Recorded medication class description",
    "component_name": "Recorded laboratory component label", "med_class": "Recorded medication class", "med_subclass": "Recorded medication subclass",
    "procedure_name": "Recorded procedure or medication label", "level": "Combined race or ethnicity category",
    "start_to_discharge_records": "Pyxis records between arrival and ED departure",
    "start_to_discharge_visits": "Visits with Pyxis record between arrival and ED departure",
}
ESCALATION_LABELS = {"2": "Local hemostatic treatment", "3": "Escalation of care"}


def dictionary_value(field, value):
    if field in ("hemostasis_level", "g1_hemostasis_level"):
        return ESCALATION_LABELS.get(value, value)
    if field == "aht_ed_category":
        return {"acute": "class A", "any_only": "class B only", "alpha1_sensitivity_only": "alpha-1 sensitivity only"}.get(value, value)
    return value


def data1(directory):
    for filename, output, final in DICTIONARIES:
        fields, rows = aggregate_rows((FINAL if final else TABLES) / filename)
        headers = tuple(DICTIONARY_HEADERS.get(field.lower(), readable(field)) for field in fields)
        data = (tuple(dictionary_value(field, row[field]) for field in fields) for row in rows)
        write_csv(directory / output, headers, data)


MAP_RULES = (
    ("Deviations are listed in the Supplement", "Supplementary Methods 9 and Table S17",
     "supplement_v1.md and tables/TableS17.csv"),
    ("removed duplicate readings across sources (Supplement)", "Supplementary Methods 3", "supplement_v1.md"),
    ("prespecified sensitivity analyses (Supplement)", "Supplementary Methods 7 and Table S8",
     "supplement_v1.md and tables/TableS8.csv"),
    ("Supplementary Figure S2", "Figure S2", "outputs/figures/submission/FigureS2.tif and FigureS2.pdf"),
    ("underwent a nasal hemostatic procedure, without a clear association with phenotype (Supplement)",
     "Table S10", "tables/TableS10.csv"),
)


def mentions(line):
    lower = line.lower()
    cursor = 0
    while (start := lower.find("supplement", cursor)) >= 0:
        stop = start + (len("supplementary") if lower.startswith("supplementary", start) else len("supplement"))
        cursor = stop
        if (start == 0 or not lower[start - 1].isalnum()) and (stop == len(lower) or not lower[stop].isalnum()):
            yield start, stop


def supplement_map(path, manuscript):
    rows = []
    for line_no, line in enumerate(manuscript.read_text(encoding="utf-8-sig").splitlines(), start=1):
        for start, stop in mentions(line):
            choices = []
            for anchor, target, artifact in MAP_RULES:
                first = line.find(anchor)
                if first >= 0 and first <= start and stop <= first + len(anchor):
                    choices.append((anchor, target, artifact))
            if len(choices) != 1:
                raise ValueError("A supplementary mention in the manuscript could not be uniquely resolved.")
            excerpt, target, artifact = choices[0]
            rows.append((line_no, excerpt, target, artifact))
    if len(rows) != len(MAP_RULES):
        raise ValueError("The manuscript supplement map does not reconcile with the reviewed mentions.")
    write_csv(path, ("Manuscript line", "Original excerpt", "Supplementary item", "Linked artifact"), rows)


def translation_path():
    places = (ROOT / "derived" / "p8_deviations_en.md", ROOT / "docs" / "deviations_en.md",
              ROOT / "release" / "docs" / "deviations_en.md")
    return next((path for path in places if path.is_file()), places[1])


def run(directory=SUPPLEMENT, manuscript=MANUSCRIPT_V2, translation=None):
    translation = translation or translation_path()
    english_deviations(translation)
    table_directory = directory / "tables"
    for builder in (table1, table2, table3, table4, table5, table6, table7, table8, table9, table10,
                    table11, table12, table13, table14, table15, table16):
        builder(table_directory)
    table17(table_directory, translation)
    data1(directory / "data1")
    map_written = manuscript.is_file()
    if map_written:
        supplement_map(directory / "supplement_map.csv", manuscript)
    return map_written


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=SUPPLEMENT)
    parser.add_argument("--manuscript-v2", type=Path, default=MANUSCRIPT_V2)
    parser.add_argument("--deviations", type=Path)
    arguments = parser.parse_args()
    try:
        include_map = run(arguments.output_dir, arguments.manuscript_v2, arguments.deviations)
    except Exception as error:
        print(f"P8 supplementary tables could not be assembled ({type(error).__name__}).", file=sys.stderr)
        raise SystemExit(1) from None
    print("P8 supplementary tables and Supplementary Data 1 completed.")
    if include_map:
        print("The read-only manuscript's supplementary mentions were mapped.")
    else:
        print("Manuscript map requires the tagged manuscript input, which is absent in the code-only release.")
