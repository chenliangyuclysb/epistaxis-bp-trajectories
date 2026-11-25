"""Fit the frozen Aim 2 analyses and write aggregate tables only."""

import csv
from dataclasses import dataclass
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.tools.sm_exceptions import PerfectSeparationError, PerfectSeparationWarning

import config as cfg


TABLES = Path("outputs/tables")
SITES = ("MC-MED", "MIMIC-IV-ED")
PHENOTYPES = ("N", "T", "P")
INPUT_COLUMNS = (
    "site", "index_variant", "visit_id", "patient_id", "bp_phenotype", "age_years", "sex", "sbp_first",
    "triage_acuity_group", "has_medication_list", "has_followup", "htn_evidence",
    "has_followup_strict", "htn_evidence_strict",
)
OUTCOME_PREDICTORS = ("age_years", "sex", "sbp_first", "bp_phenotype")
SELECTION_PREDICTORS = (*OUTCOME_PREDICTORS, "triage_acuity_group", "site")
ADJUSTED_RR = "Adjusted risk ratio"
CONTRAST = "P vs T"
ESTIMATE_HEADER = (
    "site", "analysis", "estimand", "group_or_contrast", "eligible_patients", "followup_patients",
    "numerator", "denominator", "reference_numerator", "reference_denominator",
    "estimate", "ci95_low", "ci95_high", "status",
)


def write_rows(name, header, rows):
    path = TABLES / name
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".csv.tmp")
    try:
        with temporary.open("w", encoding="utf-8-sig", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(header)
            writer.writerows(rows)
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def load_input(con):
    frame = con.execute(f"SELECT {', '.join(INPUT_COLUMNS)} FROM p3_patient_followup "
                        "ORDER BY site, index_variant, patient_id").df()
    variant_sites = {"primary": SITES, "same_day_htn": (SITES[0],), "no_prior_ed": (SITES[1],)}
    permitted = {(site, variant) for variant, sites in variant_sites.items() for site in sites}
    pairs = set(zip(frame.site, frame.index_variant))
    if not pairs <= permitted or frame[["patient_id", "visit_id"]].isna().any().any():
        raise ValueError("Unexpected P3 index keys or site/index definition.")
    if frame.duplicated(["site", "index_variant", "patient_id"]).any():
        raise ValueError("Aim 2 requires one first eligible index per patient within each definition.")
    if not frame.bp_phenotype.isin(PHENOTYPES).all():
        raise ValueError("All P3 patients require a valid C1 blood-pressure phenotype.")
    for suffix in ("", "_strict"):
        observed = frame[f"has_followup{suffix}"]
        evidence = frame[f"htn_evidence{suffix}"]
        if observed.isna().any() or not observed.astype(bool).eq(evidence.notna()).all():
            raise ValueError("HTN evidence must be observed exactly for patients with eligible follow-up.")
    if (frame.has_followup_strict & ~frame.has_followup).any():
        raise ValueError("Strict follow-up must be a subset of primary follow-up.")
    for field in ("age_years", "sbp_first"):
        values = pd.to_numeric(frame[field], errors="coerce")
        frame[field] = values.where(np.isfinite(values))
    for field in ("sex", "triage_acuity_group"):
        values = frame[field].astype("string").str.strip()
        frame[field] = values.mask(values.eq(""))
    return frame


def wilson(numerator, denominator):
    if not denominator:
        return None, None, None
    proportion = numerator / denominator
    z = cfg.WILSON_95_Z
    scale = 1 + z**2 / denominator
    center = (proportion + z**2 / (2 * denominator)) / scale
    radius = z * np.sqrt(proportion * (1 - proportion) / denominator + z**2 / (4 * denominator**2)) / scale
    return proportion, max(0, center - radius), min(1, center + radius)


def analyze_missing(frame, fields, label, analysis, model, diagnostics):
    observed = frame[list(fields)].dropna()
    n = len(frame)
    missing_over_threshold = False
    for field in fields:
        missing = int(frame[field].isna().sum())
        percent = 100 * missing / n if n else None
        requires_imputation = percent is not None and percent > cfg.P3_COVARIATE_MI_TRIGGER_PERCENT
        missing_over_threshold |= requires_imputation
        diagnostics.append((label, analysis, model, n, len(observed), field, missing, percent,
                            "MI required by SAP; estimate withheld" if requires_imputation else "Complete-case eligibility"))
    return observed.index, missing_over_threshold


def design_matrix(frame, categorical):
    columns = {"intercept": np.ones(len(frame))}
    for name in ("age_years", "sbp_first"):
        values = frame[name].to_numpy(dtype=float)
        scale = values.std()
        # Scaling only improves numerical conditioning, without changing a fitted risk ratio.
        columns[name] = (values - values.mean()) / scale if scale > 0 else values * 0.0
    for name in categorical:
        levels = sorted(frame[name].unique())
        if name == "bp_phenotype" and "T" in levels:
            levels.remove("T")
            levels.insert(0, "T")
        for level in levels[1:]:
            columns[f"{name}={level}"] = frame[name].eq(level).to_numpy(dtype=float)
    return np.column_stack(tuple(columns.values())), tuple(columns)


def fit_glm(response, design, family, weights=None):
    if len(response) <= design.shape[1] or np.linalg.matrix_rank(design) < design.shape[1]:
        return None, "Insufficient rank or residual degrees of freedom"
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            fit = sm.GLM(response, design, family=family, freq_weights=weights).fit()
        if any(issubclass(entry.category, PerfectSeparationWarning) for entry in caught):
            return None, "Complete or quasi-complete separation"
    except (ValueError, np.linalg.LinAlgError, FloatingPointError, PerfectSeparationError):
        return None, "Fit or numerical convergence failure"
    if not fit.converged or not np.isfinite(fit.params).all() or not np.isfinite(fit.fittedvalues).all():
        return None, "Fit or numerical convergence failure"
    return fit, "Estimated"


@dataclass
class SelectionModel:
    population: pd.DataFrame
    design: np.ndarray
    probability: np.ndarray
    information: np.ndarray
    scores: np.ndarray


def fit_selection(primary, diagnostics):
    complete, requires_imputation = analyze_missing(primary, SELECTION_PREDICTORS, "Both centers", "Follow-up IPW",
                                                     "Follow-up logistic regression", diagnostics)
    frame = primary.loc[complete]
    if requires_imputation:
        return None, frame, "Covariate multiple imputation required by SAP"
    if frame.has_followup.nunique() < 2 or frame.site.nunique() < 2:
        return None, frame, "Both centers and both follow-up states are required"
    design, _ = design_matrix(frame, ("sex", "triage_acuity_group", "bp_phenotype", "site"))
    response = frame.has_followup.to_numpy(dtype=float)
    fitted, status = fit_glm(response, design, sm.families.Binomial())
    if fitted is None:
        return None, frame, status
    probability = np.asarray(fitted.fittedvalues)
    if not ((probability > 0) & (probability < 1)).all():
        return None, frame, "At least one fitted follow-up probability was 0 or 1"
    information = design.T @ ((probability * (1 - probability))[:, None] * design)
    scores = (response - probability)[:, None] * design
    return SelectionModel(frame, design, probability, information, scores), frame, status


def rr_confidence_interval(outcome_design, residual, mean, weight, coefficient, fit, selection=None, positions=None):
    """Patient-level HC0 sandwich; stack the selection scores when probabilities are estimated."""
    bread = outcome_design.T @ ((weight * mean)[:, None] * outcome_design)
    outcome_scores = (weight * residual)[:, None] * outcome_design
    if selection is not None:
        # The off-diagonal block is -d(outcome score)/d(propensity coefficients).
        p = selection.probability[positions]
        cross_derivative = ((weight * residual * (1 - p))[:, None] * outcome_design).T @ selection.design[positions]
        full_scores = np.zeros((len(selection.population), outcome_design.shape[1]))
        full_scores[positions] = outcome_scores
        outcome_scores = full_scores - (cross_derivative @ np.linalg.solve(selection.information, selection.scores.T)).T
    influences = np.linalg.solve(bread, outcome_scores.T)[coefficient]
    standard_error = np.sqrt(np.square(influences).sum())
    point = fit.params[coefficient]
    if not np.isfinite(standard_error) or standard_error <= 0:
        raise FloatingPointError("The sandwich variance was degenerate.")
    estimates = np.exp([point, point - cfg.WILSON_95_Z * standard_error, point + cfg.WILSON_95_Z * standard_error])
    if not (np.isfinite(estimates) & (estimates > 0)).all():
        raise FloatingPointError("The risk-ratio interval was nonfinite.")
    return tuple(float(value) for value in estimates)


def adjusted_rr(frame, evidence, site, analysis, diagnostics, selection=None, missing_selection_status=None):
    p_vs_t = frame.loc[frame.bp_phenotype.isin(("P", "T"))]
    complete, requires_imputation = analyze_missing(p_vs_t, OUTCOME_PREDICTORS, site, analysis,
                                                     "P vs T modified Poisson", diagnostics)
    used = p_vs_t.loc[complete]
    if selection is not None:
        # A patient may be outcome-model complete while missing a selection-model predictor.
        used = used.loc[used.index.intersection(selection.population.index, sort=False)]
    counts = tuple(int(used.loc[used.bp_phenotype.eq(label), evidence].sum()) for label in ("P", "T"))
    sizes = tuple(int(used.bp_phenotype.eq(label).sum()) for label in ("P", "T"))
    empty = (None, None, None)
    if missing_selection_status:
        return (*counts, *sizes, *empty, missing_selection_status)
    if requires_imputation:
        return (*counts, *sizes, *empty, "Covariate multiple imputation required by SAP")
    if not all(sizes):
        return (*counts, *sizes, *empty, "P or T absent after eligibility checks")
    if not all(counts) or sum(counts) == sum(sizes):
        return (*counts, *sizes, *empty, "Insufficient event variation or zero events in P or T")
    design, names = design_matrix(used, ("sex", "bp_phenotype"))
    coefficient = names.index("bp_phenotype=P")
    response = used[evidence].to_numpy(dtype=float)
    positions = selection.population.index.get_indexer(used.index) if selection is not None else None
    weight = 1 / selection.probability[positions] if selection is not None else np.ones(len(used))
    fitted, status = fit_glm(response, design, sm.families.Poisson(), weight)
    if fitted is None:
        return (*counts, *sizes, *empty, status)
    mean = np.asarray(fitted.fittedvalues)
    try:
        result = rr_confidence_interval(design, response - mean, mean, weight, coefficient, fitted, selection, positions)
    except (np.linalg.LinAlgError, FloatingPointError):
        return (*counts, *sizes, *empty, "Risk-ratio sandwich variance unavailable")
    return (*counts, *sizes, *result, status)


def risk_rows(eligible, observed, evidence, label, analysis):
    rows = []
    groups = [
        ("Risk of subsequent HTN evidence", phenotype, observed.bp_phenotype.eq(phenotype), True)
        for phenotype in PHENOTYPES
    ]
    groups += [
        ("Positive predictive value", "P positive", observed.bp_phenotype.eq("P"), True),
        ("Negative predictive value", "T negative", observed.bp_phenotype.eq("T"), False),
        ("Negative predictive value audit", "T or N negative", observed.bp_phenotype.isin(("T", "N")), False),
    ]
    for measure, group, in_group, match in groups:
        y = observed.loc[in_group, evidence]
        numerator, denominator = int(y.eq(match).sum()), len(y)
        rows.append((label, analysis, measure, group, len(eligible), len(observed), numerator, denominator,
                     None, None, *wilson(numerator, denominator), "Estimated" if denominator else "No observations"))
    return rows


def rr_row(eligible, observed, evidence, label, analysis, diagnostics, selection=None, missing_selection_status=None):
    n_events_p, n_events_t, n_p, n_t, rr, low, high, status = adjusted_rr(
        observed, evidence, label, analysis, diagnostics, selection, missing_selection_status,
    )
    return (label, analysis, ADJUSTED_RR, CONTRAST, len(eligible), len(observed),
            n_events_p, n_p, n_events_t, n_t, rr, low, high, status)


def ipw_rows(primary, selection, complete, status):
    rows = []
    groups = [("Both centers", "All", primary)]
    for label in SITES:
        site_frame = primary.loc[primary.site.eq(label)]
        groups.append((label, "All", site_frame))
        groups.extend((label, phenotype, site_frame.loc[site_frame.bp_phenotype.eq(phenotype)]) for phenotype in PHENOTYPES)
    for label, phenotype, eligible in groups:
        fitted = eligible.loc[eligible.index.intersection(complete.index, sort=False)]
        observed = fitted.loc[fitted.has_followup]
        probabilities = [None] * 5
        weights = [None] * 5
        ess = None
        note = status
        if selection is not None and len(fitted):
            position = selection.population.index.get_indexer(fitted.index)
            probabilities = np.quantile(selection.probability[position], (0, .05, .5, .95, 1))
            if len(observed):
                weight = 1 / selection.probability[selection.population.index.get_indexer(observed.index)]
                weights = np.quantile(weight, (0, .05, .5, .95, 1))
                ess = (weight.sum() ** 2) / np.square(weight).sum()
            if not len(observed):
                note = "No observed follow-up in this stratum"
            elif len(observed) == len(fitted):
                note = "Only observed follow-up in this stratum"
        elif not len(fitted):
            note = "No selection-model complete cases in this stratum"
        rows.append((label, phenotype, len(eligible), len(fitted), int(eligible.has_followup.sum()),
                     len(observed) if selection is not None else None, *probabilities, *weights, ess, note))
    return rows


def run(con):
    frame = load_input(con)
    primary = frame.loc[frame.index_variant.eq("primary")]
    diagnostics, main_rows, sensitivity_rows = [], [], []
    selection, complete, ipw_status = fit_selection(primary, diagnostics)
    for label, alternative, strict_label, alternative_label in (
        (SITES[0], "same_day_htn", "Follow-up starts day 7", "Count same-day HTN diagnosis"),
        (SITES[1], "no_prior_ed", "At least 2 readings, threshold 130/80", "Exclude prior ED HTN diagnosis"),
    ):
        site_main = primary.loc[primary.site.eq(label)]
        for name, population, followup, evidence in (
            ("Primary", site_main, "has_followup", "htn_evidence"),
            (strict_label, site_main, "has_followup_strict", "htn_evidence_strict"),
            ("Medication list documented", site_main.loc[site_main.has_medication_list], "has_followup", "htn_evidence"),
            (alternative_label, frame.loc[frame.site.eq(label) & frame.index_variant.eq(alternative)],
             "has_followup", "htn_evidence"),
        ):
            observed = population.loc[population[followup]]
            rows = risk_rows(population, observed, evidence, label, name)
            rows.append(rr_row(population, observed, evidence, label, name, diagnostics))
            (main_rows if name == "Primary" else sensitivity_rows).extend(rows)
        observed = site_main.loc[site_main.has_followup]
        sensitivity_rows.append(rr_row(site_main, observed, "htn_evidence", label, "Follow-up IPW", diagnostics,
                                       selection, None if selection is not None else ipw_status))

    write_rows("p3_table3.csv", ESTIMATE_HEADER, main_rows)
    write_rows("p3_sensitivity.csv", ESTIMATE_HEADER, sensitivity_rows)
    write_rows("p3_model_diagnostics.csv", (
        "site", "analysis", "model", "eligible_rows", "complete_case_rows", "covariate", "missing_n", "missing_pct", "status",
    ), diagnostics)
    write_rows("p3_ipw_diagnostics.csv", (
        "site", "phenotype", "eligible_patients", "propensity_complete_patients", "eligible_patients_with_followup",
        "weighted_patients", "probability_min", "probability_p05", "probability_median", "probability_p95",
        "probability_max", "weight_min", "weight_p05", "weight_median", "weight_p95", "weight_max",
        "effective_sample_size", "status",
    ), ipw_rows(primary, selection, complete, ipw_status))
