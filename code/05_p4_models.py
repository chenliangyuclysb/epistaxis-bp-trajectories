"""Fit the G3-confirmed GEE models and export only aggregate Aim 3 estimates."""

import csv
from dataclasses import dataclass, field
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
from scipy.stats import chi2, norm
import statsmodels.api as sm
from statsmodels.genmod.cov_struct import Independence
from statsmodels.tools.sm_exceptions import PerfectSeparationError, PerfectSeparationWarning

import config as cfg


TABLES = Path("outputs/tables")
SITES = ("MC-MED", "MIMIC-IV-ED")
LEVELS = ("T", "P", "N")
FULL_ADJUSTERS = ("age_years", "sex", "any_antithrombotic_use", "prior_htn", "triage_acuity_group")
ADJUSTMENTS = {
    "pooled": (*FULL_ADJUSTERS, "site"),
    "mimic": FULL_ADJUSTERS,
    "mcmed": ("age_years", "sex"),
    "landmark": ("age_years", "sex", "any_antithrombotic_use", "prior_htn"),
    "age": ("age_years",),
    "none": (),
}
MODEL_METHOD = "Modified Poisson GEE, independence, patient-cluster robust"
CROSSCHECK_METHOD = "Poisson GLM, patient-cluster robust"
EFFECT_FIELDS = ("estimate", "ci95_low", "ci95_high", "log_estimate", "se", "p_value", "status")
MODEL_FIELDS = (
    "eligible_visits", "visits", "patients", "events", "parameters", "converged", "fitted_mean_gt_one_n",
    "maximum_fitted_mean", "omitted_invariant_terms", "fit_status",
)


def write_csv(name, fields, rows):
    destination = TABLES / name
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(".csv.tmp")
    try:
        with temporary.open("w", encoding="utf-8-sig", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(fields)
            writer.writerows(rows)
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)


def blank_effect(status):
    return dict.fromkeys(EFFECT_FIELDS[:-1]) | {"status": status}


def effect_summary(log_effect, variance):
    if not np.isfinite(log_effect) or not np.isfinite(variance) or variance <= 0:
        return blank_effect("Nonfinite or zero robust variance")
    se = float(np.sqrt(variance))
    with np.errstate(over="ignore", under="ignore"):
        ratio, low, high = np.exp((log_effect, log_effect - cfg.WILSON_95_Z * se, log_effect + cfg.WILSON_95_Z * se))
    if not np.all(np.isfinite([ratio, low, high])) or low <= 0:
        return blank_effect("Nonfinite risk-ratio interval")
    return dict(zip(EFFECT_FIELDS, (float(ratio), float(low), float(high), float(log_effect), se,
                                     float(2 * norm.sf(abs(log_effect) / se)), "Estimated")))


@dataclass
class Fit:
    method: str = MODEL_METHOD
    n_eligible: int = 0
    n: int = 0
    patients: int = 0
    events: int = 0
    status: str = "Not fit"
    converged: bool | None = None
    n_gt_one: int | None = None
    max_mean: float | None = None
    group_counts: dict = field(default_factory=dict)
    missing_counts: dict = field(default_factory=dict)
    omitted_terms: tuple = ()
    columns: tuple = ()
    params: np.ndarray | None = None
    covariance: np.ndarray | None = None
    cross_check: "Fit | None" = None

    def estimate(self, term_or_weights):
        weights = {term_or_weights: 1.0} if isinstance(term_or_weights, str) else term_or_weights
        if self.params is None:
            return blank_effect(self.status)
        if not set(weights) <= set(self.columns):
            return blank_effect("Contrast requires a coefficient absent from the fitted model")
        contrast = np.array([weights.get(name, 0.0) for name in self.columns])
        return effect_summary(float(contrast @ self.params), float(contrast @ self.covariance @ contrast))

    def wald(self, terms):
        if self.params is None or not set(terms) <= set(self.columns):
            return None, len(terms), None, self.status if self.params is None else "Wald terms unavailable"
        positions = [self.columns.index(term) for term in terms]
        estimate = self.params[positions]
        robust_vcov = self.covariance[np.ix_(positions, positions)]
        if np.linalg.matrix_rank(robust_vcov) < len(terms):
            return None, len(terms), None, "Rank-deficient Wald covariance"
        statistic = float(estimate @ np.linalg.solve(robust_vcov, estimate))
        if not np.isfinite(statistic) or statistic < 0:
            return None, len(terms), None, "Wald statistic unavailable"
        return statistic, len(terms), float(chi2.sf(statistic, len(terms))), "Estimated"

    def diagnostics(self):
        return (self.n_eligible, self.n, self.patients, self.events, len(self.columns), self.converged,
                self.n_gt_one, self.max_mean, "|".join(self.omitted_terms), self.status)


def normalized_frame(frame, outcome, exposure, levels, continuous, adjustment, weight_column=None):
    adjusters = ADJUSTMENTS[adjustment]
    predictors = (*continuous, *((exposure,) if exposure is not None else ()))
    columns = tuple(dict.fromkeys(("site", "patient_id", "visit_id", outcome, *predictors, *adjusters,
                                  *((weight_column,) if weight_column is not None else ()))))
    sample = frame.loc[:, columns].copy()
    if exposure is not None:
        sample = sample.loc[sample[exposure].isin(levels) | sample[exposure].isna()]
    if sample[["site", "patient_id", "visit_id"]].isna().any().any() or sample.duplicated(["site", "visit_id"]).any():
        raise ValueError("Outcome models require nonmissing and unique site/visit keys and nonmissing patient keys.")
    binary = tuple(name for name in (outcome, "any_antithrombotic_use", "prior_htn") if name in sample)
    numeric = (*binary, *continuous, *(("age_years",) if "age_years" in adjusters else ()),
               *((weight_column,) if weight_column is not None else ()))
    for column in numeric:
        values = pd.to_numeric(sample[column], errors="coerce")
        sample[column] = values.where(np.isfinite(values))
    for column in binary:
        if not sample[column].dropna().isin((0, 1)).all():
            raise ValueError("Outcome and yes/no model covariates must be binary.")
    for column in ("sex", "triage_acuity_group"):
        if column in sample:
            values = sample[column].astype("string").str.strip()
            sample[column] = values.mask(values.eq(""))
    if weight_column is not None and (sample[weight_column].dropna() <= 0).any():
        raise ValueError("GEE observation weights must be positive.")
    return sample


def available_levels(values, preferred):
    observed = set(values.unique())
    return [name for name in preferred if name in observed] + sorted(observed - set(preferred))


def model_design(sample, exposure, levels, continuous, adjustment, interaction):
    adjusters = ADJUSTMENTS[adjustment]
    values = {"intercept": np.ones(len(sample))}
    omitted = []
    for name in continuous:
        values[name] = sample[name].to_numpy(dtype=float)
    if exposure is not None:
        for level in levels[1:]:
            values[f"{exposure}={level}"] = sample[exposure].eq(level).to_numpy(dtype=float)
    if "age_years" in adjusters:
        years = sample.age_years.to_numpy(dtype=float)
        center = years - years.mean()
        values["age_centered"] = center
        if adjustment in ("pooled", "mimic", "mcmed"):
            values["age_centered_squared"] = center**2
    for name in ("any_antithrombotic_use", "prior_htn"):
        if name in adjusters:
            values[name] = sample[name].to_numpy(dtype=float)
    for name, preferences in (("sex", ("Female", "Male")), ("triage_acuity_group", ("3", "1–2", "4–5")), ("site", SITES)):
        if name in adjusters:
            categories = available_levels(sample[name], preferences)
            if len(categories) < 2:
                omitted.append(name)
            for value in categories[1:]:
                values[f"{name}={value}"] = sample[name].eq(value).to_numpy(dtype=float)
    required = {*continuous, *(f"{exposure}={level}" for level in levels[1:] if exposure is not None)}
    if interaction:
        for level in levels[1:]:
            main_term = f"{exposure}={level}"
            center_term = f"site={SITES[1]}"
            term = f"{main_term}:{center_term}"
            values[term] = values[main_term] * values[center_term]
            required.add(term)
    for name in tuple(values):
        if name != "intercept" and name not in required and not np.ptp(values[name]):
            omitted.append(name)
            del values[name]
    return np.column_stack(tuple(values.values())), tuple(values), tuple(omitted)


def finish_fit(summary, result, caught):
    summary.converged = bool(getattr(result, "converged", True))
    fitted = np.asarray(result.fittedvalues, dtype=float)
    finite = np.isfinite(fitted).all()
    summary.n_gt_one = int(np.count_nonzero(fitted > 1)) if finite else None
    summary.max_mean = float(fitted.max()) if finite else None
    covariance = np.asarray(result.cov_params())
    estimates = np.asarray(result.params)
    if any(issubclass(warning.category, PerfectSeparationWarning) for warning in caught):
        summary.status = "Separation warning"
    elif not summary.converged:
        summary.status = "Model did not converge"
    elif not finite or not np.isfinite(covariance).all() or not np.isfinite(estimates).all():
        summary.status = "Nonfinite fit or cluster-robust covariance"
    else:
        summary.status = "Estimated"
        summary.params = estimates
        summary.covariance = covariance
    return summary


def poisson_check(summary, response, design, clusters, weights):
    check = Fit(method=CROSSCHECK_METHOD, n_eligible=summary.n_eligible, n=summary.n, patients=summary.patients,
                events=summary.events, group_counts=summary.group_counts, missing_counts=summary.missing_counts,
                omitted_terms=summary.omitted_terms, columns=summary.columns)
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            model = sm.GLM(response, design, family=sm.families.Poisson(), var_weights=weights)
            result = model.fit(cov_type="cluster", cov_kwds={"groups": clusters, "use_correction": False})
        return finish_fit(check, result, caught)
    except (ValueError, np.linalg.LinAlgError, FloatingPointError, PerfectSeparationError):
        check.status = "Poisson GLM cross-check did not produce a valid fit"
        check.converged = False
        return check


def fit_model(frame, outcome, *, exposure="bp_phenotype", levels=LEVELS, continuous=(), adjustment="pooled",
              interaction=False, weight_column=None):
    if interaction and (adjustment != "pooled" or exposure is None):
        raise ValueError("The phenotype-by-center interaction requires the pooled categorical model.")
    sample = normalized_frame(frame, outcome, exposure, levels, continuous, adjustment, weight_column)
    used = sample.dropna()
    counts = {}
    if exposure is not None:
        for level in levels:
            group = used.loc[used[exposure].eq(level)]
            counts[level] = (len(group), int(group[outcome].sum()))
    summary = Fit(n_eligible=len(sample), n=len(used), patients=len(used[["site", "patient_id"]].drop_duplicates()),
                  events=int(used[outcome].sum()), group_counts=counts, missing_counts=sample.isna().sum().to_dict())
    if not len(sample):
        summary.status = "No eligible visits"
        return summary
    if any(100 * summary.missing_counts[name] / len(sample) > cfg.P4_COVARIATE_MI_TRIGGER_PERCENT
           for name in ADJUSTMENTS[adjustment]):
        summary.status = "Covariate multiple imputation required by SAP"
        return summary
    if exposure is not None and not all(n for n, _ in counts.values()):
        summary.status = "At least one prespecified exposure group is absent"
        return summary
    if summary.events in (0, summary.n) or any(not events for _, events in counts.values()):
        summary.status = "No event variation or zero events in at least one exposure group"
        return summary
    if adjustment == "pooled" and set(used.site) != set(SITES):
        summary.status = "The pooled model requires both centers"
        return summary
    design, columns, omitted = model_design(used, exposure, levels, continuous, adjustment, interaction)
    if weight_column is not None:
        summary.method += ", fixed participation weights"
    summary.columns = columns
    summary.omitted_terms = omitted
    if summary.patients <= len(columns) or np.linalg.matrix_rank(design) < len(columns):
        summary.status = "Insufficient patient clusters or singular design"
        return summary
    clusters = pd.MultiIndex.from_frame(used[["site", "patient_id"]]).factorize()[0]
    response = used[outcome].to_numpy(dtype=float)
    weights = used[weight_column].to_numpy(dtype=float) if weight_column is not None else None
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            model = sm.GEE(response, design, groups=clusters, family=sm.families.Poisson(),
                           cov_struct=Independence(), weights=weights)
            result = model.fit(cov_type="robust")
        finish_fit(summary, result, caught)
    except (ValueError, np.linalg.LinAlgError, FloatingPointError, PerfectSeparationError):
        summary.status = "GEE did not produce a valid fit"
        summary.converged = False
    if summary.params is None or summary.n_gt_one:
        summary.cross_check = poisson_check(summary, response, design, clusters, weights)
    return summary


ESTIMATE_FIELDS = (
    "endpoint", "analysis", "site", "contrast", "method", "exposed_visits", "exposed_events", "reference_visits",
    "reference_events", *EFFECT_FIELDS, *MODEL_FIELDS,
)


def result_row(endpoint, analysis, site, contrast, fit, exposure_group, reference_group, effects):
    return (endpoint, analysis, site, contrast, fit.method,
            *fit.group_counts.get(exposure_group, (None, None)), *fit.group_counts.get(reference_group, (None, None)),
            *(effects[field] for field in EFFECT_FIELDS), *fit.diagnostics())


def add_phenotype_rows(rows, endpoint, analysis, site, fit, levels=LEVELS, exposure="bp_phenotype"):
    reference = levels[0]
    for group in levels[1:]:
        effect = fit.estimate(f"{exposure}={group}")
        if not (endpoint == "7-day return" and analysis == "Main GEE" and site == "Pooled" and group == "P"):
            effect["p_value"] = None
        rows.append(result_row(endpoint, analysis, site, f"{group} vs {reference}", fit, group, reference, effect))


def e_value(risk_ratio):
    if risk_ratio is None:
        return None
    toward_harm = max(risk_ratio, 1 / risk_ratio)
    return float(toward_harm + np.sqrt(toward_harm * (toward_harm - 1)))


def meta_row(site_fits):
    effects = [fit.estimate("bp_phenotype=P") for fit in site_fits]
    fitted_meta = Fit(method="Fixed-effect inverse-variance meta-analysis", n_eligible=sum(fit.n_eligible for fit in site_fits),
                      n=sum(fit.n for fit in site_fits), patients=sum(fit.patients for fit in site_fits),
                      events=sum(fit.events for fit in site_fits),
                      group_counts={level: tuple(sum(fit.group_counts.get(level, (0, 0))[i] for fit in site_fits)
                                                for i in (0, 1)) for level in LEVELS})
    if any(effect["se"] is None for effect in effects):
        fitted_meta.status = "Both site-specific P vs T estimates required"
        combined = blank_effect(fitted_meta.status)
    else:
        precisions = np.array([1 / effect["se"]**2 for effect in effects])
        logs = np.array([effect["log_estimate"] for effect in effects])
        combined = effect_summary(float(np.average(logs, weights=precisions)), float(1 / precisions.sum()))
        fitted_meta.status = combined["status"]
        combined["p_value"] = None
    return result_row("7-day return", "Fixed-effect meta-analysis", "Pooled", "P vs T", fitted_meta, "P", "T", combined)


def collect_missingness(rows, model, endpoint, analysis, site):
    for name, missing in model.missing_counts.items():
        if name not in ("site", "visit_id", "patient_id"):
            rows.append((endpoint, analysis, site, name, model.n_eligible, missing,
                         100 * missing / model.n_eligible if model.n_eligible else None))


def run_primary_models(frame):
    estimates, crosschecks, groups, missing, e_values = [], [], [], [], []
    primary_site_fits = []
    for outcome, endpoint in (("return_7d_event", "7-day return"), ("return_30d_event", "30-day return")):
        for site, adjustment in (("Pooled", "pooled"), (SITES[0], "mcmed"), (SITES[1], "mimic")):
            subset = frame if site == "Pooled" else frame.loc[frame.site.eq(site)]
            fit = fit_model(subset, outcome, adjustment=adjustment)
            add_phenotype_rows(estimates, endpoint, "Main GEE", site, fit)
            if fit.cross_check is not None:
                reason = "GEE did not yield an estimate" if fit.params is None else "Fitted GEE mean >1"
                add_phenotype_rows(crosschecks, endpoint, reason, site, fit.cross_check)
            collect_missingness(missing, fit, endpoint, "Main GEE", site)
            for level in LEVELS:
                eligible = subset.loc[subset.bp_phenotype.eq(level)]
                observed = eligible[outcome].dropna()
                n, events = fit.group_counts.get(level, (0, 0))
                patients_n = len(eligible[["site", "patient_id"]].drop_duplicates())
                groups.append((endpoint, site, level, len(eligible), patients_n, len(observed), int(observed.sum()), n, events))
            if outcome == "return_7d_event":
                if site != "Pooled":
                    primary_site_fits.append(fit)
                else:
                    r = fit.estimate("bp_phenotype=P")
                    if r["estimate"] is None:
                        point_e, bound_e = None, None
                    else:
                        point_e = e_value(r["estimate"])
                        low, high = r["ci95_low"], r["ci95_high"]
                        bound_e = 1.0 if low <= 1 <= high else e_value(low if low > 1 else high)
                    e_values.append((site, r["estimate"], r["ci95_low"], r["ci95_high"], point_e, bound_e, r["status"]))
    estimates.append(meta_row(primary_site_fits))

    interaction = fit_model(frame, "return_7d_event", interaction=True)
    interaction_terms = tuple(f"bp_phenotype={level}:site={SITES[1]}" for level in LEVELS[1:])
    stat, df, p, int_status = interaction.wald(interaction_terms)
    interaction_rows = [("7-day return", "Phenotype × center", interaction.method, stat, df, p, int_status,
                         *interaction.diagnostics())]
    if interaction.cross_check is not None:
        check = interaction.cross_check
        check_stat, check_df, check_p, check_status = check.wald(interaction_terms)
        interaction_rows.append(("7-day return", "Phenotype × center cross-check", check.method, check_stat, check_df,
                                 check_p, check_status, *check.diagnostics()))
    collect_missingness(missing, interaction, "7-day return", "Phenotype × center", "Pooled")
    write_csv("p4_main_estimates.csv", ESTIMATE_FIELDS, estimates)
    write_csv("p4_main_crosscheck.csv", ESTIMATE_FIELDS, crosschecks)
    write_csv("p4_main_group_counts.csv", ("endpoint", "site", "phenotype", "eligible_visits", "patients", "observed_outcomes",
                                          "observed_events", "model_visits", "model_events"), groups)
    write_csv("p4_main_missingness.csv", ("endpoint", "analysis", "site", "variable", "eligible_visits", "missing_visits",
                                         "missing_pct"), missing)
    write_csv("p4_main_interaction.csv", ("endpoint", "test", "method", "wald_chi2", "df", "p_value", "status", *MODEL_FIELDS),
              interaction_rows)
    write_csv("p4_main_evalue.csv", ("site", "estimate", "ci95_low", "ci95_high", "point_e_value", "ci_limit_e_value", "status"), e_values)


def rcs_basis(sbp, knots, unit):
    """Three-knot restricted cubic spline, centered at the median knot and scaled per SBP unit."""
    k1, k2, k3 = (np.asarray(knots, dtype=float) - knots[1]) / unit
    scaled = (np.asarray(sbp, dtype=float) - knots[1]) / unit
    if not k1 < k2 < k3:
        raise ValueError("Spline knot quantiles must be distinct.")

    def nonlinear(x):
        return (np.maximum(x - k1, 0)**3 - (k3 - k1) / (k3 - k2) * np.maximum(x - k2, 0)**3
                + (k2 - k1) / (k3 - k2) * np.maximum(x - k3, 0)**3) / (k3 - k1)**2

    return np.column_stack((scaled, nonlinear(scaled) - nonlinear(0)))


def run_sbp_model(frame):
    outcome = "return_7d_event"
    terms = ("last_sbp_linear_per_10", "last_sbp_nonlinear_basis")
    eligible = normalized_frame(frame, outcome, None, (), ("sbp_last",), "pooled").dropna()
    if len(eligible):
        knots = np.quantile(eligible.sbp_last, cfg.P4_SBP_KNOT_QUANTILES)
        display = np.quantile(eligible.sbp_last, cfg.P4_SBP_CURVE_QUANTILES)
    else:
        knots = np.full(len(cfg.P4_SBP_KNOT_QUANTILES), np.nan)
        display = np.full(len(cfg.P4_SBP_CURVE_QUANTILES), np.nan)
    unit = cfg.P4_SBP_UNIT
    model_frame = frame.copy()
    if np.isfinite(knots).all() and np.all(np.diff(knots) > 0):
        raw = pd.to_numeric(model_frame.sbp_last, errors="coerce")
        model_frame[list(terms)] = rcs_basis(raw.to_numpy(dtype=float, na_value=np.nan), knots, unit)
        model = fit_model(model_frame, outcome, exposure=None, levels=(), continuous=terms)
    else:
        model = Fit(n_eligible=len(frame), n=len(eligible), events=int(eligible[outcome].sum()),
                    patients=len(eligible[["site", "patient_id"]].drop_duplicates()), status="Distinct spline knots unavailable")
    summaries, tests = [], []
    for fit in (model, model.cross_check):
        if fit is None:
            continue
        for term in fit.columns:
            effect = fit.estimate(term)
            summaries.append((fit.method, term, *(effect[field] for field in EFFECT_FIELDS if field != "p_value"), *fit.diagnostics()))
        if np.isfinite(knots).all() and np.all(np.diff(knots) > 0):
            dose = rcs_basis([knots[1] + unit], knots, unit)[0]
            effect = fit.estimate(dict(zip(terms, dose)))
        else:
            effect = blank_effect(fit.status)
        summaries.append((fit.method, f"Median SBP plus {unit} mmHg vs median SBP",
                          *(effect[field] for field in EFFECT_FIELDS if field != "p_value"), *fit.diagnostics()))
        for label, tested in (("Overall SBP spline", terms), ("Nonlinear SBP component", (terms[1],))):
            tests.append((fit.method, label, *fit.wald(tested), *fit.diagnostics()))

    curve = []
    if np.isfinite(display).all() and model.params is not None:
        grid = np.unique(np.r_[np.linspace(*display, num=cfg.P4_SBP_CURVE_POINTS), knots[1]])
        for x, values in zip(grid, rcs_basis(grid, knots, unit)):
            effect = model.estimate(dict(zip(terms, values))) if x != knots[1] else {"estimate": 1, "ci95_low": 1, "ci95_high": 1}
            curve.append((x, knots[1], effect["estimate"], effect["ci95_low"], effect["ci95_high"]))
    write_csv("p4_sbp_estimates.csv", ("method", "term_or_contrast", *(field for field in EFFECT_FIELDS if field != "p_value"),
                                       *MODEL_FIELDS), summaries)
    write_csv("p4_sbp_tests.csv", ("method", "test", "wald_chi2", "df", "p_value", "status", *MODEL_FIELDS), tests)
    write_csv("p4_sbp_curve.csv", ("sbp_last_mmhg", "reference_sbp_mmhg", "estimate", "ci95_low", "ci95_high"), curve)
    write_csv("p4_sbp_knots.csv", ("lower_knot_mmhg", "median_knot_reference_mmhg", "upper_knot_mmhg", "display_min_mmhg",
                                    "display_max_mmhg", "sbp_unit_mmhg", "model_visits", "model_status"),
              [(*knots, *display, unit, model.n, model.status)])


def run(frame):
    if "is_c2" in frame and not frame.is_c2.fillna(False).all():
        raise ValueError("Primary Aim 3 models require the selected C2 analysis frame.")
    run_primary_models(frame)
    run_sbp_model(frame)
