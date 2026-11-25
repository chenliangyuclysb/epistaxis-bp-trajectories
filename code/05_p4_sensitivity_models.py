"""Fit the G3 Aim 3 sensitivity scenarios and the SAP C1 selection analysis."""

from importlib import import_module
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.tools.sm_exceptions import PerfectSeparationWarning

import config as cfg


TABLES = Path("outputs/tables")
SITES = ("MC-MED", "MIMIC-IV-ED")
PHENOTYPES = ("N", "T", "P")
PREDICTORS = (
    "age_years", "sex", "any_antithrombotic_use", "prior_htn", "triage_acuity_group", "site",
)
SCENARIO_TITLES = {
    "S11-1": "BP threshold 160/100",
    "S11-2": "First eligible visit per patient",
    "S11-3": "Primary diagnosis only",
    "S11-4": "Exclude documented antithrombotic use",
    "S11-5": "Exclude HHT",
    "S11-6": "Exclude acute class A treatment",
    "S11-7": "Include LWBS, AMA and Observation",
    "S11-9": "Last BP at least 30 minutes before departure",
    "S11-12": "MC-MED medication list required",
    "S11-16": "Medication list required in both centers",
    "S11-13": "MC-MED count same-day HTN code as prior",
    "S11-17": "MIMIC omit prior ED HTN source",
    "S10-IPW": "C1 selection IPW",
}
ONE_SITE_CHANGES = {
    "S11-12": "MC-MED",
    "S11-13": "MC-MED",
    "S11-17": "MIMIC-IV-ED",
}
ESTIMATE_FIELDS = (
    "scenario", "analysis", "site", "outcome", "contrast", "input_visits",
    "model_visits", "model_patients", "model_events", "n_visits", "n_events",
    "t_visits", "t_events", "p_visits", "p_events", "rr", "ci95_low", "ci95_high",
    "p_value", "converged", "fitted_mean_gt_one", "status", "omitted_terms",
    "crosscheck_status", "crosscheck_rr", "crosscheck_ci95_low", "crosscheck_ci95_high",
)


def direct_home_candidates(index):
    at_home = (
        (index.site.eq("MC-MED") & index.disposition.eq(cfg.DIRECT_HOME_MCMED))
        | (index.site.eq("MIMIC-IV-ED") & index.disposition.eq(cfg.DIRECT_HOME_MIMIC))
    )
    home = index.loc[at_home & ~index.incomplete_visit].copy()
    if int(home.is_c1.sum()) != int(index.is_c2.sum()):
        raise ValueError("C1-availability target does not reproduce C2.")
    return home


def describe_c1(index, write_csv):
    wilson = import_module("02_p1_tables").wilson
    continuous = ("age_years",)
    categories = (
        "sex", "triage_acuity_group", "any_antithrombotic_use",
        "prior_htn", "hht", "has_medication_list",
    )
    rows = []
    for target, population in (("All new index visits", index),
                               ("Direct-home new index visits", direct_home_candidates(index))):
        for site, subset in population.groupby("site", sort=True):
            for c1 in (False, True):
                group = subset.loc[subset.is_c1.eq(c1)]
                rows.append((target, site, c1, "C1 participation", "", len(group),
                             group.patient_id.nunique(), *([None] * 9)))
                for variable in continuous:
                    data = group[variable].dropna()
                    q = data.quantile((0.25, 0.5, 0.75)).tolist() if len(data) else [None] * 3
                    rows.append((target, site, c1, variable, "", len(group), group.patient_id.nunique(),
                                 len(data), len(group) - len(data), None, None, None, None, *q))
                for variable in categories:
                    data = group[variable].dropna()
                    for level, count in data.value_counts().sort_index().items():
                        interval = tuple(bound / 100 for bound in wilson(int(count), len(data)))
                        rows.append((target, site, c1, variable, str(level), len(group), group.patient_id.nunique(),
                                     len(data), len(group) - len(data), int(count),
                                     float(count / len(data)), *interval, None, None, None))
    write_csv(TABLES / "p4_c1_comparison.csv",
              ("target", "site", "is_c1", "characteristic", "level", "visits", "patients",
               "observed", "missing", "count", "proportion", "wilson95_low", "wilson95_high",
               "q1", "median", "q3"), rows)


def selection_design(complete):
    age = complete.age_years.to_numpy(dtype=float)
    centered = (age - age.mean()) / age.std()
    dummies = pd.get_dummies(
        complete[["sex", "triage_acuity_group", "site"]].astype("string"),
        drop_first=True, dtype=float,
    )
    return np.column_stack((
        np.ones(len(complete)), centered, centered**2,
        complete.any_antithrombotic_use.to_numpy(dtype=float),
        complete.prior_htn.to_numpy(dtype=float), dummies.to_numpy(),
    ))


def c1_weights(index, write_csv):
    frame = direct_home_candidates(index)
    complete = frame.dropna(subset=list(PREDICTORS)).copy()
    missing_percent = (100 * frame[list(PREDICTORS)].isna().mean()).to_dict()
    status = "Estimated"
    probability = np.full(len(complete), np.nan)
    fit = None
    if any(value > cfg.P4_COVARIATE_MI_TRIGGER_PERCENT for value in missing_percent.values()):
        status = "Selection model requires multiple imputation per S10"
    elif not len(complete) or complete.is_c1.nunique() != 2:
        status = "Insufficient C1 participation variation"
    else:
        design = selection_design(complete)
        if not np.isfinite(design).all() or np.linalg.matrix_rank(design) != design.shape[1]:
            status = "Selection model rank failure"
        else:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                fit = sm.GLM(complete.is_c1.astype(int), design, family=sm.families.Binomial()).fit()
            probability = np.asarray(fit.fittedvalues)
            separated = any(issubclass(entry.category, PerfectSeparationWarning) for entry in caught)
            if separated or not fit.converged or not ((probability > 0) & (probability < 1)).all():
                status = "Selection model convergence or positivity failure"
    if status != "Estimated":
        probability[:] = np.nan
    complete["c1_probability"] = probability
    observed = complete.loc[complete.is_c1].copy()
    observed["c1_ipw"] = 1.0 / observed.c1_probability

    rows = []
    for site in ("Both centers", *SITES):
        candidate = frame if site == "Both centers" else frame.loc[frame.site.eq(site)]
        modeled = complete if site == "Both centers" else complete.loc[complete.site.eq(site)]
        c1 = observed if site == "Both centers" else observed.loc[observed.site.eq(site)]
        p = modeled.c1_probability.dropna().to_numpy()
        w = c1.c1_ipw.dropna().to_numpy()
        p_q = np.quantile(p, (0, 0.05, 0.5, 0.95, 1)).tolist() if len(p) else [None] * 5
        w_q = np.quantile(w, (0, 0.05, 0.5, 0.95, 1)).tolist() if len(w) else [None] * 5
        ess = float(w.sum() ** 2 / np.square(w).sum()) if len(w) else None
        rows.append((site, len(candidate), candidate.patient_id.nunique() if site != "Both centers" else
                     candidate[["site", "patient_id"]].drop_duplicates().shape[0],
                     len(modeled), int(candidate.is_c1.sum()), len(w), *p_q, *w_q, ess, status))
    write_csv(TABLES / "p4_c1_ipw_diagnostics.csv",
              ("site", "direct_home_index_visits", "direct_home_index_patients", "complete_propensity_visits",
               "c1_visits", "weighted_visits", "probability_min", "probability_p05", "probability_median",
               "probability_p95", "probability_max", "weight_min", "weight_p05", "weight_median",
               "weight_p95", "weight_max", "effective_sample_size", "status"), rows)
    return observed.loc[observed.c1_ipw.notna()].copy()


def sensitivity_effects(cohorts, models):
    rows = []
    for scenario in SCENARIO_TITLES:
        if scenario not in cohorts:
            raise ValueError("A prespecified sensitivity analysis is missing.")
        population = cohorts[scenario]
        if scenario in ONE_SITE_CHANGES:
            scopes = ("Both centers", ONE_SITE_CHANGES[scenario])
        else:
            scopes = ("Both centers", *SITES)
        for scope in scopes:
            subset = population if scope == "Both centers" else population.loc[population.site.eq(scope)]
            adjustment = {"Both centers": "pooled", "MC-MED": "mcmed", "MIMIC-IV-ED": "mimic"}[scope]
            kwargs = {"weight_column": "c1_ipw"} if scenario == "S10-IPW" else {}
            fitted = models.fit_model(subset, "return_7d_event", adjustment=adjustment, **kwargs)
            groups = [fitted.group_counts.get(level, (0, 0)) for level in PHENOTYPES]
            for contrast, term in (("P vs T", "bp_phenotype=P"), ("N vs T", "bp_phenotype=N")):
                estimate = fitted.estimate(term)
                check = (fitted.cross_check.estimate(term) if fitted.cross_check else
                         {"status": "Not triggered", "estimate": None, "ci95_low": None, "ci95_high": None})
                rows.append((scenario, SCENARIO_TITLES[scenario], scope, "7-day epistaxis return",
                             contrast, len(subset), fitted.n, fitted.patients, fitted.events,
                             *(value for pair in groups for value in pair),
                             estimate["estimate"], estimate["ci95_low"], estimate["ci95_high"],
                             estimate["p_value"], fitted.converged, fitted.n_gt_one, estimate["status"],
                             " | ".join(fitted.omitted_terms), check["status"], check["estimate"],
                             check["ci95_low"], check["ci95_high"]))
    return rows


def run(index, cohorts):
    p0 = import_module("01_p0_audit")
    models = import_module("05_p4_models")
    describe_c1(index, p0.write_csv)
    weighted = c1_weights(index, p0.write_csv)
    all_cohorts = dict(cohorts, **{"S10-IPW": weighted})
    rows = sensitivity_effects(all_cohorts, models)
    p0.write_csv(TABLES / "p4_sensitivity_estimates.csv", ESTIMATE_FIELDS, rows)
