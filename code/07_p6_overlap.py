"""Fit the exploratory B-class overlap-weighted return model approved at G5."""

from importlib import import_module
import sys
import warnings

import duckdb
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.tools.sm_exceptions import PerfectSeparationWarning

import config as cfg


SITES = (("mcmed", "MC-MED"), ("mimic", "MIMIC-IV-ED"))
TREATMENT = "ed_any_aht"
OUTCOME = "return_7d_event"
PROPENSITY_FIELDS = ("age_years", "sex", "sbp_first", "prior_htn", "any_antithrombotic_use", "triage_acuity_group", "site")
COHORT_COLUMNS = (
    "visit_id", "patient_id", "age_years", "sex", "prior_htn", "any_antithrombotic_use", "triage_acuity_group",
    "ed_acute_aht", TREATMENT, "ed_alpha1_sensitivity", OUTCOME,
)
BALANCE_LEVELS = {
    "sex": ("Female", "Male"),
    "triage_acuity_group": tuple(label for label, _ in cfg.TRIAGE_ACUITY_GROUPS),
    "site": tuple(site for _, site in SITES),
}


def selected_population(con):
    parts = [f"""
        SELECT '{site}' AS site, {', '.join('c.' + field for field in COHORT_COLUMNS)}, p.sbp_first
        FROM read_parquet('derived/cohort_{short}.parquet') c
        LEFT JOIN read_parquet('derived/phenotypes_{short}.parquet') p USING (visit_id)
        WHERE c.is_c2
    """ for short, site in SITES]
    return con.execute(" UNION ALL ".join(parts) + " ORDER BY site, patient_id, visit_id").df()


def complete_cases(frame):
    frame = frame.copy()
    if frame.empty or frame.duplicated(["site", "visit_id"]).any() or frame[["visit_id", "patient_id"]].isna().any().any():
        raise ValueError("C2 requires nonempty, unique site/visit keys and nonmissing patient keys.")
    for name in (TREATMENT, OUTCOME, "ed_acute_aht", "ed_alpha1_sensitivity"):
        if not frame[name].isin((True, False)).all():
            raise ValueError("The C2 exposure and outcome flags must be known and binary.")
    if (frame.ed_acute_aht & ~frame[TREATMENT]).any() or (frame.ed_alpha1_sensitivity & ~frame[TREATMENT]).any():
        raise ValueError("B-class C2 membership differs from the expected nested/alpha-1 classification.")
    for name in ("age_years", "sbp_first"):
        values = pd.to_numeric(frame[name], errors="coerce")
        frame[name] = values.where(np.isfinite(values))
    for name, labels in BALANCE_LEVELS.items():
        frame[name] = frame[name].astype("string").str.strip()
        frame[name] = frame[name].mask(frame[name].eq(""))
        if not frame[name].dropna().isin(labels).all():
            raise ValueError("The C2 data include unrecognized categorical covariate levels.")
    for name in ("prior_htn", "any_antithrombotic_use"):
        if not frame[name].dropna().isin((True, False)).all():
            raise ValueError("Historical hypertension and antithrombotic use must be binary covariates.")
    if int(frame.loc[frame[TREATMENT], OUTCOME].sum()) < cfg.P5_RARE_TREATED_EVENT_COUNT:
        raise ValueError("The B-class treatment event count did not meet the prespecified weighting threshold.")
    counts = frame[list(PROPENSITY_FIELDS)].isna().sum()
    observed = frame.dropna(subset=list(PROPENSITY_FIELDS)).copy()
    rows = [(name, len(frame), len(observed), int(missing), 100 * missing / len(frame)) for name, missing in counts.items()]
    if any(row[-1] > cfg.P4_COVARIATE_MI_TRIGGER_PERCENT for row in rows):
        raise ValueError("A required B-class propensity covariate met the SAP multiple-imputation criterion.")
    if int(observed.loc[observed[TREATMENT], OUTCOME].sum()) < cfg.P5_RARE_TREATED_EVENT_COUNT:
        raise ValueError("The B-class event count fell below the weighting threshold after eligibility checks.")
    return observed, rows


def fit_propensity(frame):
    design = {
        "intercept": np.ones(len(frame)),
        "age_years": frame.age_years.to_numpy(dtype=float) - frame.age_years.mean(),
        "sbp_first": frame.sbp_first.to_numpy(dtype=float) - frame.sbp_first.mean(),
        "prior_htn": frame.prior_htn.to_numpy(dtype=float),
        "any_antithrombotic_use": frame.any_antithrombotic_use.to_numpy(dtype=float),
    }
    for variable, levels in BALANCE_LEVELS.items():
        for level in levels[1:]:
            design[f"{variable}={level}"] = frame[variable].eq(level).to_numpy(dtype=float)
    x = np.column_stack(tuple(design.values()))
    if np.linalg.matrix_rank(x) < x.shape[1] or frame[TREATMENT].nunique() != 2:
        raise ValueError("Insufficient rank or B-class exposure variation for the propensity model.")
    treatment = frame[TREATMENT].to_numpy(dtype=float)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        fitted = sm.GLM(treatment, x, family=sm.families.Binomial()).fit()
    p = np.asarray(fitted.fittedvalues)
    if (not fitted.converged or any(issubclass(item.category, PerfectSeparationWarning) for item in caught)
            or not np.isfinite(p).all() or ((p <= 0) | (p >= 1)).any()):
        raise ValueError("The planned logistic propensity model did not yield finite interior probabilities.")
    return frame.assign(propensity=p, overlap_weight=np.where(treatment == 1, 1 - p, p), b_class_exposure=treatment)


def balance(weighted):
    treated = weighted[TREATMENT].to_numpy(dtype=bool)
    weight = weighted.overlap_weight.to_numpy()
    columns = [(name, "", weighted[name].to_numpy(dtype=float)) for name in ("age_years", "sbp_first")]
    columns += [(name, "Yes", weighted[name].to_numpy(dtype=float)) for name in ("prior_htn", "any_antithrombotic_use")]
    columns += [(name, level, weighted[name].eq(level).to_numpy(dtype=float))
                for name, levels in BALANCE_LEVELS.items() for level in levels]
    rows = []
    for name, level, values in columns:
        a, b = values[treated], values[~treated]
        means = (a.mean(), b.mean(), np.average(a, weights=weight[treated]), np.average(b, weights=weight[~treated]))
        sd = np.sqrt((a.var() + b.var()) / 2)
        unweighted_smd = (means[0] - means[1]) / sd if sd > 0 else None
        weighted_smd = (means[2] - means[3]) / sd if sd > 0 else None
        rows.append((name, level, *means, sd, unweighted_smd, weighted_smd))
    return rows


def weight_diagnostics(weighted):
    rows = []
    scopes = [("Pooled", weighted), *((site, weighted.loc[weighted.site.eq(site)]) for _, site in SITES)]
    for site, stratum in scopes:
        for arm, group in (("All", stratum), ("B class", stratum.loc[stratum[TREATMENT]]),
                           ("No B class", stratum.loc[~stratum[TREATMENT]])):
            patient_count = len(group[["site", "patient_id"]].drop_duplicates())
            if len(group):
                weight = group.overlap_weight.to_numpy()
                propensity = np.quantile(group.propensity, (0, 0.05, 0.5, 0.95, 1))
                distribution = np.quantile(weight, (0, 0.05, 0.5, 0.95, 1))
                total = float(weight.sum())
                ess = total**2 / np.square(weight).sum()
            else:
                propensity, distribution, total, ess = [None] * 5, [None] * 5, 0, None
            rows.append((site, arm, len(group), patient_count, int(group[OUTCOME].sum()), total, ess, *propensity, *distribution))
    return rows


def tables_for(frame, models):
    observed, missingness = complete_cases(frame)
    weighted = fit_propensity(observed)
    fitted = models.fit_model(weighted, OUTCOME, exposure=None, levels=(), continuous=("b_class_exposure",),
                              adjustment="none", weight_column="overlap_weight")
    effect = fitted.estimate("b_class_exposure")
    if effect["status"] != "Estimated" or fitted.n != len(weighted):
        raise RuntimeError("The planned weighted clustered modified Poisson model was not estimable on the full eligible sample.")
    exposed = weighted.loc[weighted[TREATMENT]]
    unexposed = weighted.loc[~weighted[TREATMENT]]
    risks = tuple(np.average(group[OUTCOME].to_numpy(dtype=float), weights=group.overlap_weight) for group in (exposed, unexposed))
    if not np.isclose(risks[0] / risks[1], effect["estimate"]):
        raise RuntimeError("The treatment-only weighted risk ratio did not reconcile to the two weighted risks.")
    estimate = ("Pooled C2", "B class vs no B class", fitted.n, fitted.patients, len(exposed), int(exposed[OUTCOME].sum()),
                len(unexposed), int(unexposed[OUTCOME].sum()), *(100 * risk for risk in risks),
                *(effect[key] for key in ("estimate", "ci95_low", "ci95_high", "log_estimate", "se")),
                f"{models.MODEL_METHOD}, overlap weights treated as fixed", True, fitted.converged, fitted.n_gt_one,
                effect["status"])
    return {
        "p6_bclass_overlap_estimate.csv": ((
            "population", "treatment", "visits", "patients", "treated_visits", "treated_events", "untreated_visits", "untreated_events",
            "weighted_treated_risk_percent", "weighted_untreated_risk_percent", "estimate", "ci95_low", "ci95_high", "log_estimate", "se",
            "method", "propensity_converged", "outcome_converged", "outcome_fitted_mean_gt_one_n", "status",
        ), [estimate]),
        "p6_bclass_overlap_balance.csv": ((
            "covariate", "level", "unweighted_treated_mean", "unweighted_untreated_mean", "weighted_treated_mean", "weighted_untreated_mean",
            "unweighted_pooled_sd", "unweighted_smd", "weighted_smd",
        ), balance(weighted)),
        "p6_bclass_overlap_weights.csv": ((
            "site", "treatment_arm", "visits", "patients", "events", "sum_weights", "effective_sample_size",
            "propensity_min", "propensity_p05", "propensity_median", "propensity_p95", "propensity_max",
            "weight_min", "weight_p05", "weight_median", "weight_p95", "weight_max",
        ), weight_diagnostics(weighted)),
        "p6_bclass_overlap_missingness.csv": ((
            "covariate", "eligible_visits", "modeled_visits", "missing_visits", "missing_pct",
        ), missingness),
    }


def run(con):
    models = import_module("05_p4_models")
    for name, (header, rows) in tables_for(selected_population(con), models).items():
        models.write_csv(name, header, rows)


def main():
    try:
        import_module("02_p1_setup").verify_project_runtime()
        with duckdb.connect() as con:
            run(con)
    except Exception as error:
        # Exceptions can contain source data and are never formatted as full messages.
        print(f"P6 B-class overlap weighting stopped ({type(error).__name__}).", file=sys.stderr)
        return 1
    print("P6 B-class overlap weighting completed; aggregate results written.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
