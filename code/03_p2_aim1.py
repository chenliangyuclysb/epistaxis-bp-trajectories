"""Estimate Aim 1 from the outcome-free C1 phenotype files."""

import math
import warnings
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import t as student_t
from statsmodels.tools.sm_exceptions import PerfectSeparationError

import config as cfg


ROOT = Path(__file__).resolve().parents[1]
SITES = (("mcmed", "MC-MED"), ("mimic", "MIMIC-IV-ED"))
PHENOTYPES = ("N", "T", "P")
PRIMARY = "Primary"
HTN_GROUPS = (False, True)
INPUT_FIELDS = (
    "visit_id", "patient_id", "bp_phenotype",
    "sbp_first", "sbp_last", "time_first_minutes", "time_last_minutes",
    "first_normal_minutes", "arrival_elevated", "arrival_severe", "last_normal",
    "prior_htn", "has_medication_list", "prior_htn_medication",
    "ed_acute_aht", "ed_first_acute_aht_minutes",
)


def load_site(con, derived, short):
    fields = INPUT_FIELDS
    if short == "mimic":
        fields += ("prior_htn_hosp_diagnosis", "prior_htn_ed_diagnosis")
    path = derived / f"phenotypes_{short}.parquet"
    return con.execute(f"SELECT {', '.join(fields)} FROM read_parquet(?)", [str(path)]).df()


def validate_site(frame):
    required = (
        "visit_id", "patient_id", "bp_phenotype", "sbp_first", "sbp_last", "arrival_elevated",
        "arrival_severe", "last_normal", "prior_htn", "has_medication_list", "prior_htn_medication",
        "ed_acute_aht", "time_first_minutes", "time_last_minutes",
    )
    if frame[list(required)].isna().any().any() or not frame["visit_id"].is_unique:
        raise ValueError("Missing required C1 inputs or non-unique C1 visit identifiers.")
    first = frame["time_first_minutes"]
    last = frame["time_last_minutes"]
    normal = frame["first_normal_minutes"]
    drug = frame["ed_first_acute_aht_minutes"]
    if (
        not np.isfinite(frame[["sbp_first", "sbp_last", "time_first_minutes", "time_last_minutes"]]).all().all()
        or (first < 0).any() or (last < first).any()
        or (normal.notna() & ((normal < first) | (normal > last))).any()
        or (drug < 0).any() or not drug.notna().equals(frame["ed_acute_aht"].astype(bool))
    ):
        raise ValueError("Invalid blood-pressure or treatment time fields.")
    expected = np.where(~frame["arrival_elevated"], "N", np.where(frame["last_normal"], "T", "P"))
    if (
        not np.array_equal(frame["bp_phenotype"], expected)
        or (frame["arrival_severe"] & ~frame["arrival_elevated"]).any()
        or (frame["bp_phenotype"].eq("T") & normal.isna()).any()
    ):
        raise ValueError("Phenotypes disagree with the arrival and final blood-pressure flags.")


def scenarios(frame, short):
    yield PRIMARY, frame
    yield "Medication list available (S11-16)", frame[frame["has_medication_list"]]
    yield "Exclude acute class A treatment (S10)", frame[~frame["ed_acute_aht"]]
    if short == "mimic":
        htn_without_ed = frame["prior_htn_medication"] | frame["prior_htn_hosp_diagnosis"]
        if (
            frame[["prior_htn_hosp_diagnosis", "prior_htn_ed_diagnosis"]].isna().any().any()
            or not (htn_without_ed | frame["prior_htn_ed_diagnosis"]).equals(frame["prior_htn"])
        ):
            raise ValueError("MIMIC hypertension source flags do not reconcile with the primary definition.")
        yield "Ignore prior ED hypertension diagnosis (S11-17)", frame.assign(prior_htn=htn_without_ed)


def wilson(count, total):
    if total == 0:
        return math.nan, math.nan, math.nan
    p = count / total
    z = cfg.WILSON_95_Z
    divisor = 1.0 + z**2 / total
    center = (p + z**2 / (2.0 * total)) / divisor
    half = z * math.sqrt(p * (1.0 - p) / total + z**2 / (4.0 * total**2)) / divisor
    return 100.0 * p, 100.0 * max(0.0, center - half), 100.0 * min(1.0, center + half)


def describe(site, analysis, frame):
    elevated = frame["arrival_elevated"]
    measures = (
        ("Arrival elevated", int(elevated.sum()), len(frame)),
        ("Arrival severe", int(frame["arrival_severe"].sum()), len(frame)),
        ("Last reading normal", int(frame["last_normal"].sum()), len(frame)),
        ("Last reading normal among arrival elevated", int((elevated & frame["last_normal"]).sum()), int(elevated.sum())),
        *((label, int((frame["bp_phenotype"] == label).sum()), len(frame)) for label in PHENOTYPES),
    )
    for label, count, total in measures:
        estimate, lo, hi = wilson(count, total)
        yield {
            "site": site, "analysis": analysis, "measure": label, "visits": len(frame),
            "patients": frame["patient_id"].nunique(), "numerator": count, "denominator": total,
            "percent": estimate, "wilson95_low_pct": lo, "wilson95_high_pct": hi,
        }


def phenotype_htn_table(site, analysis, frame):
    for phenotype in PHENOTYPES:
        group = frame[frame["bp_phenotype"] == phenotype]
        documented = int(group["prior_htn"].sum())
        estimate, lo, hi = wilson(documented, len(group))
        yield {
            "site": site, "analysis": analysis, "phenotype": phenotype,
            "visits": len(group), "patients": group["patient_id"].nunique(),
            "prior_htn_visits": documented, "no_documented_htn_visits": len(group) - documented,
            "prior_htn_pct": estimate, "wilson95_low_pct": lo, "wilson95_high_pct": hi,
        }


def fit_htn_ratio(site, analysis, frame):
    pt = frame[frame["bp_phenotype"].isin(("P", "T"))]
    persistent = pt["bp_phenotype"].eq("P")
    clusters = pt["patient_id"].nunique()
    p_n = int(persistent.sum())
    t_n = len(pt) - p_n
    p_htn = int(pt.loc[persistent, "prior_htn"].sum())
    t_htn = int(pt.loc[~persistent, "prior_htn"].sum())
    result = {
        "site": site, "analysis": analysis, "contrast": "P versus T",
        "p_visits": p_n, "t_visits": t_n, "patient_clusters": clusters,
        "p_prior_htn": p_htn, "t_prior_htn": t_htn,
        "prevalence_ratio": math.nan, "ci95_low": math.nan, "ci95_high": math.nan,
        "cluster_robust_log_se": math.nan, "status": "No positive events in one or both groups",
    }
    if not p_htn or not t_htn:
        return result
    result["prevalence_ratio"] = (p_htn / p_n) / (t_htn / t_n)
    if clusters < 2 or len(pt) <= 2:
        result["status"] = "Too few observations or patient clusters for a robust interval"
        return result
    if p_htn == p_n and t_htn == t_n:
        result["status"] = "No variation in prior hypertension"
        return result
    design = np.column_stack([np.ones(len(pt)), persistent.to_numpy(dtype=float)])
    try:
        with warnings.catch_warnings(record=True) as fit_warnings:
            warnings.simplefilter("always")
            fitted = sm.GLM(pt["prior_htn"].to_numpy(dtype=float), design, family=sm.families.Poisson()).fit(
                cov_type="cluster", cov_kwds={"groups": pt["patient_id"], "use_correction": True}
            )
    except (np.linalg.LinAlgError, PerfectSeparationError, ValueError, ZeroDivisionError):
        result["status"] = "Clustered Poisson fit unavailable"
        return result
    beta, se = fitted.params[1], fitted.bse[1]
    if not fitted.converged or not np.isfinite([beta, se]).all() or se <= 0:
        result["status"] = "Non-convergence or nonpositive robust variance"
        return result
    low, high = np.exp([beta - cfg.WILSON_95_Z * se, beta + cfg.WILSON_95_Z * se])
    result.update(
        prevalence_ratio=math.exp(beta), ci95_low=low, ci95_high=high,
        cluster_robust_log_se=se, status="OK with fit warnings" if fit_warnings else "OK",
    )
    return result


def recovery_observations(frame):
    elevated = frame[frame["arrival_elevated"]]
    before_a = elevated["ed_first_acute_aht_minutes"].isna() | (
        elevated["time_first_minutes"] < elevated["ed_first_acute_aht_minutes"]
    )
    eligible = elevated[before_a]
    last = eligible["time_last_minutes"].to_numpy(dtype=float)
    drug = eligible["ed_first_acute_aht_minutes"].fillna(np.inf).to_numpy(dtype=float)
    normal = eligible["first_normal_minutes"].fillna(np.inf).to_numpy(dtype=float)
    recovered = normal < drug
    censor_at_a = (~recovered) & (drug <= last)
    end = np.where(recovered, normal, np.minimum(last, drug))
    return elevated, eligible, recovered, censor_at_a, end


def km_at(entry, stop, recovered, targets):
    """Evaluate the product-limit estimate with delayed entry and Greenwood bounds."""
    entry = np.sort(entry)
    exit_ = np.sort(stop)
    event_times, counts = np.unique(stop[recovered], return_counts=True)

    def number_at_risk(time):
        # Same-time entry precedes recovery; all removals follow that time's risk set.
        return np.searchsorted(entry, time, side="right") - np.searchsorted(exit_, time, side="left")

    risk_at_events = number_at_risk(event_times)
    if np.any(risk_at_events < counts):
        raise ValueError("Recovery count exceeds the concurrent risk set.")
    survival = np.concatenate(([1.0], np.cumprod(1.0 - counts / risk_at_events)))
    contribution = np.divide(
        counts, risk_at_events * (risk_at_events - counts),
        out=np.zeros_like(counts, dtype=float), where=risk_at_events > counts,
    )
    greenwood = np.concatenate(([0.0], np.cumsum(contribution)))
    index = np.searchsorted(event_times, targets, side="right")
    probability = survival[index]
    lo = probability.copy()
    hi = probability.copy()
    bounded = (probability > 0) & (probability < 1)
    log_inverse = -np.log(probability[bounded])
    envelope = cfg.WILSON_95_Z * np.sqrt(greenwood[index[bounded]]) / log_inverse
    lo[bounded] = np.exp(-np.exp(np.log(log_inverse) + envelope))
    hi[bounded] = np.exp(-np.exp(np.log(log_inverse) - envelope))
    supported = (targets >= entry[0]) & (targets <= exit_[-1]) if len(entry) else np.zeros(len(targets), dtype=bool)
    at_risk = number_at_risk(targets)
    reached = np.flatnonzero(survival[1:] <= 0.5)
    median = event_times[reached[0]] if len(reached) else math.nan
    return (
        at_risk,
        np.where(supported, 100.0 * (1.0 - probability), np.nan),
        np.where(supported, 100.0 * (1.0 - hi), np.nan),
        np.where(supported, 100.0 * (1.0 - lo), np.nan),
        median,
    )


def recovery(site, analysis, documented, frame, targets):
    elevated, eligible, recovered, censor_at_a, stop = recovery_observations(frame)
    entry = eligible["time_first_minutes"].to_numpy(dtype=float)
    target_min = np.asarray(targets, dtype=float)
    at_risk, prob, lo, hi, median = km_at(entry, stop, recovered, target_min)
    keys = {"site": site, "analysis": analysis, "prior_htn": documented}
    summary = {
        **keys, "arrival_elevated_visits": len(elevated),
        "arrival_elevated_patients": elevated["patient_id"].nunique(),
        "excluded_without_pre_a_elevated_baseline": len(elevated) - len(eligible),
        "km_visits": len(eligible), "km_patients": eligible["patient_id"].nunique(),
        "recovered_visits": int(recovered.sum()), "censored_at_acute_a": int(censor_at_a.sum()),
        "censored_at_last_bp": len(eligible) - int(recovered.sum()) - int(censor_at_a.sum()),
        "median_recovery_minutes": median,
    }
    points = [
        {**keys, "minutes_from_arrival": minute, "at_risk_visits": n,
         "recovered_pct": p, "pointwise95_low_pct": lower, "pointwise95_high_pct": upper}
        for minute, n, p, lower, upper in zip(target_min, at_risk, prob, lo, hi)
    ]
    return summary, points


def paired_mimic_change(analysis, documented, frame):
    count = len(frame)
    clusters = frame["patient_id"].nunique()
    delta = (frame["sbp_last"] - frame["sbp_first"]).to_numpy(dtype=float)
    q1, median, q3 = np.quantile(delta, [0.25, 0.50, 0.75]) if count else (math.nan,) * 3
    mean = float(delta.mean()) if count else math.nan
    se = low = high = math.nan
    if clusters >= 2 and count >= 2:
        fitted = sm.OLS(delta, np.ones((count, 1))).fit(
            cov_type="cluster", cov_kwds={"groups": frame["patient_id"], "use_correction": True}
        )
        se = fitted.bse[0]
        q = student_t.ppf(0.975, clusters - 1)
        low, high = mean - q * se, mean + q * se
    return {
        "site": "MIMIC-IV-ED", "analysis": analysis, "prior_htn": documented,
        "visits": count, "patients": clusters, "mean_change_mmhg": mean,
        "sd_change_mmhg": delta.std(ddof=1) if count > 1 else math.nan,
        "cluster_robust_mean95_low_mmhg": low, "cluster_robust_mean95_high_mmhg": high,
        "median_change_mmhg": median, "q25_change_mmhg": q1, "q75_change_mmhg": q3,
    }


def analyze(site_frames):
    descriptive, cross, ratios, km_summary, km_points, km_grid, paired = [], [], [], [], [], [], []
    plot_times = np.arange(0, max(cfg.AIM1_KM_TIMEPOINTS_MINUTES) + 1, cfg.AIM1_KM_PLOT_STEP_MINUTES)
    report_times = sorted(set(cfg.AIM1_KM_TIMEPOINTS_MINUTES) | set(cfg.AIM1_KM_LANDMARKS_MINUTES))
    for short, site in SITES:
        frame = site_frames[short]
        validate_site(frame)
        for name, included in scenarios(frame, short):
            descriptive.extend(describe(site, name, included))
            cross.extend(phenotype_htn_table(site, name, included))
            ratios.append(fit_htn_ratio(site, name, included))
            for documented in HTN_GROUPS:
                group = included[included["prior_htn"] == documented]
                targets = sorted(set(report_times) | set(plot_times)) if name == PRIMARY else report_times
                summary, values = recovery(site, name, documented, group, targets)
                km_summary.append(summary)
                km_points.extend(point for point in values if point["minutes_from_arrival"] in report_times)
                if name == PRIMARY:
                    km_grid.extend(point for point in values if point["minutes_from_arrival"] in plot_times)
                if short == "mimic":
                    paired.append(paired_mimic_change(name, documented, group))
    return {
        "p2_aim1_proportions": descriptive,
        "p2_aim1_phenotype_htn": cross,
        "p2_aim1_prevalence_ratio": ratios,
        "p2_aim1_km_summary": km_summary,
        "p2_aim1_km_timepoints": km_points,
        "p2_aim1_km_grid": km_grid,
        "p2_aim1_mimic_paired_sbp": paired,
    }


def run(derived=ROOT / "derived", tables=ROOT / "outputs" / "tables"):
    with duckdb.connect() as con:
        site_frames = {short: load_site(con, derived, short) for short, _ in SITES}
    results = analyze(site_frames)
    tables.mkdir(parents=True, exist_ok=True)
    for name, rows in results.items():
        pd.DataFrame(rows).to_csv(tables / f"{name}.csv", index=False, float_format="%.6f")


if __name__ == "__main__":
    try:
        run()
    except Exception as exc:
        raise SystemExit(f"Outcome-free Aim 1 failed ({type(exc).__name__}).") from None
