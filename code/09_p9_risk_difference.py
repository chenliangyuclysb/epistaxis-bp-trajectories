"""Estimate absolute differences in seven-day return risk from the frozen P4 model."""

from collections import Counter
import csv
from importlib import import_module
from pathlib import Path
import sys

import duckdb
import numpy as np
import pandas as pd
from statsmodels.stats.proportion import confint_proportions_2indep

import config as cfg


DERIVED = Path("derived")
TABLES = Path("outputs/tables")
OUTCOME = "return_7d_event"
EXPOSURE = "bp_phenotype"
LEVELS = ("T", "P", "N")
ADJUSTMENT = "pooled"
CLUSTER_KEYS = ["site", "patient_id"]


def cluster_counts(frame):
    return (
        len(frame),
        len(frame[CLUSTER_KEYS].drop_duplicates()),
        int(frame[OUTCOME].sum()),
    )


def frozen_p4_results():
    with (TABLES / "p4_main_group_counts.csv").open(encoding="utf-8-sig", newline="") as stream:
        groups = {
            row["phenotype"]: tuple(int(row[key]) for key in ("model_visits", "patients", "model_events"))
            for row in csv.DictReader(stream)
            if row["endpoint"] == "7-day return" and row["site"] == "Pooled"
        }
    with (TABLES / "p4_main_estimates.csv").open(encoding="utf-8-sig", newline="") as stream:
        estimates = {
            row["contrast"]: row
            for row in csv.DictReader(stream)
            if row["endpoint"] == "7-day return" and row["site"] == "Pooled" and row["analysis"] == "Main GEE"
        }
    return groups, estimates


def validate_point_model(fit, pooled, counts, risks):
    frozen_counts, frozen_estimates = frozen_p4_results()
    if (counts != frozen_counts or (fit.n, fit.patients, fit.events) != pooled
            or not fit.converged or risks is None):
        raise ValueError("The frozen P4 return cohort or model could not be reproduced.")
    for index, level in enumerate(LEVELS[1:], start=1):
        source = frozen_estimates[f"{level} vs T"]
        estimate = fit.estimate(f"{EXPOSURE}={level}")["estimate"]
        if (tuple(int(source[key]) for key in ("visits", "patients", "events")) != pooled
                or not np.isclose(estimate, float(source["estimate"]), rtol=1e-10, atol=0)
                or not np.isclose(risks[index] / risks[0], estimate, rtol=1e-10, atol=0)):
            raise ValueError("The risk ratios do not reproduce the reviewed pooled P4 estimates.")


def standardized_risks(models, sample, fit):
    if fit.params is None:
        return None
    design, columns, _ = models.model_design(sample, EXPOSURE, LEVELS, (), ADJUSTMENT, False)
    if columns != fit.columns:
        raise ValueError("The fitted and prediction design columns do not agree.")
    positions = [columns.index(f"{EXPOSURE}={level}") for level in LEVELS[1:]]
    design[:, positions] = 0
    with np.errstate(over="ignore", invalid="ignore"):
        # All visits receive the same phenotype indicator in each prediction.
        reference = np.exp(design @ fit.params).mean()
        risks = np.r_[reference, reference * np.exp(fit.params[positions])]
    return risks if np.isfinite(risks).all() and ((0 <= risks) & (risks <= 1)).all() else None


def cluster_bootstrap(models, sample):
    codes = pd.MultiIndex.from_frame(sample[CLUSTER_KEYS]).factorize()[0]
    counts = np.bincount(codes)
    blocks = np.split(np.argsort(codes, kind="stable"), counts.cumsum()[:-1])
    clusters = len(blocks)
    rng = np.random.default_rng(cfg.RANDOM_SEED)
    estimates = []
    failures = Counter()
    for _ in range(cfg.P9_RD_BOOTSTRAPS):
        draw = rng.integers(clusters, size=clusters)
        sampled = sample.take(np.concatenate([blocks[code] for code in draw])).copy()
        # Duplicate selections of a patient are independent bootstrap clusters.
        sampled["patient_id"] = np.repeat(np.arange(clusters), counts[draw])
        sampled["visit_id"] = np.arange(len(sampled))
        fit = models.fit_model(sampled, OUTCOME, exposure=EXPOSURE, levels=LEVELS, adjustment=ADJUSTMENT)
        risks = standardized_risks(models, sampled, fit)
        if risks is not None:
            estimates.append(risks)
        else:
            failures[fit.status if fit.params is None else "Invalid standardized risk"] += 1
    if len(estimates) < 2:
        raise ValueError("Insufficient successful patient-cluster replicates for percentile intervals.")
    return np.stack(estimates), failures


def result_rows(models, sample):
    pooled = cluster_counts(sample)
    counts = {level: cluster_counts(sample.loc[sample[EXPOSURE].eq(level)]) for level in LEVELS}
    fit = models.fit_model(sample, OUTCOME, exposure=EXPOSURE, levels=LEVELS, adjustment=ADJUSTMENT)
    risks = standardized_risks(models, sample, fit)
    validate_point_model(fit, pooled, counts, risks)
    bootstrap, failures = cluster_bootstrap(models, sample)
    alpha = 1 - cfg.P9_RD_CONFIDENCE_LEVEL
    rows = []
    for index, level in enumerate(LEVELS[1:], start=1):
        visits, patients, events = counts[level]
        ref_visits, ref_patients, ref_events = counts[LEVELS[0]]
        standard = {
            "endpoint": "7-day return",
            "site": "Pooled",
            "contrast": f"{level} vs T",
            "group": level,
            "reference": LEVELS[0],
            "group_visits": visits,
            "group_patients": patients,
            "group_events": events,
            "reference_visits": ref_visits,
            "reference_patients": ref_patients,
            "reference_events": ref_events,
            "pooled_visits": pooled[0],
            "pooled_patients": pooled[1],
            "pooled_events": pooled[2],
        }
        low, high = np.quantile(bootstrap[:, index] - bootstrap[:, 0], [alpha / 2, 1 - alpha / 2])
        rows.append(standard | {
            "analysis": "Adjusted",
            "standardization_visits": pooled[0],
            "group_risk_percent": 100 * risks[index],
            "reference_risk_percent": 100 * risks[0],
            "risk_difference_pp": 100 * (risks[index] - risks[0]),
            "ci95_low_pp": 100 * low,
            "ci95_high_pp": 100 * high,
            "fit_method": fit.method + "; marginal standardization",
            "ci_method": "Pooled patient-cluster percentile bootstrap; standardize within each resample",
            "bootstrap_requested": cfg.P9_RD_BOOTSTRAPS,
            "bootstrap_successful": len(bootstrap),
            "bootstrap_failed": sum(failures.values()),
            "bootstrap_failure_reasons": "; ".join(f"{reason} ({count})" for reason, count in sorted(failures.items())),
        })
        observed, ref_observed = events / visits, ref_events / ref_visits
        low, high = confint_proportions_2indep(events, visits, ref_events, ref_visits,
                                              method="newcomb", compare="diff", alpha=alpha)
        rows.append(standard | {
            "analysis": "Crude",
            "standardization_visits": None,
            "group_risk_percent": 100 * observed,
            "reference_risk_percent": 100 * ref_observed,
            "risk_difference_pp": 100 * (observed - ref_observed),
            "ci95_low_pp": 100 * low,
            "ci95_high_pp": 100 * high,
            "fit_method": "Observed visit risks",
            "ci_method": "Newcombe hybrid score (unclustered)",
            "bootstrap_requested": None,
            "bootstrap_successful": None,
            "bootstrap_failed": None,
            "bootstrap_failure_reasons": "",
        })
    return rows


def main():
    try:
        import_module("02_p1_setup").verify_project_runtime()
        models = import_module("05_p4_models")
        prepare_index = import_module("05_p4_aim3").prepare_index
        with duckdb.connect(str(DERIVED / "epistaxis.duckdb"), read_only=True) as con:
            index = prepare_index(con)
        frame = index.loc[index.is_c2]
        sample = models.normalized_frame(frame, OUTCOME, EXPOSURE, LEVELS, (), ADJUSTMENT)
        if len(frame) != len(sample) or sample.isna().any().any():
            raise ValueError("Every frozen return visit must qualify for the pooled model.")
        rows = result_rows(models, sample)
        fields = tuple(rows[0])
        models.write_csv("p9_risk_difference.csv", fields, [tuple(row[field] for field in fields) for row in rows])
    except Exception as error:
        print(f"P9 risk-difference analysis stopped ({type(error).__name__}).", file=sys.stderr)
        return 1
    print("P9 pooled risk-difference estimates exported.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
