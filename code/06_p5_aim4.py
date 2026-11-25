"""Summarize ED treatment, treatment selection and 7-day returns for P5."""

import csv
from importlib import import_module
from math import exp, log, sqrt
from pathlib import Path
import sys

import duckdb
import numpy as np
from scipy.stats import fisher_exact

import config as cfg


DEFINITIONS = (
    ("A", "ed_acute_aht", "peak_sbp_for_a"),
    ("B", "ed_any_aht", "peak_sbp_for_b"),
    ("B + alpha-1", "ed_any_or_alpha1_aht", "peak_sbp_for_b_alpha1"),
)
OUTCOME = "return_7d_event"
TERMS = {
    "peak_sbp_units": "Peak SBP, per 10 mmHg",
    "prior_htn": "Prior hypertension recorded",
    "age_years": "Age, per year",
    "is_mimic": "MIMIC-IV-ED vs MC-MED",
}


def populations(frame, sites):
    return (("Pooled", frame), *((site, frame.loc[frame.site.eq(site)]) for site in sites))


def patient_count(frame):
    return len(frame[["site", "patient_id"]].drop_duplicates())


def validate(frame):
    if frame[["site", "visit_id", "patient_id"]].isna().any().any() or frame.duplicated(["site", "visit_id"]).any():
        raise ValueError("Index visits require unique site/visit keys and nonmissing patient keys.")
    flags = [flag for _, flag, _ in DEFINITIONS]
    fields = [*flags, "ed_anxiolytic", "is_c1", "is_c2"]
    if any(not frame[name].isin((True, False)).all() for name in fields):
        raise ValueError("Treatment and cohort flags must be known and binary.")
    if (frame[flags[0]] & ~frame[flags[1]]).any() or (frame[flags[1]] & ~frame[flags[2]]).any():
        raise ValueError("Treatment definitions must be nested.")
    if not frame.loc[frame.is_c2, OUTCOME].isin((True, False)).all():
        raise ValueError("Each C2 visit must have a binary 7-day outcome.")
    c1 = frame.loc[frame.is_c1]
    expected = c1.sbp_first.ge(cfg.BP_ELEVATED_SBP) | c1.dbp_first.ge(cfg.BP_ELEVATED_DBP)
    if not c1.arrival_elevated.eq(expected).all() or not c1.arrival_elevated.isin((True, False)).all():
        raise ValueError("Arrival elevation must agree with the frozen C1 definition.")


def treatment_features(frame, models, wilson):
    features = []
    estimates = []
    missingness = []
    continuous = ("peak_sbp_units", "prior_htn", "age_years", "is_mimic")
    for definition, flag, peak in DEFINITIONS:
        for site, stratum in populations(frame, models.SITES):
            for treated, group in stratum.groupby(flag, sort=True):
                population = (definition, site, "Treated" if treated else "Untreated", len(group), patient_count(group))
                for variable, label in (("age_years", "Age (years)"), ("peak_sbp_full", "Full ED peak SBP (mmHg)"),
                                        (peak, "Model peak SBP (mmHg)")):
                    values = group[variable].dropna()
                    quartiles = (values.quantile(cfg.P5_TIMING_QUANTILES)
                                 if len(values) >= cfg.P5_TIMING_MIN_SUMMARY_N else (None, None, None))
                    features.append((*population, label, len(values), *quartiles, None, None, None, None))
                for values, label in ((group.prior_htn, "Prior hypertension recorded"),
                                      (group.ed_anxiolytic, "ED anxiolytic recorded"),
                                      (group.site.eq(models.SITES[1]), "MIMIC-IV-ED")):
                    observed = int(values.notna().sum())
                    count = int(values.eq(True).sum())
                    rate = 100 * count / observed if observed else None
                    features.append((*population, label, observed, None, None, None, count, rate, *wilson(count, observed)))

        inputs = frame.assign(peak_sbp_units=frame[peak] / cfg.P5_TREATMENT_SBP_UNIT,
                              is_mimic=frame.site.eq(models.SITES[1]))
        for name in continuous:
            unavailable = int(inputs[name].isna().sum())
            missingness.append((definition, TERMS[name], len(inputs), unavailable,
                                100 * unavailable / len(inputs) if len(inputs) else None))
        if any(100 * inputs[name].isna().mean() > cfg.P4_COVARIATE_MI_TRIGGER_PERCENT for name in continuous):
            raise ValueError("Treatment model needs the prespecified covariate imputation procedure.")
        complete = inputs.dropna(subset=list(continuous))
        fit = models.fit_model(inputs, flag, exposure=None, continuous=continuous, adjustment="none")
        if fit.n != len(complete) or fit.events != int(complete[flag].sum()):
            raise ValueError("Treatment model sample counts do not reconcile.")
        for candidate in (fit, fit.cross_check):
            if candidate is None:
                continue
            for term, label in TERMS.items():
                effect = candidate.estimate(term)
                estimates.append((definition, label, candidate.method, candidate.events, candidate.n - candidate.events,
                                  *(effect[name] for name in models.EFFECT_FIELDS), *candidate.diagnostics()))

    models.write_csv("p5_treatment_features.csv", (
        "treatment_definition", "site", "group", "visits", "patients", "characteristic", "observed",
        "q1", "median", "q3", "count", "percent", "ci95_percent_low", "ci95_percent_high",
    ), features)
    models.write_csv("p5_treatment_selection_models.csv", (
        "treatment_definition", "predictor", "method", "treated_visits", "untreated_visits",
        *models.EFFECT_FIELDS, *models.MODEL_FIELDS,
    ), estimates)
    models.write_csv("p5_treatment_selection_missingness.csv", (
        "treatment_definition", "predictor", "eligible_visits", "missing_visits", "missing_percent",
    ), missingness)


def crude_risk_ratio(treated_events, treated_n, untreated_events, untreated_n):
    if not treated_n or not untreated_n:
        return None, None, None, None, "An exposure group is absent"
    treated_risk = treated_events / treated_n
    untreated_risk = untreated_events / untreated_n
    probability = float(fisher_exact(((treated_events, treated_n - treated_events),
                                     (untreated_events, untreated_n - untreated_events)), alternative="two-sided").pvalue)
    if not untreated_events:
        return None, None, None, probability, "Reference risk zero; RR and Katz interval unavailable"
    ratio = treated_risk / untreated_risk
    if not treated_events:
        return ratio, None, None, probability, "Treated risk zero; Katz interval unavailable"
    variance = 1 / treated_events - 1 / treated_n + 1 / untreated_events - 1 / untreated_n
    half_width = cfg.WILSON_95_Z * sqrt(variance)
    return ratio, exp(log(ratio) - half_width), exp(log(ratio) + half_width), probability, "Estimated"


def treatment_returns(frame, models):
    rows = []
    for definition, flag, _ in DEFINITIONS:
        for site, stratum in populations(frame, models.SITES):
            treated = stratum.loc[stratum[flag]]
            untreated = stratum.loc[~stratum[flag]]
            treated_events = int(treated[OUTCOME].sum())
            untreated_events = int(untreated[OUTCOME].sum())
            rows.append((definition, site, len(stratum), patient_count(stratum), len(treated), patient_count(treated),
                         treated_events, 100 * treated_events / len(treated) if len(treated) else None, len(untreated),
                         patient_count(untreated), untreated_events,
                         100 * untreated_events / len(untreated) if len(untreated) else None,
                         treated_events < cfg.P5_RARE_TREATED_EVENT_COUNT,
                         *crude_risk_ratio(treated_events, len(treated), untreated_events, len(untreated))))
    models.write_csv("p5_treatment_return_7d.csv", (
        "treatment_definition", "site", "visits", "patients", "treated_visits", "treated_patients", "treated_events",
        "treated_risk_percent", "untreated_visits", "untreated_patients", "untreated_events", "untreated_risk_percent",
        "treated_events_below_sap_threshold", "crude_rr", "katz_ci95_low", "katz_ci95_high", "fisher_two_sided_p",
        "status",
    ), rows)


def source_metadata(p0, label):
    original = p0.fingerprint("p5_" + label)
    ed_files = sorted(Path("data_raw/mimic-iv-ed").rglob("*"))
    ed = []
    for path in ed_files:
        if path.is_file():
            metadata = path.stat()
            ed.append((path.as_posix(), metadata.st_size, metadata.st_mtime_ns))
    return original, ed


def verify_inputs(p0, before):
    after = source_metadata(p0, "after")
    if after != before:
        raise RuntimeError("P5 source size or modification time changed.")
    with Path("outputs/raw_fingerprint_p4_after.csv").open(encoding="utf-8-sig", newline="") as stream:
        previous = list(csv.reader(stream))[1:]
    if previous != [[str(value) for value in row] for row in after[0]]:
        raise RuntimeError("Original source metadata differ from P4 closeout.")
    p0.write_csv(Path("outputs/tables/p5_source_audit.csv"), ("source_scope", "files", "size_and_mtime_unchanged"),
                 (("Original sources, P4 closeout to P5", len(after[0]), True),
                  ("Extracted ED files, during P5", len(after[1]), True)))


def main():
    step = "environment and input metadata"
    try:
        import_module("02_p1_setup").verify_project_runtime()
        p0 = import_module("01_p0_audit")
        models = import_module("05_p4_models")
        before = source_metadata(p0, "before")
        stage_4 = import_module("05_p4_aim3")
        exposure = import_module("06_p5_exposure")
        wilson = import_module("02_p1_tables").wilson
        with duckdb.connect(str(stage_4.DERIVED / "epistaxis.duckdb"), read_only=True) as con:
            con.execute(f"SET temp_directory = {p0.literal(stage_4.DERIVED / 'duckdb_tmp')}")
            step = "frozen cohorts and ED treatment definitions"
            index = stage_4.prepare_index(con)
            readings = exposure.load_readings(con, index)
            frame = exposure.prepare(con, index, readings)
            validate(frame)
            exposure.write_tables(frame)
            step = "C1 treatment selection"
            treatment_features(frame.loc[frame.is_c1 & frame.arrival_elevated.fillna(False)], models, wilson)
            step = "C1 risk-set matching for post-dose blood pressure"
            import_module("06_p5_matching").run(frame, readings)
            step = "C2 7-day crude descriptions"
            treatment_returns(frame.loc[frame.is_c2], models)
        step = "source metadata verification"
        verify_inputs(p0, before)
    except Exception as error:
        # Scanner and model errors may carry row content.
        print(f"P5 stopped at {step} ({type(error).__name__}).", file=sys.stderr)
        return 1
    print("P5 treatment and 7-day aggregate tables complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
