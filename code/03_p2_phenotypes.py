"""Build outcome-free P2 phenotypes, dense metrics and sparse-sampling summaries."""

from importlib import import_module
from pathlib import Path
import sys

import duckdb
import numpy as np
import pandas as pd

import config as cfg


DERIVED = Path("derived")
TABLES = Path("outputs/tables")
SITES = (("mcmed", "MC-MED"), ("mimic", "MIMIC-IV-ED"))
PHENOTYPES = ("N", "T", "P")
SOURCE_ORDER = {"Triage": 0, "Monitor": 1, "Chart": 2, "Vitalsign": 2}
BASE_COLUMNS = (
    "visit_id", "patient_id", "age_years", "sex", "has_medication_list",
    "prior_htn", "prior_htn_diagnosis", "prior_htn_medication",
    "prior_htn_hosp_diagnosis", "prior_htn_ed_diagnosis",
    "ed_acute_aht", "ed_first_acute_aht_minutes", "is_c1",
)
DENSE_CONTINUOUS = ("twa_sbp", "hypertension_load_sbp", "sbp_arv", "peak_sbp")
MODEL_INPUT = (
    "visit_id", "patient_id", "minutes_from_arrival", "sbp", "prior_htn", "age_years", "sex",
    "has_medication_list", "ed_acute_aht", "ed_first_acute_aht_minutes",
)
HEAD = (
    "site", "phenotype", "metric", "n", "denominator", "pct", "wilson95_lo_pct",
    "wilson95_hi_pct", "q1", "median", "q3",
)


def elevated(sbp, dbp):
    return (sbp >= cfg.BP_ELEVATED_SBP) | (dbp >= cfg.BP_ELEVATED_DBP)


def severe(sbp, dbp):
    return (sbp >= cfg.BP_SEVERE_SBP) | (dbp >= cfg.BP_SEVERE_DBP)


def load_base(con, site):
    # The analysis projects only allowlisted columns from the audited cohort.
    return con.execute(f"""
        SELECT {", ".join(BASE_COLUMNS)},
            epoch(departure_ts - arrival_ts) / 60.0 AS ed_duration_minutes
        FROM read_parquet('derived/cohort_{site}.parquet')
        ORDER BY visit_id
    """).df()


def load_reads(con, site, base):
    con.register("p2_analysis_visits", base[["visit_id"]])
    try:
        return con.execute(f"""
            SELECT b.visit_id, b.source, b.minutes_from_arrival, b.sbp, b.dbp
            FROM read_parquet('derived/bp_long_{site}.parquet') b
            SEMI JOIN p2_analysis_visits c USING (visit_id)
            ORDER BY visit_id, minutes_from_arrival, source, sbp, dbp
        """).df()
    finally:
        con.unregister("p2_analysis_visits")


def ordered_groups(readings):
    groups = readings.assign(priority=readings.source.map(SOURCE_ORDER))
    if groups.priority.isna().any():
        raise ValueError("Unexpected BP source.")
    result = {}
    for visit_id, group in groups.groupby("visit_id", sort=True):
        group = group.sort_values(["minutes_from_arrival", "priority", "source", "sbp", "dbp"])
        result[visit_id] = group[["minutes_from_arrival", "sbp", "dbp"]].to_numpy(dtype=float)
    return result


def positive_area(times, values, threshold):
    """Integrate the positive part of the interpolated line above a threshold."""
    left = values[:-1] - threshold
    right = values[1:] - threshold
    span = np.diff(times)
    both = (left >= 0) & (right >= 0)
    crossing = ((left > 0) & (right < 0)) | ((left < 0) & (right > 0))
    positive = np.maximum(left[crossing], right[crossing])
    return float(
        np.sum(span[both] * (left[both] + right[both]) / 2)
        + np.sum(span[crossing] * positive**2 / (2 * np.abs(left[crossing] - right[crossing])))
    )


def dense_metrics(sequence):
    time, sbp, dbp = sequence.T
    length = len(time)
    span = time[-1] - time[0] if length else 0
    return {
        "n_bp": length,
        "twa_sbp": float(np.trapezoid(sbp, time) / span) if span > 0 else np.nan,
        "hypertension_load_sbp": positive_area(time, sbp, cfg.BP_ELEVATED_SBP) / span if span > 0 else np.nan,
        "sbp_arv": float(np.abs(np.diff(sbp)).mean()) if length >= cfg.BP_MIN_ARV_READINGS else np.nan,
        "peak_sbp": float(sbp.max()) if length else np.nan,
        "ever_severe": bool(severe(sbp, dbp).any()) if length else None,
    }


def first_normal(sequence):
    if not len(sequence) or not elevated(*sequence[0, 1:]):
        return np.nan
    normal_positions = np.flatnonzero(~elevated(sequence[:, 1], sequence[:, 2]))
    return float(sequence[normal_positions[0], 0]) if len(normal_positions) else np.nan


def summarize_visit(visit_id, sequence, include_dense):
    times, sbp, dbp = sequence.T
    first = sequence[0]
    last = sequence[np.searchsorted(times, times[-1])]
    arrival_elevated = bool(elevated(*first[1:]))
    last_normal = not bool(elevated(*last[1:]))
    label = ("T" if last_normal else "P") if arrival_elevated else "N"
    anytime_normal = first_normal(sequence)
    row = {
        "visit_id": visit_id,
        "bp_phenotype": label,
        "time_first_minutes": first[0],
        "sbp_first": first[1],
        "dbp_first": first[2],
        "time_last_minutes": last[0],
        "sbp_last": last[1],
        "dbp_last": last[2],
        "arrival_elevated": arrival_elevated,
        "arrival_severe": bool(severe(*first[1:])),
        "last_normal": last_normal,
        "first_normal_minutes": anytime_normal,
    }
    if include_dense:
        row.update(dense_metrics(sequence))
        row["persistent_subtype"] = ("P-rebound" if np.isfinite(anytime_normal) else "P-never") if label == "P" else None
    return row


def build_phenotypes(base, grouped, site):
    c1 = base.loc[base.is_c1].copy()
    if set(c1.visit_id) != set(grouped):
        raise ValueError("C1 visits and BP groups differ.")
    metrics = pd.DataFrame(
        [summarize_visit(visit_id, values, site == "mcmed") for visit_id, values in grouped.items()]
    )
    result = c1.merge(metrics, on="visit_id", validate="one_to_one").sort_values("visit_id")
    result["minutes_last_to_departure"] = result.ed_duration_minutes - result.time_last_minutes
    if result.minutes_last_to_departure.lt(-1e-7).any():
        raise ValueError("An eligible last BP follows departure.")
    return result


def copy_frame(con, frame, name, literal):
    path = DERIVED / name
    pending = path.with_suffix(".parquet.tmp")
    con.register("p2_parquet_export", frame)
    try:
        con.execute(f"COPY p2_parquet_export TO {literal(pending)} (FORMAT PARQUET, COMPRESSION ZSTD)")
        pending.replace(path)
    finally:
        con.unregister("p2_parquet_export")
        pending.unlink(missing_ok=True)


def alternate_mc_reads(con):
    """Reapply S4 on the fixed MC-MED C1 visits with the wider B/C window."""
    con.execute("""
        CREATE TEMP TABLE p1_mcmed_index AS
        SELECT c.visit_id AS CSN, c.arrival_ts, c.departure_ts, c.disposition AS ED_dispo,
            v.Triage_SBP, v.Triage_DBP
        FROM read_parquet('derived/cohort_mcmed.parquet') c
        JOIN mcmed_visits v ON v.CSN = c.visit_id
        WHERE c.is_c1
    """)
    con.execute("""
        CREATE TEMP TABLE p1_mimic_index AS
        SELECT NULL::VARCHAR AS stay_id, NULL::TIMESTAMP AS arrival_ts,
            NULL::TIMESTAMP AS departure_ts, NULL::VARCHAR AS disposition
        WHERE false
    """)
    core = import_module("02_p1_cohort_core")
    core.prepare_bp_slots(con)
    core.deduplicate_bp(con, near_minutes=cfg.MC_DENSE_ALT_DUP_MINUTES)
    return con.execute("""
        SELECT visit_id, source, minutes_from_arrival, sbp, dbp
        FROM p1_bp_mcmed
        ORDER BY visit_id, minutes_from_arrival, source, sbp, dbp
    """).df()


def alternate_dense_metrics(phenotypes, main_groups, alternate_groups):
    visits = []
    for row in phenotypes.itertuples(index=False):
        for window, groups in (
            (cfg.G0_BP_DUP_MINUTES, main_groups),
            (cfg.MC_DENSE_ALT_DUP_MINUTES, alternate_groups),
        ):
            values = groups.get(row.visit_id, np.empty((0, 3)))
            observations = dense_metrics(values)
            times = values[:, 0]
            observed_recovery = first_normal(values)
            recovery = (
                observed_recovery if len(times) and np.isfinite(observed_recovery)
                and (not row.ed_acute_aht or observed_recovery < row.ed_first_acute_aht_minutes)
                else np.nan
            )
            visits.append({
                "visit_id": row.visit_id,
                "window_minutes": window,
                "time_first_minutes": times[0] if len(times) else np.nan,
                "time_last_minutes": times[-1] if len(times) else np.nan,
                "recorded_first_normal_minutes": observed_recovery,
                "pretreatment_first_normal_minutes": recovery,
                "would_meet_c1": len(times) >= cfg.C1_MIN_BP_READINGS
                and times[-1] - times[0] >= cfg.C1_MIN_SPAN_MINUTES,
                **observations,
            })
    return pd.DataFrame(visits)


def coefficient_kappa(matrix):
    n = matrix.sum()
    observed = np.trace(matrix) / n
    expected = (matrix.sum(axis=0) @ matrix.sum(axis=1)) / n**2
    return float((observed - expected) / (1 - expected)) if expected < 1 else float("nan")


def monte_carlo(mc_groups, mimic_groups):
    donor = np.asarray(
        [sum(sequence[:, 0] > 0) for sequence in mimic_groups.values()],
        dtype=int,
    )
    generator = np.random.default_rng(cfg.RANDOM_SEED)
    plans = []
    result = []
    confusion = np.zeros((3, 3), dtype=np.int64)
    for sequence in mc_groups.values():
        if sequence[0, 0] != 0:
            raise ValueError("MC-MED C1 is missing its arrival anchor for sparse sampling.")
        times, sbp, dbp = sequence.T
        is_high = elevated(sbp, dbp)
        post = np.flatnonzero(times > 0)
        available = donor[(donor >= cfg.SPARSE_MIN_POST_TRIAGE_READINGS) & (donor <= len(post))]
        if not len(available):
            raise ValueError("No MIMIC post-triage counts fit one MC-MED visit.")
        last = np.searchsorted(times, times[-1])
        real = (2 if is_high[last] else 1) if is_high[0] else 0
        full_recovery = bool(is_high[0] and (~is_high[post]).any())
        plans.append((times, is_high, post, available, real, full_recovery))
    full_recoveries = sum(plan[-1] for plan in plans)
    for replication in range(cfg.SPARSE_BP_REPETITIONS):
        matrix = np.zeros((3, 3), dtype=np.int64)
        detected = 0
        for times, is_high, post, available, real, full_recovery in plans:
            k = generator.choice(available)
            chosen = generator.choice(post, size=k, replace=False)
            last = chosen[times[chosen] == times[chosen].max()].min()
            sparse = (2 if is_high[last] else 1) if is_high[0] else 0
            matrix[real, sparse] += 1
            detected += full_recovery and bool((~is_high[chosen]).any())
        confusion += matrix
        result.append((replication + 1, coefficient_kappa(matrix), int(detected), full_recoveries))
    return result, confusion, donor


def number_stats(values):
    series = pd.Series(values).dropna().astype(float)
    if not len(series):
        return 0, None, None, None
    q1, median, q3 = series.quantile([0.25, 0.5, 0.75])
    return len(series), float(q1), float(median), float(q3)


def write_tables(p0, full, base, variants, sparse):
    p1 = import_module("02_p1_tables")
    def proportion(site, group, metric, n, denominator):
        return (site, group, metric, n, denominator, p1.percent(n, denominator),
                *p1.wilson(n, denominator), None, None, None)

    table2, htn, alignment = [], [], []
    for short, site in SITES:
        frame = full[short]
        for metric in ("arrival_elevated", "arrival_severe", "last_normal"):
            table2.append(proportion(site, "All C1", metric, int(frame[metric].sum()), len(frame)))
        for label in PHENOTYPES:
            table2.append(proportion(site, "All C1", label, int(frame.bp_phenotype.eq(label).sum()), len(frame)))
        if short == "mcmed":
            persistent = frame.loc[frame.bp_phenotype.eq("P")]
            for subtype in ("P-never", "P-rebound"):
                n = int(persistent.persistent_subtype.eq(subtype).sum())
                table2.append(proportion(site, "Persistent", subtype, n, len(persistent)))
                table2.append(proportion(site, "All C1", subtype, n, len(frame)))
            for label in ("All C1", *PHENOTYPES):
                group = frame if label == "All C1" else frame.loc[frame.bp_phenotype.eq(label)]
                for metric in (*DENSE_CONTINUOUS, "n_bp"):
                    n, q1, median, q3 = number_stats(group[metric])
                    table2.append((site, label, metric, n, len(group), None, None, None, q1, median, q3))
                table2.append(proportion(site, label, "ever_severe", int(group.ever_severe.sum()), len(group)))
        minutes = frame.minutes_last_to_departure
        values = minutes.quantile([0.05, 0.25, 0.5, 0.75, 0.95]).tolist()
        alignment.append((site, len(frame), *values))

        sources = (
            {"PMH diagnosis": "prior_htn_diagnosis", "Home medication": "prior_htn_medication"}
            if short == "mcmed" else {
                "Medrecon medication": "prior_htn_medication",
                "Prior hospital diagnosis": "prior_htn_hosp_diagnosis",
                "Prior ED diagnosis": "prior_htn_ed_diagnosis",
            }
        )
        for population, cohort in (("Index", base[short]), ("C1", frame)):
            counts = cohort[list(sources.values())].fillna(False).astype(bool)
            if not np.array_equal(counts.any(axis=1).to_numpy(), cohort.prior_htn.to_numpy()):
                raise ValueError("HTN union and source flags disagree.")
            diagnosed = int(cohort.prior_htn.sum())
            for label, field in sources.items():
                n = int(counts[field].sum())
                only = int((counts[field] & counts.sum(axis=1).eq(1)).sum())
                htn.append((site, population, label, len(cohort), diagnosed, n, p1.percent(n, len(cohort)),
                            *p1.wilson(n, len(cohort)), only, p1.percent(only, len(cohort)),
                            p1.percent(only, diagnosed)))
            htn.append((site, population, "Any documented source", len(cohort), diagnosed, diagnosed,
                        p1.percent(diagnosed, len(cohort)), *p1.wilson(diagnosed, len(cohort)), None, None, None))
    p0.write_csv(TABLES / "p2_table2.csv", HEAD, table2)
    p0.write_csv(TABLES / "p2_htn_sources.csv",
                 ("site", "cohort", "source", "visits", "prior_htn_visits", "source_visits", "source_pct",
                  "source_wilson95_lo_pct", "source_wilson95_hi_pct", "source_only_visits",
                  "source_only_pct_of_all", "source_only_pct_of_htn"), htn)
    p0.write_csv(TABLES / "p2_last_bp_to_departure.csv", ("site", "visits", "p05_minutes", "q1_minutes",
                                                         "median_minutes", "q3_minutes", "p95_minutes"), alignment)

    dense_rows = []
    for window, group in variants.groupby("window_minutes"):
        for metric in ("n_bp", *DENSE_CONTINUOUS, "recorded_first_normal_minutes", "pretreatment_first_normal_minutes"):
            n, q1, median, q3 = number_stats(group[metric])
            dense_rows.append((window, len(group), metric, n, q1, median, q3, None))
        for metric in ("would_meet_c1", "ever_severe"):
            dense_rows.append((window, len(group), metric, int(group[metric].count()), None, None, None,
                               int(group[metric].sum())))
    p0.write_csv(TABLES / "p2_dense_dedup_sensitivity.csv", ("window_minutes", "fixed_c1_visits", "metric", "n",
                                                            "q1", "median", "q3", "true_count"), dense_rows)

    repetitions, confusion, donor = sparse
    kappa = np.array([row[1] for row in repetitions])
    sensitivities = np.array([row[2] / row[3] if row[3] else np.nan for row in repetitions])
    sparse_summary = []
    for label, values in (("Cohen kappa", kappa), ("Sensitivity of ever-normal detection", sensitivities)):
        low, high = np.quantile(values, cfg.SPARSE_BP_INTERVAL_QUANTILES)
        sparse_summary.append((label, len(values), cfg.RANDOM_SEED, float(np.mean(values)), float(low), float(high)))
    p0.write_csv(TABLES / "p2_sparse_summary.csv", ("metric", "replications", "seed", "mean", "p025", "p975"),
                 sparse_summary)
    p0.write_csv(TABLES / "p2_sparse_repetitions.csv", ("replication", "kappa", "detected_ever_normal",
                                                        "full_data_ever_normal"), repetitions)
    p0.write_csv(TABLES / "p2_sparse_confusion.csv", ("full_data", "sparse", "pooled_count", "mean_per_replication",
                                                     "row_pct"), [
        (true, predicted, confusion[i, j], confusion[i, j] / len(repetitions),
         p1.percent(confusion[i, j], confusion[i, :].sum()))
        for i, true in enumerate(PHENOTYPES) for j, predicted in enumerate(PHENOTYPES)
    ])
    keys, counts = np.unique(donor, return_counts=True)
    p0.write_csv(TABLES / "p2_sparse_mimic_post_triage_counts.csv", ("post_triage_readings", "stays"),
                 zip(keys.tolist(), counts.tolist()))


def write_dense_km_sensitivity(p0, frame, variants):
    recovery = import_module("03_p2_aim1").recovery
    common = frame[["visit_id", "patient_id", "arrival_elevated", "prior_htn", "ed_first_acute_aht_minutes"]]
    summaries, points = [], []
    for window, group in variants.groupby("window_minutes"):
        selection = group.rename(columns={"recorded_first_normal_minutes": "first_normal_minutes"})
        selection = selection.merge(common, on="visit_id", validate="one_to_one")
        for documented in (False, True):
            subset = selection.loc[selection.prior_htn.eq(documented)]
            summary, estimates = recovery(
                "MC-MED", f"±{window}-minute deduplication", documented, subset,
                cfg.AIM1_KM_LANDMARKS_MINUTES,
            )
            summaries.append(summary)
            points.extend(estimates)
    for filename, rows in (("p2_dense_km_sensitivity.csv", summaries),
                           ("p2_dense_km_sensitivity_points.csv", points)):
        p0.write_csv(TABLES / filename, tuple(rows[0]), [tuple(row.values()) for row in rows])


def run():
    p0 = import_module("01_p0_audit")
    full, base, groups, readings = {}, {}, {}, {}
    with duckdb.connect("derived/epistaxis.duckdb", read_only=True) as con:
        con.execute(f"SET temp_directory = {p0.literal(DERIVED / 'duckdb_tmp')}")
        for short, _ in SITES:
            base[short] = load_base(con, short)
            readings[short] = load_reads(con, short, base[short].loc[base[short].is_c1])
            groups[short] = ordered_groups(readings[short])
            full[short] = build_phenotypes(base[short], groups[short], short)
        alternate = alternate_mc_reads(con)
        variants = alternate_dense_metrics(full["mcmed"], groups["mcmed"], ordered_groups(alternate))
        sparse = monte_carlo(groups["mcmed"], groups["mimic"])
        for short, _ in SITES:
            copy_frame(con, full[short], f"phenotypes_{short}.parquet", p0.literal)
        copy_frame(con, variants, "p2_mcmed_dense_dedup_sensitivity.parquet", p0.literal)
        model = readings["mcmed"].merge(full["mcmed"], on="visit_id", validate="many_to_one")
        copy_frame(con, model[list(MODEL_INPUT)], "p2_mcmed_model_input.parquet", p0.literal)
        write_tables(p0, full, base, variants, sparse)
        write_dense_km_sensitivity(p0, full["mcmed"], variants)


if __name__ == "__main__":
    try:
        import_module("02_p1_setup").verify_project_runtime()
        if not (TABLES / "p2_g1_updates.csv").is_file():
            raise RuntimeError("Apply the G1 update before Aim 1.")
        run()
        print("P2 outcome-free phenotypes, dense metrics and Monte Carlo summaries complete.")
    except Exception as error:
        # Avoid emitting parser details that may contain source data values.
        print(f"P2 phenotype build stopped ({type(error).__name__}).", file=sys.stderr)
        raise SystemExit(1) from None
