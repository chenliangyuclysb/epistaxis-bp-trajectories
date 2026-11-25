"""Reconstruct P5 ED antihypertensive exposure and its aggregate timing audit."""

from importlib import import_module
from pathlib import Path

import numpy as np
import pandas as pd

import config as cfg


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "outputs" / "tables"
SITES = (("mcmed", "MC-MED"), ("mimic", "MIMIC-IV-ED"))
EXPOSURES = (
    ("Class A", "a", "ed_acute_aht", "ed_first_acute_aht_minutes", "o.aht_acute"),
    ("Class B", "b", "ed_any_aht", "ed_first_any_aht_minutes", "o.aht_any"),
    ("Alpha-1 blockers", "alpha1", "ed_alpha1_sensitivity", "ed_first_alpha1_minutes", "o.aht_alpha1_sensitivity"),
    ("Class B + alpha-1", "b_alpha1", "ed_any_or_alpha1_aht", "ed_first_any_or_alpha1_minutes",
     "(o.aht_any OR o.aht_alpha1_sensitivity)"),
)
PREVIOUS_FLAGS = ("ed_acute_aht", "ed_any_aht", "ed_alpha1_sensitivity", "ed_anxiolytic", "ed_opioid")
PREVIOUS_TIMES = ("ed_first_acute_aht_minutes", "ed_first_any_aht_minutes")
ORDER_QA_FIELDS = (
    "placed_ed", "placed_ed_no_documented_admin", "placed_ed_invalid_admin_timestamp", "placed_ed_admin_outside_ed",
)


def load_readings(con, frame):
    p2 = import_module("03_p2_phenotypes")
    return {
        site: p2.ordered_groups(p2.load_reads(con, short, frame.loc[frame.site.eq(site)]))
        for short, site in SITES
    }


def stage_dictionaries(con):
    dictionary = import_module("02_p1_dictionaries")
    dictionary.map_labels(
        con, "mcmed_orders", ("Order_type", "Procedure_name"), "p5_mcmed_medication_labels", dictionary.ED_FIELDS,
        lambda kind, name: dictionary.classify_aht_ed(name, kind), None,
    )
    dictionary.map_labels(
        con, "mimic_pyxis", ("name",), "p5_mimic_medication_labels", dictionary.ED_FIELDS,
        lambda name: dictionary.classify_aht_ed(name, mimic=True), None,
    )


def stage_medication_records(con, short):
    flags = "d.aht_acute, d.aht_any, d.aht_alpha1_sensitivity, d.anxiolytic, d.opioid"
    recognized = "d.aht_any OR d.aht_alpha1_sensitivity OR d.anxiolytic OR d.opioid"
    if short == "mcmed":
        source = f"""
            SELECT o.CSN AS visit_id, {flags}, try_cast(o.Order_time AS TIMESTAMP) AS order_ts,
                try_cast(o.First_admin_time AS TIMESTAMP) AS admin_ts,
                nullif(trim(o.First_admin_time), '') IS NULL AS admin_blank
            FROM mcmed_orders o
            JOIN p5_mcmed_medication_labels d ON o.Order_type IS NOT DISTINCT FROM d.Order_type
                AND o.Procedure_name IS NOT DISTINCT FROM d.Procedure_name
            WHERE upper(trim(o.Order_type)) = 'MEDICATIONS' AND ({recognized})
        """
    else:
        source = f"""
            SELECT o.stay_id AS visit_id, {flags}, try_cast(o.charttime AS TIMESTAMP) AS admin_ts
            FROM mimic_pyxis o JOIN p5_mimic_medication_labels d USING (name)
            WHERE {recognized}
        """
    con.execute(f"""
        CREATE OR REPLACE TEMP TABLE p5_{short}_medication_records AS
        WITH labels AS ({source})
        SELECT i.visit_id, i.arrival_ts, i.departure_ts, t.* EXCLUDE (visit_id)
        FROM labels t JOIN p5_exposure_index i USING (visit_id)
    """)


def collect_exposures(con, short):
    admin = "o.admin_ts BETWEEN i.arrival_ts AND i.departure_ts"
    metrics = [
        f"count(*) FILTER (WHERE ({admin}) AND o.anxiolytic) > 0 AS ed_anxiolytic",
        f"count(*) FILTER (WHERE ({admin}) AND o.opioid) > 0 AS ed_opioid",
    ]
    for _, code, flag, relative_time, class_filter in EXPOSURES:
        eligible_dose = f"({admin}) AND {class_filter}"
        metrics.extend((
            f"count(*) FILTER (WHERE {eligible_dose}) > 0 AS {flag}",
            f"min(epoch(o.admin_ts-i.arrival_ts)/60.0) FILTER (WHERE {eligible_dose}) AS {relative_time}",
        ))
        if short == "mcmed":
            for name, condition in zip(ORDER_QA_FIELDS, (
                "true", "o.admin_blank", "NOT o.admin_blank AND o.admin_ts IS NULL",
                f"o.admin_ts IS NOT NULL AND NOT ({admin})",
            )):
                metrics.append(f"""
                    count(*) FILTER (WHERE {class_filter} AND (o.order_ts BETWEEN i.arrival_ts AND i.departure_ts)
                        AND ({condition})) AS {code}_{name}
                """)
    return con.execute(f"""
        SELECT i.visit_id, {', '.join(metrics)}
        FROM p5_exposure_index i
        LEFT JOIN p5_{short}_medication_records o USING (visit_id)
        GROUP BY i.visit_id
        ORDER BY i.visit_id
    """).df()


def verify_prior_exposures(con, short, current):
    fields = (*PREVIOUS_FLAGS, *PREVIOUS_TIMES)
    prior = con.execute(f"SELECT visit_id, {', '.join(fields)} FROM read_parquet('derived/cohort_{short}.parquet')").df()
    paired = current.merge(prior, how="left", on="visit_id", suffixes=("", "_p1"), validate="one_to_one")
    for name in PREVIOUS_FLAGS:
        if paired[f"{name}_p1"].isna().any() or not paired[name].eq(paired[f"{name}_p1"]).all():
            raise ValueError("The reconstructed ED medication exposure disagrees with the P1 dictionary output.")
    for name in PREVIOUS_TIMES:
        if not np.allclose(paired[name], paired[f"{name}_p1"], equal_nan=True, rtol=0, atol=1e-9):
            raise ValueError("The reconstructed ED first-dose time disagrees with the P1 output.")


def verify_inclusion(frame):
    invalid = frame.ed_acute_aht & ~frame.ed_any_aht
    invalid |= frame.ed_any_or_alpha1_aht.ne(frame.ed_any_aht | frame.ed_alpha1_sensitivity)
    duration = (frame.departure_ts - frame.arrival_ts).dt.total_seconds() / 60.0
    for _, _, flag, relative_time, _ in EXPOSURES:
        time = frame[relative_time]
        invalid |= frame[flag].ne(time.notna()) | time.lt(0) | time.gt(duration)
    invalid |= frame.ed_acute_aht & frame.ed_first_any_aht_minutes.gt(frame.ed_first_acute_aht_minutes)
    union_time = frame[["ed_first_any_aht_minutes", "ed_first_alpha1_minutes"]].min(axis=1)
    if invalid.any() or not np.allclose(frame.ed_first_any_or_alpha1_minutes, union_time, equal_nan=True, rtol=0, atol=1e-9):
        raise ValueError("P5 medication flags or first-dose inclusion conditions failed.")


def bp_features(frame, readings):
    p2 = import_module("03_p2_phenotypes")
    features = []
    times = [exposure[3] for exposure in EXPOSURES]
    fields = ["site", "visit_id", "arrival_ts", "departure_ts", *times]
    for row in frame[fields].itertuples(index=False):
        seq = readings[row.site].get(row.visit_id, np.empty((0, 3), dtype=float))
        valid_minutes = min((row.departure_ts - row.arrival_ts).total_seconds() / 60, cfg.ED_MAX_HOURS * 60)
        if not np.isfinite(seq).all() or (seq[:, 0] < 0).any() or (seq[:, 0] > valid_minutes).any():
            raise ValueError("P5 BP values or cleaned observation windows did not reconcile.")
        full_peak = float(seq[:, 1].max()) if len(seq) else np.nan
        arrival = seq[0] if len(seq) else (np.nan, np.nan, np.nan)
        record = {
            "site": row.site, "visit_id": row.visit_id, "peak_sbp_full": full_peak,
            "p5_first_bp_minutes": arrival[0], "p5_first_sbp": arrival[1], "p5_first_dbp": arrival[2],
            "arrival_elevated": bool(p2.elevated(arrival[1], arrival[2])) if len(seq) else None,
        }
        for _, code, _, time, _ in EXPOSURES:
            dose = getattr(row, time)
            pre = seq[seq[:, 0] < dose] if np.isfinite(dose) else np.empty((0, 3), dtype=float)
            peak = float(pre[:, 1].max()) if len(pre) else np.nan
            # Preserve the first reading in the latest timestamp, as in P2.
            last_position = int(np.searchsorted(pre[:, 0], pre[-1, 0])) if len(pre) else 0
            last = pre[last_position] if len(pre) else (np.nan, np.nan, np.nan)
            record.update({
                f"peak_sbp_before_{code}": peak,
                f"peak_sbp_for_{code}": peak if np.isfinite(dose) else full_peak,
                f"last_sbp_before_{code}": last[1],
                f"last_dbp_before_{code}": last[2],
                f"last_bp_to_dose_minutes_{code}": dose - last[0],
                f"severe_bp_before_{code}": bool(p2.severe(pre[:, 1], pre[:, 2]).any()) if len(pre) else None,
                f"bp_tied_to_dose_{code}": bool(np.any(seq[:, 0] == dose)),
                f"last_pre_dose_timestamp_tied_{code}": len(pre) - last_position > 1,
            })
        features.append(record)
    features = pd.DataFrame(features)
    features["arrival_elevated"] = features.arrival_elevated.astype("boolean")
    replaced = set(features.columns) & (set(frame.columns) - {"site", "visit_id"})
    return frame.drop(columns=sorted(replaced)).merge(features, how="left", on=["site", "visit_id"], validate="one_to_one")


def prepare(con, frame, readings=None):
    """Return an in-memory index frame with the common P5 doses and BP summaries."""
    step = "input eligibility"
    try:
        known_sites = [site for _, site in SITES]
        if (not frame.site.isin(known_sites).all() or frame.duplicated(["site", "visit_id"]).any()
                or frame[["arrival_ts", "departure_ts"]].isna().any().any()):
            raise ValueError("P5 input requires unique complete index time bounds.")
        if not frame.departure_ts.gt(frame.arrival_ts).all():
            raise ValueError("P5 input contains invalid ED time bounds.")
        step = "local D-AHT-ED dictionary"
        stage_dictionaries(con)
        reconstructed = []
        for short, site in SITES:
            selected = frame.loc[frame.site.eq(site)]
            step = f"{site} medication times"
            con.register("p5_exposure_index", selected[["visit_id", "arrival_ts", "departure_ts"]])
            try:
                stage_medication_records(con, short)
                current = collect_exposures(con, short)
            finally:
                con.unregister("p5_exposure_index")
            verify_prior_exposures(con, short, current)
            removed = set(current.columns) & (set(selected.columns) - {"visit_id"})
            reconstructed.append(selected.drop(columns=sorted(removed)).merge(
                current, on="visit_id", how="left", validate="one_to_one",
            ))
        step = "dose and BP quality checks"
        frame = pd.concat(reconstructed, ignore_index=True)
        verify_inclusion(frame)
        frame = bp_features(frame, load_readings(con, frame) if readings is None else readings)
        c1 = frame.loc[frame.is_c1]
        if (c1.arrival_elevated.isna().any() or not np.allclose(c1.p5_first_sbp, c1.sbp_first, rtol=0, atol=1e-9)
                or not np.allclose(c1.p5_first_dbp, c1.dbp_first, rtol=0, atol=1e-9)):
            raise ValueError("The P5 C1 arrival blood pressure disagrees with P2 phenotypes.")
    except Exception:
        raise RuntimeError(f"P5 exposure preparation failed at {step}.") from None
    return frame


def quartiles(values):
    if len(values) < cfg.P5_TIMING_MIN_SUMMARY_N:
        return (None,) * len(cfg.P5_TIMING_QUANTILES)
    return tuple(values.quantile(cfg.P5_TIMING_QUANTILES))


def frequency_and_timing(frame):
    wilson = import_module("03_p2_aim1").wilson
    frequency, timing, undelivered = [], [], []
    for _, site in SITES:
        all_index = frame.loc[frame.site.eq(site)]
        for population, sample in (("Index", all_index), ("C1", all_index.loc[all_index.is_c1])):
            for label, code, flag, time, _ in EXPOSURES:
                doses = sample.loc[sample[flag]]
                percent, ci_low, ci_high = wilson(len(doses), len(sample))
                frequency.append((site, population, label, len(sample), sample.patient_id.nunique(), len(doses),
                                  doses.patient_id.nunique(), percent, ci_low, ci_high))
                sbp = doses[f"last_sbp_before_{code}"].dropna()
                gap = doses[f"last_bp_to_dose_minutes_{code}"].dropna()
                severe = doses[f"severe_bp_before_{code}"]
                if len(sbp) != len(gap) or len(sbp) != severe.notna().sum():
                    raise ValueError("P5 pre-dose BP and severe-elevation denominators differ.")
                n_severe = int(severe.eq(True).sum())
                severe_percent, severe_low, severe_high = wilson(n_severe, len(sbp))
                timing.append((site, population, label, len(doses), *quartiles(doses[time]), len(sbp),
                               len(doses) - len(sbp), int(doses[f"bp_tied_to_dose_{code}"].sum()),
                               int(doses[f"last_pre_dose_timestamp_tied_{code}"].sum()),
                               *quartiles(sbp), *quartiles(gap), n_severe, severe_percent, severe_low, severe_high))
                if site == "MC-MED":
                    placed, undelivered_col, bad_time, outside = (f"{code}_{name}" for name in ORDER_QA_FIELDS)
                    any_undelivered = sample[undelivered_col].gt(0)
                    undelivered.append((site, population, label, len(sample), int(sample[placed].sum()),
                                        int(sample[undelivered_col].sum()), int(any_undelivered.sum()),
                                        sample.loc[any_undelivered].patient_id.nunique(),
                                        int((any_undelivered & ~sample[flag]).sum()),
                                        int((any_undelivered & sample[flag]).sum()),
                                        int(sample[bad_time].sum()), int(sample[outside].sum())))
    return frequency, timing, undelivered


def write_tables(frame, output_dir=TABLES):
    """Export aggregate frequency and pre-dose timing by center and index/C1 denominator."""
    write_csv = import_module("01_p0_audit").write_csv
    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    frequency, timing, undelivered = frequency_and_timing(frame)
    keys = ("site", "population", "exposure_definition")
    write_csv(directory / "p5_exposure_frequency.csv", (*keys,
        "eligible_visits", "eligible_patients", "treated_visits", "treated_patients", "treated_percent",
        "wilson_ci95_low", "wilson_ci95_high"), frequency)
    write_csv(directory / "p5_exposure_timing.csv", (*keys,
        "treated_visits", "first_dose_minutes_q25", "first_dose_minutes_median", "first_dose_minutes_q75",
        "strict_pre_dose_bp_visits", "no_strict_pre_dose_bp_visits", "same_time_bp_visits", "last_pre_dose_timestamp_tied_visits",
        "last_pre_dose_sbp_q25", "last_pre_dose_sbp_median", "last_pre_dose_sbp_q75",
        "last_pre_dose_bp_gap_minutes_q25", "last_pre_dose_bp_gap_minutes_median", "last_pre_dose_bp_gap_minutes_q75",
        "severe_before_dose_visits", "severe_before_dose_percent_observed", "severe_wilson_ci95_low",
        "severe_wilson_ci95_high"), timing)
    write_csv(directory / "p5_exposure_mcmed_no_administration.csv", (*keys,
        "eligible_visits", "qualifying_medication_orders_placed_in_ed", "orders_without_documented_administration",
        "visits_with_orders_without_documented_administration", "patients_with_orders_without_documented_administration",
        "visits_with_no_admin_order_and_no_in_class_ed_dose", "visits_with_no_admin_order_and_in_class_ed_dose",
        "orders_with_unparseable_administration_timestamp", "orders_with_administration_outside_ed"), undelivered)
