"""Run P4 landmark, inpatient and exploratory secondary analyses."""

from importlib import import_module
from pathlib import Path

import numpy as np
import pandas as pd

import config as cfg


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "outputs" / "tables"
MCMED = "MC-MED"
MIMIC = "MIMIC-IV-ED"
LANDMARK_PHENOTYPE = "landmark_phenotype"
ARRIVAL_SBP = "arrival_sbp_per_10"
LANDMARK_EVENT = "landmark_event"
INPATIENT_EVENT = "inpatient_ent_hemostasis_event"
RETURN_EVENT = "return_7d_event"
PHENOTYPE_LEVELS = ("T", "P", "N")
LANDMARK_LEVELS = ("T60", "P60", "N60")
BURDEN_LEVELS = ("Tertile 1", "Tertile 2", "Tertile 3")
HEMOSTASIS_LABEL = "Landmark hemostasis escalation"
INPATIENT_LABEL = "Inpatient ENT hemostasis procedure"
RETURN_LABEL = "7-day epistaxis return"


def events(frame, outcome):
    return int(frame[outcome].eq(True).sum())


def require_binary(frame, outcome):
    if not frame[outcome].isin((False, True)).all():
        raise ValueError("A secondary-outcome risk set contains an unobserved or invalid outcome.")


def add_andexanet(con, base):
    """Apply the G2 reversal addition to the v1.0 endpoint without modifying the G1 payload."""
    con.execute("""
        CREATE OR REPLACE TEMP TABLE p4_mcmed_andexanet_orders AS
        SELECT CSN AS visit_id, upper(trim(Order_type)) = $medications AS medication_order,
            try_cast(Order_time AS TIMESTAMP) AS order_ts
        FROM mcmed_orders
        WHERE regexp_matches(upper(trim(Procedure_name)), $andexanet)
    """, {
        "medications": cfg.G1_ESC_REVERSAL_MEDICATION_ORDER,
        "andexanet": cfg.G0_ESC_REVERSAL_PATTERNS["andexanet"],
    })
    n_labels, n_medications = con.execute("""
        SELECT count(*), count(*) FILTER (WHERE medication_order) FROM p4_mcmed_andexanet_orders
    """).fetchall()[0]
    con.register("p4_mcmed_secondary_index", base[["visit_id", "arrival_ts", "departure_ts"]])
    try:
        delta = con.execute("""
            SELECT i.visit_id,
                count(*) AS andexanet_orders_in_ed,
                count(*) FILTER (WHERE o.order_ts <= i.arrival_ts + $landmark * INTERVAL '1 minute') > 0
                    AS andexanet_before_or_at_60min,
                count(*) FILTER (WHERE o.order_ts > i.arrival_ts + $landmark * INTERVAL '1 minute') > 0
                    AS andexanet_after_60min
            FROM p4_mcmed_andexanet_orders o JOIN p4_mcmed_secondary_index i USING (visit_id)
            WHERE o.medication_order AND o.order_ts BETWEEN i.arrival_ts AND i.departure_ts
            GROUP BY i.visit_id
        """, {"landmark": cfg.ESC_LANDMARK_MINUTES}).df()
    finally:
        con.unregister("p4_mcmed_secondary_index")
    result = base.merge(delta, on="visit_id", how="left", validate="one_to_one")
    for name in ("andexanet_before_or_at_60min", "andexanet_after_60min"):
        result[name] = result[name].eq(True)
    return result, {
        "g2_named_order_records": int(n_labels),
        "g2_named_medication_order_records": int(n_medications),
        "g2_medication_orders_in_index_ed": int(delta.andexanet_orders_in_ed.sum()),
        "g2_index_visits_with_in_ed_orders": len(delta),
        "g2_index_visits_with_orders_by_landmark": int(delta.andexanet_before_or_at_60min.sum()),
        "g2_index_visits_with_orders_after_landmark": int(delta.andexanet_after_60min.sum()),
    }


def derive_landmark_exposures(con, base):
    p2 = import_module("03_p2_phenotypes")
    readings = p2.load_reads(con, "mcmed", base)
    groups = p2.ordered_groups(readings)
    rows = []
    for visit_id, sequence in groups.items():
        early = sequence[(sequence[:, 0] >= 0) & (sequence[:, 0] <= cfg.ESC_LANDMARK_MINUTES)]
        if len(early) == 0:
            continue
        first = early[0]
        after_arrival = early[early[:, 0] > 0]
        label = None
        if len(after_arrival):
            first_elevated = p2.elevated(first[1], first[2])
            last = after_arrival[-1]
            last_elevated = p2.elevated(last[1], last[2])
            label = ("P60" if last_elevated else "T60") if first_elevated else "N60"
        rows.append((visit_id, first[1] / cfg.P4_SBP_UNIT, label, first[0] > 0, len(early) == 1 and first[0] > 0))
    exposures = pd.DataFrame(rows, columns=(
        "visit_id", ARRIVAL_SBP, LANDMARK_PHENOTYPE, "first_bp_after_zero", "single_positive_time_bp",
    ))
    return base.merge(exposures, on="visit_id", how="left", validate="one_to_one")


def landmark_risk_set(base, suffix):
    old_eligibility = base[f"landmark_eligible{suffix}"].eq(True)
    old_end = f"hemostasis_level3_after_60min{suffix}"
    if base[f"landmark_eligible{suffix}"].isna().any():
        raise ValueError("The original MC-MED landmark eligibility is missing.")
    require_binary(base[old_eligibility], old_end)
    new_eligibility = old_eligibility & ~base.andexanet_before_or_at_60min
    risk = base[new_eligibility].copy()
    risk[LANDMARK_EVENT] = risk[old_end].eq(True) | risk.andexanet_after_60min
    require_binary(risk, LANDMARK_EVENT)
    phenotype = risk[risk[LANDMARK_PHENOTYPE].notna()].copy()
    arrival = risk[np.isfinite(risk[ARRIVAL_SBP])].copy()
    audit = {
        "index_visits": len(base),
        "g1_landmark_risk_set_visits": int(old_eligibility.sum()),
        "g1_landmark_events": events(base[old_eligibility], old_end),
        "g2_additional_early_event_exclusions": int((old_eligibility & ~new_eligibility).sum()),
        "g2_additional_late_events_in_revised_risk_set": int((risk.andexanet_after_60min & ~risk[old_end].eq(True)).sum()),
        "landmark_risk_set_visits": len(risk),
        "landmark_risk_set_patients": risk.patient_id.nunique(),
        "landmark_events": events(risk, LANDMARK_EVENT),
        "phenotype_observed_visits": len(phenotype),
        "phenotype_observed_events": events(phenotype, LANDMARK_EVENT),
        "phenotype_missing_visits": len(risk) - len(phenotype),
        "phenotype_missing_events": events(risk, LANDMARK_EVENT) - events(phenotype, LANDMARK_EVENT),
        "arrival_sbp_observed_visits": len(arrival),
        "arrival_sbp_observed_events": events(arrival, LANDMARK_EVENT),
        "arrival_first_reading_after_minute_zero_visits": int(arrival.first_bp_after_zero.eq(True).sum()),
        "single_positive_time_reading_visits": int(arrival.single_positive_time_bp.eq(True).sum()),
    }
    if not phenotype.visit_id.isin(arrival.visit_id).all():
        raise ValueError("A classified early phenotype lacks an arrival reading.")
    return risk, phenotype, arrival, audit


def model_rows(fit_model, *, site, analysis, endpoint, population, data, outcome, exposure=None,
               levels=(), continuous=(), adjustment="none", with_exposure=None,
               excluded_incomplete_visits=None):
    """Record common risk ratio outputs, original and fitted denominators and group events."""
    require_binary(population, outcome)
    fit = fit_model(data, outcome, exposure=exposure, levels=levels, continuous=continuous, adjustment=adjustment)
    group_exposure = exposure is not None
    reference = levels[0] if group_exposure else None
    terms = tuple((label, f"{exposure}={label}") for label in levels[1:]) if group_exposure else tuple(
        (f"Per {cfg.P4_SBP_UNIT} mmHg", term) for term in continuous
    )
    if with_exposure is None:
        with_exposure = data
    base = {
        "site": site, "analysis": analysis, "outcome": endpoint,
        "exposure": exposure or continuous[0], "adjustment": adjustment,
        "raw_risk_set_visits": len(population), "raw_risk_set_patients": population.patient_id.nunique(),
        "raw_risk_set_events": events(population, outcome),
        "excluded_incomplete_visits_after_exposure_criteria": excluded_incomplete_visits,
        "exposure_observed_visits": len(with_exposure), "exposure_observed_patients": with_exposure.patient_id.nunique(),
        "exposure_observed_events": events(with_exposure, outcome),
        "model_eligible_visits": len(data), "model_eligible_events": events(data, outcome),
        "modeled_visits": fit.n, "modeled_patients": fit.patients, "modeled_events": fit.events,
        "fit_method": fit.method, "fit_status": fit.status, "converged": fit.converged,
        "fitted_means_above_one": fit.n_gt_one, "maximum_fitted_mean": fit.max_mean,
        "omitted_terms": "; ".join(fit.omitted_terms),
    }
    rows = []
    for contrast, term in terms:
        fitted_reference = fit.group_counts.get(reference, (None, None))
        fitted_contrast = fit.group_counts.get(contrast, (None, None))
        raw_reference = data[data[exposure].eq(reference)] if group_exposure else None
        raw_contrast = data[data[exposure].eq(contrast)] if group_exposure else None
        effect = fit.estimate(term)
        effect.pop("p_value", None)
        check = fit.cross_check
        check_effect = check.estimate(term) if check is not None else {}
        rows.append(base | {
            "contrast": contrast, "reference": reference,
            "raw_contrast_visits": len(raw_contrast) if group_exposure else None,
            "raw_contrast_events": events(raw_contrast, outcome) if group_exposure else None,
            "raw_reference_visits": len(raw_reference) if group_exposure else None,
            "raw_reference_events": events(raw_reference, outcome) if group_exposure else None,
            "fitted_contrast_visits": fitted_contrast[0], "fitted_contrast_events": fitted_contrast[1],
            "fitted_reference_visits": fitted_reference[0], "fitted_reference_events": fitted_reference[1],
            "cross_check_method": check.method if check is not None else None,
            "cross_check_status": check_effect.get("status"), "cross_check_estimate": check_effect.get("estimate"),
            "cross_check_ci95_low": check_effect.get("ci95_low"), "cross_check_ci95_high": check_effect.get("ci95_high"),
        } | effect)
    return rows


def run_landmark(con, base, fit_model):
    base, reversal_audit = add_andexanet(con, base)
    base = derive_landmark_exposures(con, base)
    rows, audit_rows = [], []
    definitions = (("Primary landmark", ""), ("Any admission landmark sensitivity", "_any_admission_diagnostic"))
    for name, suffix in definitions:
        risk, phenotype, arrival, audit = landmark_risk_set(base, suffix)
        audit_rows.append({"analysis": name} | reversal_audit | audit)
        shared = {
            "site": MCMED, "analysis": name, "endpoint": HEMOSTASIS_LABEL, "population": risk,
            "outcome": LANDMARK_EVENT, "adjustment": "landmark",
        }
        rows.extend(model_rows(fit_model, data=phenotype, exposure=LANDMARK_PHENOTYPE, levels=LANDMARK_LEVELS, **shared))
        rows.extend(model_rows(fit_model, data=arrival, continuous=(ARRIVAL_SBP,), **shared))
    return rows, audit_rows


def run_inpatient(base, fit_model):
    eligible = base[base.inpatient_ent_outcome_eligible.eq(True)].copy()
    require_binary(eligible, INPATIENT_EVENT)
    require_binary(eligible, "incomplete_visit")
    with_phenotype = eligible[eligible.is_c1.eq(True) & eligible.bp_phenotype.isin(PHENOTYPE_LEVELS)]
    model = with_phenotype[with_phenotype.incomplete_visit.eq(False)].copy()
    rows = []
    for adjustment in ("none", "age"):
        rows.extend(model_rows(
            fit_model, site=MIMIC, analysis="Exploratory inpatient", endpoint=INPATIENT_LABEL,
            population=eligible, data=model, outcome=INPATIENT_EVENT, exposure="bp_phenotype",
            levels=PHENOTYPE_LEVELS, adjustment=adjustment,
            with_exposure=with_phenotype, excluded_incomplete_visits=len(with_phenotype) - len(model),
        ))
    return rows


def run_exploratory(base, fit_model):
    c2 = base[base.is_c2.eq(True)].copy()
    require_binary(c2, RETURN_EVENT)
    persistent = c2[c2.bp_phenotype.eq("P")].copy()
    subtypes = ("P-never", "P-rebound")
    subtype_data = persistent[persistent.persistent_subtype.isin(subtypes)].copy()
    rows = model_rows(
        fit_model, site=MCMED, analysis="Exploratory P subtype", endpoint=RETURN_LABEL,
        population=persistent, data=subtype_data, outcome=RETURN_EVENT,
        exposure="persistent_subtype", levels=subtypes,
    )
    measured = c2[np.isfinite(c2.hypertension_load_sbp)].copy()
    if measured.empty or (measured.hypertension_load_sbp < 0).any():
        raise ValueError("MC-MED C2 has no valid non-negative SBP burden values.")
    cutpoints = measured.hypertension_load_sbp.quantile(cfg.P4_BURDEN_TERTILE_QUANTILES).to_numpy()
    assignments = np.searchsorted(cutpoints, measured.hypertension_load_sbp, side="left")
    measured["sbp_burden_tertile"] = np.array(BURDEN_LEVELS)[assignments]
    rows.extend(model_rows(
        fit_model, site=MCMED, analysis="Exploratory SBP burden", endpoint=RETURN_LABEL,
        population=c2, data=measured, outcome=RETURN_EVENT, exposure="sbp_burden_tertile",
        levels=BURDEN_LEVELS,
    ))
    cutpoint_audit = [{
        "tertile": level,
        "strict_lower_mmhg": cutpoints[position - 1] if position > 0 else None,
        "inclusive_upper_mmhg": cutpoints[position] if position < len(cutpoints) else None,
        "visits": int((assignments == position).sum()),
        "events": events(measured[assignments == position], RETURN_EVENT),
    } for position, level in enumerate(BURDEN_LEVELS)]
    return rows, cutpoint_audit


def run(con, frame, output_dir=TABLES):
    """Write the authorized P4 aggregate secondary and exploratory results."""
    model = import_module("05_p4_models")
    p0 = import_module("01_p0_audit")
    mcmed = frame[frame.site.eq(MCMED)].copy()
    mimic = frame[frame.site.eq(MIMIC)].copy()
    if mcmed.empty or mimic.empty or mcmed.visit_id.duplicated().any() or mimic.visit_id.duplicated().any():
        raise ValueError("The P4 site inputs are empty or have duplicate index visits.")
    stage = "MC-MED landmark"
    try:
        secondary, landmark_audit = run_landmark(con, mcmed, model.fit_model)
        stage = "MIMIC inpatient"
        secondary.extend(run_inpatient(mimic, model.fit_model))
        stage = "MC-MED exploratory"
        exploratory, burden_audit = run_exploratory(mcmed, model.fit_model)
        stage = "aggregate exports"
        if output_dir is not None:
            directory = Path(output_dir)
            directory.mkdir(parents=True, exist_ok=True)
            for stem, rows in (
                ("p4_secondary_estimates", secondary), ("p4_landmark_audit", landmark_audit),
                ("p4_exploratory_estimates", exploratory), ("p4_exploratory_burden_cutpoints", burden_audit),
            ):
                fields = tuple(rows[0])
                p0.write_csv(directory / f"{stem}.csv", fields, ([row[key] for key in fields] for row in rows))
    except Exception:
        raise RuntimeError(f"P4 secondary analysis failed during {stage}.") from None
    return secondary, exploratory
