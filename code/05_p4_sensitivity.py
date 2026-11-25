"""Build the frozen-plan P4 sensitivity cohorts without fitting outcome models."""

from importlib import import_module
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd

import config as cfg


TABLES = Path("outputs/tables")
SITES = (
    ("mcmed", "MC-MED", "CSN", "MRN", "ED_dispo"),
    ("mimic", "MIMIC-IV-ED", "stay_id", "subject_id", "disposition"),
)
KEY = ["site", "visit_id"]
SCENARIOS = {
    "S11-1": "Reclassify at 160/100",
    "S11-2": "First eligible epistaxis visit",
    "S11-3": "Primary diagnosis only for index and return",
    "S11-4": "Exclude recorded antithrombotic use",
    "S11-5": "Exclude HHT",
    "S11-6": "Exclude acute class A antihypertensive use",
    "S11-7": "Include LWBS, AMA, and observation",
    "S11-9": "Last BP at least 30 minutes before departure",
    "S11-12": "Medication list required in MC-MED",
    "S11-13": "Count same-day PMH hypertension in MC-MED",
    "S11-16": "Medication list required in both centers",
    "S11-17": "Ignore prior ED hypertension diagnosis in MIMIC-IV-ED",
}
REBUILDS = {"S11-2": "first", "S11-3": "dx"}


def prepare_index_sets(con):
    core = import_module("02_p1_cohort_core")
    core.prepare_candidates(con)
    core.prepare_episode_cohorts(con)
    for short, _, visit, patient, _ in SITES:
        con.execute(f"""
            CREATE OR REPLACE TEMP TABLE p4_s11_{short}_history AS
            SELECT {visit} AS visit_id, {patient} AS patient_id, arrival_ts, dx_match
            FROM p1_{short}_adult_episodes
        """)
        con.execute(f"""
            CREATE OR REPLACE TEMP TABLE p4_s11_{short}_first AS
            SELECT * FROM p1_{short}_eligible_visits WHERE is_first_eligible_visit
        """)
        con.execute(f"""
            CREATE OR REPLACE TEMP TABLE p1_{short}_candidates AS
            SELECT * FROM p1_{short}_candidates WHERE dx_match
        """)
    core.prepare_episode_cohorts(con)
    for short, _, _, _, _ in SITES:
        con.execute(f"""
            CREATE OR REPLACE TEMP TABLE p4_s11_{short}_dx AS SELECT * FROM p1_{short}_index
        """)


def missing_index_rows(con, p4_index):
    con.register("p4_s11_current_keys", p4_index[KEY])
    try:
        for short, site, visit, _, _ in SITES:
            con.execute(f"""
                CREATE OR REPLACE TEMP TABLE p1_{short}_index AS
                WITH either AS (
                    SELECT *, 0 AS priority FROM p4_s11_{short}_first
                    UNION ALL
                    SELECT *, 1 AS priority FROM p4_s11_{short}_dx
                )
                SELECT a.* EXCLUDE (priority)
                FROM either AS a
                ANTI JOIN p4_s11_current_keys AS b ON b.site = $site AND b.visit_id = a.{visit}
                QUALIFY row_number() OVER (PARTITION BY a.{visit} ORDER BY a.priority) = 1
            """, {"site": site}).fetchall()
    finally:
        con.unregister("p4_s11_current_keys")
    return sum(
        con.execute(f"SELECT count(*) FROM p1_{short}_index").fetchall()[0][0]
        for short, _, _, _, _ in SITES
    )


def stage_extra_covariates(con):
    covariates = import_module("02_p1_covariates")
    d = import_module("02_p1_dictionaries")
    for source, subset, id_field, keys, lookup, fields, classify in (
        ("mcmed_meds", "mcmed", "MRN", ("Med_class", "Med_subclass"),
         "mcmed_meds", d.HOME_FIELDS, d.classify_home_acap),
        ("mimic_medrecon", "mimic", "stay_id", ("etcdescription",),
         "mimic_medrecon", d.HOME_FIELDS, lambda name: d.classify_home_acap(name, mimic=True)),
        ("mcmed_orders", "mcmed", "CSN", ("Order_type", "Procedure_name"),
         "mcmed_orders", d.ED_FIELDS + d.ESCALATION_FIELDS,
         lambda kind, name: d.classify_aht_ed(name, kind) | d.classify_escalation(kind, name)),
        ("mimic_pyxis", "mimic", "stay_id", ("name",),
         "mimic_pyxis", d.ED_FIELDS, lambda name: d.classify_aht_ed(name, mimic=True)),
    ):
        source_subset = (
            f"(SELECT {', '.join('r.' + col for col in keys)} FROM {source} AS r "
            f"SEMI JOIN p1_{subset}_index AS i USING ({id_field}))"
        )
        d.map_labels(con, source_subset, keys, f"p1_dict_{lookup}", fields, classify, None)
    covariates.create_mcmed_history(con)
    covariates.create_mimic_history(con)
    covariates.create_mcmed_g1_history(con)
    covariates.create_mimic_g1_history(con)
    covariates.create_order_events(con)
    pairs = [(label, level) for label, levels in cfg.TRIAGE_ACUITY_GROUPS for level in levels]
    con.execute(
        "CREATE OR REPLACE TEMP TABLE p4_s11_acuity AS SELECT * FROM (VALUES "
        + ", ".join(["(?, ?)"] * len(pairs)) + ") g(triage_acuity_group, acuity_level)",
        [item for pair in pairs for item in pair],
    )


def extra_baseline(con):
    return con.execute("""
        SELECT 'MC-MED' AS site, i.CSN AS visit_id, i.MRN AS patient_id,
            i.arrival_ts, i.departure_ts, i.ED_dispo AS disposition, i.age_years,
            i.dx_match, i.cc_match, i.is_first_eligible_visit, i.prior_epistaxis_365d,
            i.incomplete_visit, i.is_c1, i.is_c2, i.bp_readings, i.bp_span_minutes,
            CASE upper(trim(i.Gender)) WHEN 'F' THEN 'Female' WHEN 'FEMALE' THEN 'Female'
                WHEN 'M' THEN 'Male' WHEN 'MALE' THEN 'Male' END AS sex,
            a.triage_acuity_group, m.has_medication_list, h.hht, h.ckd, h.liver_disease,
            m.anticoagulant_use, m.antiplatelet_use,
            m.anticoagulant_use OR m.antiplatelet_use AS any_antithrombotic_use,
            m.prior_htn_medication, d.prior_htn_diagnosis,
            d.prior_htn_diagnosis OR m.prior_htn_medication AS prior_htn,
            NULL::BOOLEAN AS prior_htn_hosp_diagnosis, NULL::BOOLEAN AS prior_htn_ed_diagnosis,
            d.prior_htn_diagnosis OR d.prior_htn_same_day OR m.prior_htn_medication
                AS prior_htn_same_day_sensitivity,
            coalesce(e.ed_acute_aht, false) AS ed_acute_aht, e.ed_first_acute_aht_minutes
        FROM p1_mcmed_index AS i
        JOIN p1_mcmed_pmh_history AS d ON d.visit_id = i.CSN
        JOIN p1_mcmed_home_meds AS m ON m.visit_id = i.CSN
        JOIN p2_mcmed_g1_history AS h ON h.visit_id = i.CSN
        LEFT JOIN p1_mcmed_ed_meds AS e ON e.visit_id = i.CSN
        LEFT JOIN p4_s11_acuity AS a ON a.acuity_level = coalesce(
            try_cast(i.Triage_acuity AS DOUBLE),
            map_extract_value(map($acuity_codes, $acuity_levels), upper(trim(i.Triage_acuity)))
        )
        UNION ALL
        SELECT 'MIMIC-IV-ED', i.stay_id, i.subject_id,
            i.arrival_ts, i.departure_ts, i.disposition, i.age_years,
            i.dx_match, i.cc_match, i.is_first_eligible_visit, i.prior_epistaxis_365d,
            i.incomplete_visit, i.is_c1, i.is_c2, i.bp_readings, i.bp_span_minutes,
            CASE upper(trim(i.gender)) WHEN 'F' THEN 'Female' WHEN 'FEMALE' THEN 'Female'
                WHEN 'M' THEN 'Male' WHEN 'MALE' THEN 'Male' END,
            a.triage_acuity_group, m.has_medication_list, h.hht, h.ckd, h.liver_disease,
            m.anticoagulant_use, m.antiplatelet_use, m.anticoagulant_use OR m.antiplatelet_use,
            m.prior_htn_medication, h.prior_htn_hosp_diagnosis OR h.prior_htn_ed_diagnosis,
            h.prior_htn_hosp_diagnosis OR h.prior_htn_ed_diagnosis OR m.prior_htn_medication,
            h.prior_htn_hosp_diagnosis, h.prior_htn_ed_diagnosis, NULL::BOOLEAN,
            coalesce(e.ed_acute_aht, false), e.ed_first_acute_aht_minutes
        FROM p1_mimic_index AS i
        JOIN p1_mimic_home_meds AS m ON m.visit_id = i.stay_id
        JOIN p2_mimic_g1_history AS h ON h.visit_id = i.stay_id
        LEFT JOIN p1_mimic_ed_meds AS e ON e.visit_id = i.stay_id
        LEFT JOIN (
            SELECT t.stay_id, min(try_cast(t.acuity AS DOUBLE)) AS acuity_level
            FROM mimic_triage AS t
            SEMI JOIN p1_mimic_index AS x USING (stay_id)
            WHERE isfinite(try_cast(t.acuity AS DOUBLE))
            GROUP BY t.stay_id
            HAVING count(DISTINCT try_cast(t.acuity AS DOUBLE)) = 1
        ) AS t USING (stay_id)
        LEFT JOIN p4_s11_acuity AS a ON t.acuity_level = a.acuity_level
    """, {"acuity_codes": list(cfg.MCMED_TRIAGE_ACUITY_CODES),
          "acuity_levels": list(cfg.MCMED_TRIAGE_ACUITY_CODES.values())}).df()


def enrich_index_sets(con, p4_index):
    expected = missing_index_rows(con, p4_index)
    if not expected:
        return p4_index
    core = import_module("02_p1_cohort_core")
    phenotype = import_module("03_p2_phenotypes")
    core.prepare_bp_slots(con)
    core.deduplicate_bp(con)
    stage_extra_covariates(con)
    extras = extra_baseline(con)
    if len(extras) != expected or extras.duplicated(KEY).any():
        raise ValueError("Rebuilt index covariates do not match the added visits.")
    results = []
    for short, site, visit, _, _ in SITES:
        reads = con.execute(f"""
            SELECT b.visit_id, b.source, b.minutes_from_arrival, b.sbp, b.dbp
            FROM p1_bp_{short} AS b
            JOIN p1_{short}_index AS i ON i.{visit} = b.visit_id
            WHERE i.is_c1
        """).df()
        ordered = phenotype.ordered_groups(reads)
        selected = extras[extras["site"].eq(site)]
        expected_c1 = set(selected.loc[selected["is_c1"], "visit_id"])
        if set(ordered) != expected_c1:
            raise ValueError("Rebuilt BP records do not match the added C1 visits.")
        if ordered:
            metrics = pd.DataFrame(
                phenotype.summarize_visit(key, values, short == "mcmed") for key, values in ordered.items()
            )
            selected = selected.merge(metrics, on="visit_id", how="left", validate="one_to_one")
        results.append(selected)
    return pd.concat([p4_index, *results], ignore_index=True)


def rebuilt_selection(con, base, kind):
    rows = []
    for short, site, visit, _, _ in SITES:
        selected = con.execute(f"""
            SELECT $site AS site, {visit} AS visit_id, is_first_eligible_visit, prior_epistaxis_365d
            FROM p4_s11_{short}_{kind}
        """, {"site": site}).df()
        shared = [col for col in selected if col not in KEY and col in base]
        merged = selected.merge(base.drop(columns=shared), on=KEY, how="left",
                                validate="one_to_one", indicator=True)
        if not merged["_merge"].eq("both").all():
            raise ValueError("A rebuilt visit is missing from the assembled sensitivity inputs.")
        rows.append(merged.drop(columns="_merge"))
    return pd.concat(rows, ignore_index=True)


def recompute_returns(con, cohort, diagnosis_only=False):
    results = []
    for short, site, _, _, _ in SITES:
        sample = cohort[cohort["site"].eq(site)].copy()
        if not len(sample):
            results.append(sample)
            continue
        con.register("p4_s11_outcome_index", sample[["visit_id", "patient_id", "departure_ts"]])
        try:
            rows = con.execute(f"""
                WITH next_valid_return AS (
                    SELECT i.visit_id, i.departure_ts, min(r.arrival_ts) AS return_ts
                    FROM p4_s11_outcome_index AS i
                    LEFT JOIN p4_s11_{short}_history AS r
                      ON r.patient_id = i.patient_id AND r.visit_id <> i.visit_id
                     AND (NOT $dx_only OR r.dx_match)
                     AND r.arrival_ts > i.departure_ts
                     AND r.arrival_ts <= i.departure_ts + $long_hours * INTERVAL '1 hour'
                    GROUP BY i.visit_id, i.departure_ts
                )
                SELECT visit_id, true AS return_outcome_eligible,
                       coalesce(return_ts <= departure_ts + $short_hours * INTERVAL '1 hour', false)
                           AS return_7d_event,
                       return_ts IS NOT NULL AS return_30d_event,
                       epoch(return_ts - departure_ts) / 60.0 AS return_minutes_within_30d
                FROM next_valid_return
            """, {"dx_only": diagnosis_only, "long_hours": cfg.RETURN_30D_HOURS,
                  "short_hours": cfg.RETURN_7D_HOURS}).df()
        finally:
            con.unregister("p4_s11_outcome_index")
        replace = [col for col in rows if col != "visit_id" and col in sample]
        results.append(sample.drop(columns=replace).merge(rows, on="visit_id", validate="one_to_one"))
    return pd.concat(results, ignore_index=True)


def classify_endpoints(frame, systolic, diastolic):
    first_elevated = (frame["sbp_first"] >= systolic) | (frame["dbp_first"] >= diastolic)
    end_normal = (frame["sbp_last"] < systolic) & (frame["dbp_last"] < diastolic)
    return frame.assign(bp_phenotype=np.where(first_elevated, np.where(end_normal, "T", "P"), "N"))


def earlier_endpoints(con, c2):
    phenotype = import_module("03_p2_phenotypes")
    results = []
    for short, site, _, _, _ in SITES:
        sample = c2[c2["site"].eq(site)]
        con.register("p4_s11_earlier_index", sample[["visit_id", "arrival_ts", "departure_ts"]])
        try:
            reads = con.execute(f"""
                SELECT b.visit_id, b.source, b.minutes_from_arrival, b.sbp, b.dbp
                FROM read_parquet('derived/bp_long_{short}.parquet') AS b
                JOIN p4_s11_earlier_index AS i USING (visit_id)
                WHERE cast(round(b.minutes_from_arrival * 60000000.0) AS BIGINT) <=
                    date_diff('microsecond', i.arrival_ts, i.departure_ts - ? * INTERVAL '1 minute')
            """, [cfg.P4_PREDEPARTURE_MINUTES]).df()
        finally:
            con.unregister("p4_s11_earlier_index")
        rows = []
        for key, order in phenotype.ordered_groups(reads).items():
            time, sbp, dbp = order[np.searchsorted(order[:, 0], order[-1, 0])]
            rows.append((key, time, sbp, dbp))
        fields = ["visit_id", "time_last_minutes", "sbp_last", "dbp_last"]
        alternate = pd.DataFrame(rows, columns=fields)
        if len(alternate):
            results.append(sample.drop(columns=fields[1:]).merge(alternate, on="visit_id", validate="one_to_one"))
        else:
            results.append(sample.iloc[:0].copy())
    return classify_endpoints(pd.concat(results, ignore_index=True), cfg.BP_ELEVATED_SBP, cfg.BP_ELEVATED_DBP)


def validate_frames(main, scenarios):
    required = ["visit_id", "patient_id", "site", "return_7d_event", "return_30d_event", "bp_phenotype",
                "is_c1", "sbp_first", "dbp_first", "sbp_last", "dbp_last", "prior_htn"]
    for frame in (main, *scenarios.values()):
        if (
            frame[required].isna().any().any() or frame.duplicated(KEY).any()
            or not frame["is_c1"].all()
            or not frame["bp_phenotype"].isin(("N", "T", "P")).all()
            or (frame["return_7d_event"].eq(True) & frame["return_30d_event"].ne(True)).any()
        ):
            raise ValueError("Invalid P4 sensitivity model cohort.")
    if scenarios["S11-2"].duplicated(["site", "patient_id"]).any():
        raise ValueError("First-eligible-visit sensitivity includes a repeated patient.")
    if not scenarios["S11-3"]["dx_match"].all():
        raise ValueError("Primary-diagnosis sensitivity includes a chief-complaint-only visit.")
    main_return = main.set_index(KEY)["return_7d_event"]
    for key in ("S11-2", "S11-3", "S11-7"):
        alternative = scenarios[key].set_index(KEY)["return_7d_event"]
        frozen, rebuilt = main_return.align(alternative, join="inner")
        if key == "S11-3":
            invalid = (rebuilt & ~frozen).any()
        else:
            invalid = (rebuilt != frozen).any()
        if invalid:
            raise ValueError("Rebuilt return flags disagree unexpectedly with the frozen main endpoint.")


def write_flow(tables, main, scenarios, stages):
    rows = []
    for key, frame in scenarios.items():
        source, with_c1, with_disposition = stages[key]
        populations = (("MC-MED",), ("MIMIC-IV-ED",), ("MC-MED", "MIMIC-IV-ED"))
        for sites in populations:
            subset = frame[frame["site"].isin(sites)]
            base = main[main["site"].isin(sites)]
            merged = base.merge(subset, on=KEY, suffixes=("_main", "_sensitivity"), validate="one_to_one")
            rows.append({
                "analysis": key, "definition": SCENARIOS[key],
                "site": sites[0] if len(sites) == 1 else "Pooled",
                "index_source_visits": source["site"].isin(sites).sum(),
                "after_c1_visits": with_c1["site"].isin(sites).sum(),
                "after_disposition_visits": with_disposition["site"].isin(sites).sum(),
                "analysis_visits": len(subset),
                "analysis_patients": subset[["site", "patient_id"]].drop_duplicates().shape[0],
                "return_7d_events": subset["return_7d_event"].sum(),
                "phenotype_n": subset["bp_phenotype"].eq("N").sum(),
                "phenotype_t": subset["bp_phenotype"].eq("T").sum(),
                "phenotype_p": subset["bp_phenotype"].eq("P").sum(),
                "added_visits_vs_main": len(subset) - len(merged),
                "removed_visits_vs_main": len(base) - len(merged),
                "changed_phenotype_on_shared_visits": (
                    merged["bp_phenotype_main"].ne(merged["bp_phenotype_sensitivity"]).sum()
                ),
                "changed_prior_htn_on_shared_visits": (
                    merged["prior_htn_main"].ne(merged["prior_htn_sensitivity"]).sum()
                ),
                "changed_7d_outcome_on_shared_visits": (
                    merged["return_7d_event_main"].ne(merged["return_7d_event_sensitivity"]).sum()
                ),
            })
    tables.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(tables / "p4_sensitivity_flow.csv", index=False)


def build(con, p4_index, tables=TABLES):
    """Return pooled-ready frames keyed by S11 number and write the aggregate cohort flow."""
    try:
        main = p4_index[p4_index["is_c2"]].copy()
        mcmed = main["site"].eq("MC-MED")
        mimic = main["site"].eq("MIMIC-IV-ED")
        sens13 = main.copy()
        sens13.loc[mcmed, "prior_htn"] = main.loc[mcmed, "prior_htn_same_day_sensitivity"]
        sens17 = main.copy()
        sens17.loc[mimic, "prior_htn"] = (
            main.loc[mimic, "prior_htn_medication"] | main.loc[mimic, "prior_htn_hosp_diagnosis"]
        )
        variants = {
            "S11-1": classify_endpoints(main, cfg.P4_S11_HIGH_SBP, cfg.P4_S11_HIGH_DBP),
            "S11-4": main[~main["any_antithrombotic_use"]],
            "S11-5": main[~main["hht"]],
            "S11-6": main[~main["ed_acute_aht"]],
            "S11-9": earlier_endpoints(con, main),
            "S11-12": main[~mcmed | main["has_medication_list"]],
            "S11-13": sens13,
            "S11-16": main[main["has_medication_list"]],
            "S11-17": sens17,
        }
        stages = {key: (p4_index, p4_index[p4_index["is_c1"]], main) for key in variants}
        prepare_index_sets(con)
        base = enrich_index_sets(con, p4_index)
        for key, kind in REBUILDS.items():
            source = rebuilt_selection(con, base, kind)
            selected = source[source["is_c2"]]
            variants[key] = recompute_returns(con, selected, diagnosis_only=kind == "dx")
            stages[key] = (source, source[source["is_c1"]], selected)
        all_c1 = p4_index[p4_index["is_c1"]]
        allowed = all_c1["is_c2"] | all_c1["disposition"].str.strip().str.upper().isin(cfg.P4_S11_EXTRA_DISPOSITIONS)
        valid_time = all_c1["departure_ts"].gt(all_c1["arrival_ts"])
        combined = all_c1[allowed & valid_time]
        variants["S11-7"] = recompute_returns(con, combined)
        stages["S11-7"] = (p4_index, all_c1, combined)
        variants = {key: variants[key] for key in SCENARIOS}
        validate_frames(main, variants)
        write_flow(tables, main, variants, stages)
        return variants
    except duckdb.Error as exc:
        raise RuntimeError(f"P4 sensitivity cohort construction failed ({type(exc).__name__}).") from None
