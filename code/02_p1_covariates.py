"""Build index-level covariates and timed event flags from the approved G0 mappings."""

import csv
from pathlib import Path

import config as cfg


def diagnosis_patterns():
    return {
        "htn10": cfg.MCMED_HTN_ICD10_PATTERN,
        "htn9": cfg.MCMED_HTN_ICD9_PATTERN,
        "ckd10": cfg.CKD_DX10_PATTERN,
        "ckd9": cfg.CKD_DX9_PATTERN,
        "liver10": cfg.LIVER_DX10_PATTERN,
        "liver9": cfg.LIVER_DX9_PATTERN,
        "hht10": cfg.HHT_DX10_PATTERN,
        "hht9": cfg.HHT_DX9_PATTERN,
    }


def create_mcmed_history(con):
    con.execute(
        """
        CREATE OR REPLACE TEMP TABLE p1_mcmed_pmh_history AS
        WITH relevant_people AS (SELECT DISTINCT MRN FROM p1_mcmed_index),
        coded AS (
            SELECT p.MRN, try_cast(p.Noted_date AS DATE) AS noted_date,
                   regexp_matches(trim(p.Code), $htn10) OR regexp_matches(trim(p.Code), $htn9) AS htn,
                   regexp_matches(trim(p.Code), $ckd10) OR regexp_matches(trim(p.Code), $ckd9) AS ckd,
                   regexp_matches(trim(p.Code), $liver10) OR regexp_matches(trim(p.Code), $liver9) AS liver_disease,
                   regexp_matches(trim(p.Code), $hht10) OR regexp_matches(trim(p.Code), $hht9) AS hht
            FROM mcmed_pmh p JOIN relevant_people USING (MRN)
        )
        SELECT i.CSN AS visit_id,
               count(*) FILTER (WHERE p.htn AND p.noted_date < i.arrival_ts::DATE) > 0 AS prior_htn_diagnosis,
               count(*) FILTER (WHERE p.htn AND p.noted_date = i.arrival_ts::DATE) > 0 AS prior_htn_same_day,
               count(*) FILTER (WHERE p.ckd AND p.noted_date < i.arrival_ts::DATE) > 0 AS ckd,
               count(*) FILTER (WHERE p.liver_disease AND p.noted_date < i.arrival_ts::DATE) > 0 AS liver_disease,
               count(*) FILTER (WHERE p.hht AND p.noted_date < i.arrival_ts::DATE) > 0 AS hht
        FROM p1_mcmed_index i LEFT JOIN coded p USING (MRN)
        GROUP BY i.CSN
        """,
        diagnosis_patterns(),
    )
    con.execute(
        """
        CREATE OR REPLACE TEMP TABLE p1_mcmed_home_meds AS
        WITH relevant_people AS (SELECT DISTINCT MRN FROM p1_mcmed_index),
        meds AS (
            SELECT m.MRN, try_cast(m.Entry_date AS DATE) AS entry_date,
                   try_cast(m.Start_date AS DATE) AS start_date, try_cast(m.End_date AS DATE) AS end_date,
                   nullif(trim(m.Start_date), '') IS NULL AS start_missing,
                   nullif(trim(m.End_date), '') IS NULL AS end_missing,
                   d.aht_home, d.aht_home_sensitivity_only, d.anticoagulant, d.antiplatelet
            FROM mcmed_meds m JOIN relevant_people USING (MRN)
            LEFT JOIN p1_dict_mcmed_meds d ON m.Med_class IS NOT DISTINCT FROM d.Med_class
                                     AND m.Med_subclass IS NOT DISTINCT FROM d.Med_subclass
        ), dates AS (
            SELECT i.CSN, m.MRN,
                   m.entry_date < i.arrival_ts::DATE
                       AND (m.start_missing OR m.start_date <= i.arrival_ts::DATE)
                       AND (m.end_missing OR m.end_date >= i.arrival_ts::DATE) AS active_at_index,
                   m.aht_home, m.aht_home_sensitivity_only, m.anticoagulant, m.antiplatelet
            FROM p1_mcmed_index i LEFT JOIN meds m USING (MRN)
        )
        SELECT CSN AS visit_id, count(MRN) > 0 AS has_medication_list,
               count(*) FILTER (WHERE active_at_index AND aht_home) > 0 AS prior_htn_medication,
               count(*) FILTER (WHERE active_at_index AND aht_home_sensitivity_only) > 0 AS prior_htn_home_sensitivity_only,
               count(*) FILTER (WHERE active_at_index AND anticoagulant) > 0 AS anticoagulant_use,
               count(*) FILTER (WHERE active_at_index AND antiplatelet) > 0 AS antiplatelet_use
        FROM dates GROUP BY CSN
        """
    )


def create_mimic_history(con):
    # For each condition an earliest eligible diagnosis time is sufficient for every later index.
    # A previous ED stay has discharged by arrival; hosp uses the SAP's admission timestamp.
    con.execute(
        """
        CREATE OR REPLACE TEMP TABLE p1_mimic_diagnosis_history AS
        WITH relevant_people AS (SELECT DISTINCT subject_id FROM p1_mimic_index),
        diagnoses AS (
            SELECT ed.subject_id, try_cast(ed.outtime AS TIMESTAMP) AS code_ts,
                   'ed' AS source, trim(d.icd_code) AS code, trim(d.icd_version) AS code_version
            FROM relevant_people p
            JOIN mimic_edstays ed USING (subject_id)
            JOIN mimic_diagnosis d ON d.stay_id = ed.stay_id AND d.subject_id = ed.subject_id
            WHERE try_cast(ed.outtime AS TIMESTAMP) > try_cast(ed.intime AS TIMESTAMP)
            UNION ALL
            SELECT a.subject_id, try_cast(a.admittime AS TIMESTAMP),
                   'hosp', trim(d.icd_code), trim(d.icd_version)
            FROM relevant_people p
            JOIN mimic_hosp_admissions a USING (subject_id)
            JOIN mimic_hosp_diagnoses_icd d ON d.hadm_id = a.hadm_id AND d.subject_id = a.subject_id
        ), coded AS (
            SELECT subject_id, code_ts, source,
                   CASE code_version WHEN '10' THEN regexp_matches(code, $htn10)
                                     WHEN '9' THEN regexp_matches(code, $htn9) END AS htn,
                   CASE code_version WHEN '10' THEN regexp_matches(code, $ckd10)
                                     WHEN '9' THEN regexp_matches(code, $ckd9) END AS ckd,
                   CASE code_version WHEN '10' THEN regexp_matches(code, $liver10)
                                     WHEN '9' THEN regexp_matches(code, $liver9) END AS liver_disease,
                   CASE code_version WHEN '10' THEN regexp_matches(code, $hht10)
                                     WHEN '9' THEN regexp_matches(code, $hht9) END AS hht
            FROM diagnoses
        ), first_record AS (
            SELECT subject_id,
                   min(code_ts) FILTER (WHERE htn AND source = 'ed') AS ed_htn,
                   min(code_ts) FILTER (WHERE htn AND source = 'hosp') AS hosp_htn,
                   min(code_ts) FILTER (WHERE ckd AND source = 'ed') AS ed_ckd,
                   min(code_ts) FILTER (WHERE ckd AND source = 'hosp') AS hosp_ckd,
                   min(code_ts) FILTER (WHERE liver_disease AND source = 'ed') AS ed_liver,
                   min(code_ts) FILTER (WHERE liver_disease AND source = 'hosp') AS hosp_liver,
                   min(code_ts) FILTER (WHERE hht AND source = 'ed') AS ed_hht,
                   min(code_ts) FILTER (WHERE hht AND source = 'hosp') AS hosp_hht
            FROM coded GROUP BY subject_id
        )
        SELECT i.stay_id AS visit_id,
               coalesce(ed_htn <= i.arrival_ts, false) OR coalesce(hosp_htn < i.arrival_ts, false) AS prior_htn_diagnosis,
               coalesce(ed_ckd <= i.arrival_ts, false) OR coalesce(hosp_ckd < i.arrival_ts, false) AS ckd,
               coalesce(ed_liver <= i.arrival_ts, false) OR coalesce(hosp_liver < i.arrival_ts, false) AS liver_disease,
               coalesce(ed_hht <= i.arrival_ts, false) OR coalesce(hosp_hht < i.arrival_ts, false) AS hht
        FROM p1_mimic_index i LEFT JOIN first_record USING (subject_id)
        """,
        diagnosis_patterns(),
    )
    con.execute(
        """
        CREATE OR REPLACE TEMP TABLE p1_mimic_home_meds AS
        SELECT i.stay_id AS visit_id, count(m.stay_id) > 0 AS has_medication_list,
               count(*) FILTER (WHERE d.aht_home) > 0 AS prior_htn_medication,
               count(*) FILTER (WHERE d.aht_home_sensitivity_only) > 0 AS prior_htn_home_sensitivity_only,
               count(*) FILTER (WHERE d.anticoagulant) > 0 AS anticoagulant_use,
               count(*) FILTER (WHERE d.antiplatelet) > 0 AS antiplatelet_use
        FROM p1_mimic_index i LEFT JOIN mimic_medrecon m USING (stay_id)
        LEFT JOIN p1_dict_mimic_medrecon d ON m.etcdescription = d.etcdescription
        GROUP BY i.stay_id
        """
    )


def create_mcmed_g1_history(con):
    con.execute(
        """
        CREATE OR REPLACE TEMP TABLE p2_mcmed_g1_history AS
        WITH people AS (SELECT DISTINCT MRN FROM p1_mcmed_index),
        coded_pmh AS (
            SELECT p.MRN, try_cast(p.Noted_date AS DATE) AS noted_date,
                   regexp_matches(trim(p.Code), $ckd10) OR regexp_matches(trim(p.Code), $ckd9) AS ckd,
                   regexp_matches(trim(p.Code), $liver10) OR regexp_matches(trim(p.Code), $liver9) AS liver_disease,
                   regexp_matches(trim(p.Code), $hht10) OR regexp_matches(trim(p.Code), $hht9) AS hht
            FROM mcmed_pmh p JOIN people USING (MRN)
        ), history AS (
            SELECT i.CSN, count(*) FILTER (WHERE p.ckd AND p.noted_date < i.arrival_ts::DATE) > 0 AS ckd_before,
                   count(*) FILTER (WHERE p.liver_disease AND p.noted_date < i.arrival_ts::DATE) > 0 AS liver_before,
                   count(*) FILTER (WHERE p.hht AND p.noted_date < i.arrival_ts::DATE) > 0 AS hht_before,
                   count(*) FILTER (WHERE p.ckd AND p.noted_date <= i.arrival_ts::DATE) > 0 AS ckd_through_date,
                   count(*) FILTER (WHERE p.liver_disease AND p.noted_date <= i.arrival_ts::DATE) > 0 AS liver_through_date,
                   count(*) FILTER (WHERE p.hht) > 0 AS hht_any_pmh
            FROM p1_mcmed_index i LEFT JOIN coded_pmh p USING (MRN)
            GROUP BY i.CSN
        ), hht_ed AS (
            SELECT v.MRN FROM mcmed_visits v JOIN people USING (MRN)
            WHERE regexp_matches(trim(v.Dx_ICD10), $hht10) OR regexp_matches(trim(v.Dx_ICD9), $hht9)
            GROUP BY v.MRN
        )
        SELECT i.CSN AS visit_id, h.ckd_before, h.liver_before, h.hht_before,
               h.ckd_through_date OR coalesce(regexp_matches(trim(i.Dx_ICD10), $ckd10), false)
                   OR coalesce(regexp_matches(trim(i.Dx_ICD9), $ckd9), false) AS ckd,
               h.liver_through_date OR coalesce(regexp_matches(trim(i.Dx_ICD10), $liver10), false)
                   OR coalesce(regexp_matches(trim(i.Dx_ICD9), $liver9), false) AS liver_disease,
               h.hht_any_pmh OR e.MRN IS NOT NULL AS hht
        FROM p1_mcmed_index i JOIN history h USING (CSN) LEFT JOIN hht_ed e USING (MRN)
        """,
        {name: pattern for name, pattern in diagnosis_patterns().items() if not name.startswith("htn")},
    )


def create_mimic_g1_history(con):
    con.execute(
        """
        CREATE OR REPLACE TEMP TABLE p2_mimic_g1_history AS
        WITH people AS (SELECT DISTINCT subject_id FROM p1_mimic_index),
        diagnoses AS (
            SELECT d.subject_id, d.stay_id,
                   CASE WHEN try_cast(e.outtime AS TIMESTAMP) > try_cast(e.intime AS TIMESTAMP)
                        THEN try_cast(e.outtime AS TIMESTAMP) END AS code_ts,
                   'ed' AS source, trim(d.icd_code) AS code, trim(d.icd_version) AS code_version
            FROM mimic_diagnosis d JOIN people USING (subject_id)
            LEFT JOIN mimic_edstays e ON e.subject_id = d.subject_id AND e.stay_id = d.stay_id
            UNION ALL
            SELECT d.subject_id, NULL::VARCHAR, try_cast(a.admittime AS TIMESTAMP),
                   'hosp', trim(d.icd_code), trim(d.icd_version)
            FROM mimic_hosp_diagnoses_icd d JOIN people USING (subject_id)
            LEFT JOIN mimic_hosp_admissions a ON a.subject_id = d.subject_id AND a.hadm_id = d.hadm_id
        ), coded AS (
            SELECT subject_id, stay_id, code_ts, source,
                   CASE code_version WHEN '10' THEN regexp_matches(code, $htn10)
                                     WHEN '9' THEN regexp_matches(code, $htn9) END AS htn,
                   CASE code_version WHEN '10' THEN regexp_matches(code, $ckd10)
                                     WHEN '9' THEN regexp_matches(code, $ckd9) END AS ckd,
                   CASE code_version WHEN '10' THEN regexp_matches(code, $liver10)
                                     WHEN '9' THEN regexp_matches(code, $liver9) END AS liver_disease,
                   CASE code_version WHEN '10' THEN regexp_matches(code, $hht10)
                                     WHEN '9' THEN regexp_matches(code, $hht9) END AS hht
            FROM diagnoses
        ), first_record AS (
            SELECT subject_id,
                   min(code_ts) FILTER (WHERE htn AND source = 'ed') AS ed_htn,
                   min(code_ts) FILTER (WHERE htn AND source = 'hosp') AS hosp_htn,
                   min(code_ts) FILTER (WHERE ckd AND source = 'ed') AS ed_ckd,
                   min(code_ts) FILTER (WHERE ckd AND source = 'hosp') AS hosp_ckd,
                   min(code_ts) FILTER (WHERE liver_disease AND source = 'ed') AS ed_liver,
                   min(code_ts) FILTER (WHERE liver_disease AND source = 'hosp') AS hosp_liver,
                   min(code_ts) FILTER (WHERE hht AND source = 'ed') AS ed_hht,
                   min(code_ts) FILTER (WHERE hht AND source = 'hosp') AS hosp_hht,
                   count(*) FILTER (WHERE hht) > 0 AS hht_any_time
            FROM coded GROUP BY subject_id
        ), at_index AS (
            SELECT stay_id, count(*) FILTER (WHERE ckd) > 0 AS ckd_any_position,
                   count(*) FILTER (WHERE liver_disease) > 0 AS liver_any_position
            FROM coded WHERE source = 'ed'
            GROUP BY stay_id
        ), flags AS (
            SELECT i.stay_id AS visit_id,
                   coalesce(hosp_htn < i.arrival_ts, false) AS prior_htn_hosp_diagnosis,
                   coalesce(ed_htn <= i.arrival_ts, false) AS prior_htn_ed_diagnosis,
                   coalesce(ed_ckd <= i.arrival_ts, false) OR coalesce(hosp_ckd < i.arrival_ts, false) AS ckd_before,
                   coalesce(ed_liver <= i.arrival_ts, false) OR coalesce(hosp_liver < i.arrival_ts, false) AS liver_before,
                   coalesce(ed_hht <= i.arrival_ts, false) OR coalesce(hosp_hht < i.arrival_ts, false) AS hht_before,
                   coalesce(x.ckd_any_position, false) AS index_ckd,
                   coalesce(x.liver_any_position, false) AS index_liver,
                   coalesce(hht_any_time, false) AS hht
            FROM p1_mimic_index i LEFT JOIN first_record USING (subject_id) LEFT JOIN at_index x USING (stay_id)
        )
        SELECT visit_id, prior_htn_hosp_diagnosis, prior_htn_ed_diagnosis, ckd_before, liver_before, hht_before,
               ckd_before OR index_ckd AS ckd, liver_before OR index_liver AS liver_disease, hht
        FROM flags
        """,
        diagnosis_patterns(),
    )


def create_triage(con):
    con.execute(
        """
        CREATE OR REPLACE TEMP TABLE p1_mcmed_triage_pain AS
        WITH parsed AS (
            SELECT i.CSN AS visit_id, try_cast(n."Time" AS TIMESTAMP) AS pain_time,
                   try_cast(n.Value AS DOUBLE) AS pain, i.arrival_ts, i.departure_ts
            FROM mcmed_numerics n JOIN p1_mcmed_index i USING (CSN)
            WHERE list_contains($pain_measures, upper(trim(n.Measure)))
        ), first_clear AS (
            SELECT visit_id, pain_time, arrival_ts, min(pain) AS pain
            FROM parsed
            WHERE isfinite(pain) AND pain BETWEEN $pain_min AND $pain_max
                AND pain_time BETWEEN arrival_ts AND departure_ts
                AND pain_time <= arrival_ts + $lookback * INTERVAL '1 minute'
            GROUP BY visit_id, pain_time, arrival_ts
            HAVING count(DISTINCT pain) = 1
        )
        SELECT visit_id, pain AS triage_pain_score, epoch(pain_time - arrival_ts) / 60 AS triage_pain_minutes
        FROM first_clear
        QUALIFY row_number() OVER (PARTITION BY visit_id ORDER BY pain_time) = 1
        """,
        {
            "pain_measures": list(cfg.TRIAGE_PAIN_MEASURES),
            "pain_min": cfg.TRIAGE_PAIN_MIN,
            "pain_max": cfg.TRIAGE_PAIN_MAX,
            "lookback": cfg.TRIAGE_PAIN_LOOKUP_MINUTES,
        },
    )
    con.execute(
        """
        CREATE OR REPLACE TEMP TABLE p1_mimic_triage_values AS
        WITH parsed AS (
            SELECT i.stay_id, try_cast(t.acuity AS DOUBLE) AS acuity,
                   try_cast(t.pain AS DOUBLE) AS pain, try_cast(t.heartrate AS DOUBLE) AS heart_rate
            FROM p1_mimic_index i LEFT JOIN mimic_triage t USING (stay_id)
        )
        SELECT stay_id AS visit_id,
               CASE WHEN count(DISTINCT acuity) FILTER (WHERE isfinite(acuity)) = 1
                    THEN min(acuity) FILTER (WHERE isfinite(acuity)) END AS triage_level,
               CASE WHEN count(DISTINCT pain) FILTER (WHERE isfinite(pain) AND pain BETWEEN $pain_min AND $pain_max) = 1
                    THEN min(pain) FILTER (WHERE isfinite(pain) AND pain BETWEEN $pain_min AND $pain_max) END AS triage_pain_score,
               CASE WHEN count(DISTINCT heart_rate) FILTER (WHERE isfinite(heart_rate)) = 1
                    THEN min(heart_rate) FILTER (WHERE isfinite(heart_rate)) END AS triage_heart_rate
        FROM parsed GROUP BY stay_id
        """,
        {"pain_min": cfg.TRIAGE_PAIN_MIN, "pain_max": cfg.TRIAGE_PAIN_MAX},
    )


def create_labs(con):
    con.execute(
        """
        CREATE OR REPLACE TEMP TABLE p1_mcmed_lab_values AS
        WITH parsed AS (
            SELECT i.CSN AS visit_id, d.lab_measure, d.lab_priority,
                   try_cast(l.Order_time AS TIMESTAMP) AS order_ts,
                   try_cast(l.Result_time AS TIMESTAMP) AS result_ts,
                   try_cast(l.Component_value AS DOUBLE) AS value_numeric,
                   try_cast(l.Component_result AS DOUBLE) AS result_numeric,
                   i.arrival_ts, i.departure_ts
            FROM mcmed_labs l JOIN p1_mcmed_index i USING (CSN)
            JOIN p1_dict_mcmed_labs d USING (Component_name)
            WHERE d.lab_measure IS NOT NULL
        ), numeric_results AS (
            SELECT visit_id, lab_measure, lab_priority, order_ts, result_ts,
                   coalesce(CASE WHEN isfinite(value_numeric) THEN value_numeric END,
                            CASE WHEN isfinite(result_numeric) THEN result_numeric END) AS lab_value
            FROM parsed
            WHERE order_ts BETWEEN arrival_ts AND departure_ts AND result_ts BETWEEN order_ts AND departure_ts
        ), without_conflicts AS (
            SELECT visit_id, lab_measure, lab_priority, order_ts, result_ts, min(lab_value) AS lab_value
            FROM numeric_results WHERE lab_value IS NOT NULL
            GROUP BY visit_id, lab_measure, lab_priority, order_ts, result_ts
            HAVING count(DISTINCT lab_value) = 1
        ), first_result AS (
            SELECT visit_id, lab_measure, lab_priority, lab_value FROM without_conflicts
            QUALIFY row_number() OVER (PARTITION BY visit_id, lab_measure ORDER BY lab_priority, order_ts, result_ts) = 1
        )
        SELECT visit_id,
               max(lab_value) FILTER (WHERE lab_measure = 'hemoglobin') AS hemoglobin,
               max(lab_value) FILTER (WHERE lab_measure = 'platelet_count') AS platelet_count,
               max(lab_value) FILTER (WHERE lab_measure = 'inr') AS inr,
               max(lab_priority = $fallback) FILTER (WHERE lab_measure = 'hemoglobin') AS hemoglobin_from_poc,
               max(lab_priority = $fallback) FILTER (WHERE lab_measure = 'inr') AS inr_from_poc
        FROM first_result GROUP BY visit_id
        """,
        {"fallback": max(priority for _, priority in cfg.G0_LAB_LABELS.values())},
    )


def create_mcmed_order_events(con, g1=False):
    reversal = "coalesce(r.reversal_confirmed, false)" if g1 else "false"
    reversal_join = """
        LEFT JOIN p2_dict_mcmed_reversal r ON o.Order_type IS NOT DISTINCT FROM r.Order_type
                                         AND o.Procedure_name IS NOT DISTINCT FROM r.Procedure_name
    """ if g1 else ""
    con.execute(
        f"""
        CREATE OR REPLACE TEMP TABLE p1_mcmed_order_events AS
        SELECT i.CSN AS visit_id, i.arrival_ts, i.departure_ts,
               upper(trim(o.Order_type)) = 'MEDICATIONS' AS medication_order,
               try_cast(o.First_admin_time AS TIMESTAMP) AS admin_ts,
               nullif(trim(o.First_admin_time), '') IS NULL AS admin_blank,
               try_cast(o.Order_time AS TIMESTAMP) AS order_ts,
               d.aht_acute, d.aht_any, d.aht_alpha1_sensitivity, d.anxiolytic, d.opioid,
               d.hemostasis_level, d.txa_route, d.ent_referral, d.reversal_candidate,
               {reversal} AS reversal_confirmed
        FROM mcmed_orders o JOIN p1_mcmed_index i USING (CSN)
        JOIN p1_dict_mcmed_orders d ON o.Order_type IS NOT DISTINCT FROM d.Order_type
                                   AND o.Procedure_name IS NOT DISTINCT FROM d.Procedure_name
        {reversal_join}
        WHERE d.aht_acute OR d.aht_any OR d.aht_alpha1_sensitivity OR d.anxiolytic OR d.opioid
           OR d.hemostasis_level IS NOT NULL OR d.ent_referral OR d.reversal_candidate
        """
    )


def create_order_events(con):
    create_mcmed_order_events(con)
    con.execute(
        """
        CREATE OR REPLACE TEMP TABLE p1_mimic_order_events AS
        SELECT i.stay_id AS visit_id, i.arrival_ts, i.departure_ts,
               try_cast(p.charttime AS TIMESTAMP) AS admin_ts,
               d.aht_acute, d.aht_any, d.aht_alpha1_sensitivity, d.anxiolytic, d.opioid
        FROM mimic_pyxis p JOIN p1_mimic_index i USING (stay_id)
        JOIN p1_dict_mimic_pyxis d ON p.name = d.name
        WHERE d.aht_acute OR d.aht_any OR d.aht_alpha1_sensitivity OR d.anxiolytic OR d.opioid
        """
    )
    for site, med_filter in (("mcmed", "AND medication_order"), ("mimic", "")):
        con.execute(
            f"""
            CREATE OR REPLACE TEMP TABLE p1_{site}_ed_meds AS
            SELECT visit_id,
                   count(*) FILTER (WHERE aht_acute) > 0 AS ed_acute_aht,
                   count(*) FILTER (WHERE aht_any) > 0 AS ed_any_aht,
                   count(*) FILTER (WHERE aht_alpha1_sensitivity) > 0 AS ed_alpha1_sensitivity,
                   count(*) FILTER (WHERE anxiolytic) > 0 AS ed_anxiolytic,
                   count(*) FILTER (WHERE opioid) > 0 AS ed_opioid,
                   min(epoch(admin_ts - arrival_ts) / 60) FILTER (WHERE aht_acute) AS ed_first_acute_aht_minutes,
                   min(epoch(admin_ts - arrival_ts) / 60) FILTER (WHERE aht_any) AS ed_first_any_aht_minutes
            FROM p1_{site}_order_events
            WHERE admin_ts BETWEEN arrival_ts AND departure_ts {med_filter}
            GROUP BY visit_id
            """
        )


def create_mcmed_escalation(con, g1=False):
    riskset = "in_main_riskset" if g1 else "departure_ts > landmark_ts"
    g1_fields = """,
               reversal_confirmed AS escalation_reversal_confirmed,
               in_broad_riskset AS landmark_eligible_any_admission_diagnostic,
               coalesce(first_broad_escalation <= landmark_ts, false) AS hemostasis_level3_by_60min_any_admission_diagnostic,
               CASE WHEN in_broad_riskset THEN coalesce(last_broad_escalation > landmark_ts, false) END
                   AS hemostasis_level3_after_60min_any_admission_diagnostic
    """ if g1 else ""
    con.execute(
        f"""
        CREATE OR REPLACE TEMP TABLE p1_mcmed_hemostasis AS
        WITH in_ed_orders AS (
            SELECT visit_id,
                   count(*) FILTER (WHERE hemostasis_level = 2) > 0 AS hemostasis_level2,
                   count(*) FILTER (WHERE hemostasis_level = 2 AND txa_route = 'topical_or_nebulized') > 0 AS level2_txa_local,
                   count(*) FILTER (WHERE hemostasis_level = 2 AND txa_route = 'intravenous_or_oral') > 0 AS level2_txa_systemic,
                   count(*) FILTER (WHERE hemostasis_level = 2 AND txa_route = 'unknown') > 0 AS level2_txa_route_unknown,
                   count(*) FILTER (WHERE ent_referral) > 0 AS ent_referral,
                   count(*) FILTER (WHERE reversal_candidate) > 0 AS escalation_reversal_candidate,
                   count(*) FILTER (WHERE reversal_confirmed) > 0 AS reversal_confirmed,
                   min(order_ts) FILTER (WHERE hemostasis_level = 3 OR ($g1 AND reversal_confirmed)) AS first_order_escalation,
                   max(order_ts) FILTER (WHERE hemostasis_level = 3 OR ($g1 AND reversal_confirmed)) AS last_order_escalation
            FROM p1_mcmed_order_events
            WHERE order_ts BETWEEN arrival_ts AND departure_ts
            GROUP BY visit_id
        ), admissions AS (
            SELECT CSN, dx_match, arrival_ts,
                   arrival_ts + $landmark * INTERVAL '1 minute' AS landmark_ts, departure_ts,
                   CASE WHEN list_contains($admission_dispositions, trim(ED_dispo))
                        THEN try_cast(Dispo_time AS TIMESTAMP) END AS admission_ts
            FROM p1_mcmed_index
        ), valid_escalation_times AS (
            SELECT a.CSN AS visit_id, a.arrival_ts, a.departure_ts, a.landmark_ts,
                   coalesce(o.hemostasis_level2, false) AS hemostasis_level2,
                   coalesce(o.level2_txa_local, false) AS level2_txa_local,
                   coalesce(o.level2_txa_systemic, false) AS level2_txa_systemic,
                   coalesce(o.level2_txa_route_unknown, false) AS level2_txa_route_unknown,
                   coalesce(o.ent_referral, false) AS ent_referral,
                   coalesce(o.escalation_reversal_candidate, false) AS escalation_reversal_candidate,
                   coalesce(o.reversal_confirmed, false) AS reversal_confirmed,
                   CASE WHEN a.admission_ts BETWEEN a.arrival_ts AND a.departure_ts THEN a.admission_ts END AS any_admission_ts,
                   CASE WHEN a.dx_match AND a.admission_ts BETWEEN a.arrival_ts AND a.departure_ts
                        THEN a.admission_ts END AS epistaxis_admission_ts,
                   o.first_order_escalation, o.last_order_escalation
            FROM admissions a LEFT JOIN in_ed_orders o ON a.CSN = o.visit_id
        ), summary AS (
            SELECT *, least(first_order_escalation, epistaxis_admission_ts) AS first_escalation,
                   greatest(last_order_escalation, epistaxis_admission_ts) AS last_escalation,
                   least(first_order_escalation, any_admission_ts) AS first_broad_escalation,
                   greatest(last_order_escalation, any_admission_ts) AS last_broad_escalation
            FROM valid_escalation_times
        ), risksets AS (
            SELECT *, departure_ts > landmark_ts AND NOT coalesce(first_escalation <= landmark_ts, false) AS in_main_riskset,
                   departure_ts > landmark_ts AND NOT coalesce(first_broad_escalation <= landmark_ts, false) AS in_broad_riskset
            FROM summary
        )
        SELECT visit_id, hemostasis_level2, level2_txa_local, level2_txa_systemic, level2_txa_route_unknown, ent_referral,
               first_escalation IS NOT NULL AS hemostasis_level3,
               coalesce(first_escalation <= landmark_ts, false) AS hemostasis_level3_by_60min,
               CASE WHEN {riskset} THEN coalesce(last_escalation > landmark_ts, false) END AS hemostasis_level3_after_60min,
               {riskset} AS landmark_eligible,
               first_order_escalation IS NOT NULL OR any_admission_ts IS NOT NULL AS hemostasis_level3_any_admission_diagnostic,
               epistaxis_admission_ts IS NOT NULL AS escalation_admission_epistaxis_primary,
               any_admission_ts IS NOT NULL AS escalation_admission_any,
               escalation_reversal_candidate,
               epoch(first_escalation - arrival_ts) / 60 AS first_escalation_minutes
               {g1_fields}
        FROM risksets
        """,
        {"landmark": cfg.ESC_LANDMARK_MINUTES, "admission_dispositions": list(cfg.ESC_ADMISSION_DISPOSITIONS), "g1": g1},
    )


def assemble_covariates(con):
    groups = [(label, level) for label, levels in cfg.TRIAGE_ACUITY_GROUPS for level in levels]
    con.execute(
        "CREATE OR REPLACE TEMP TABLE p1_triage_groups AS SELECT * FROM (VALUES "
        + ", ".join(["(?, ?)"] * len(groups))
        + ") AS v(acuity_group, acuity_level)",
        [value for pair in groups for value in pair],
    )
    con.execute(
        """
        CREATE OR REPLACE TEMP TABLE p1_mcmed_covariates AS
        SELECT i.CSN AS visit_id,
               CASE upper(trim(i.Gender)) WHEN 'F' THEN 'Female' WHEN 'FEMALE' THEN 'Female'
                                          WHEN 'M' THEN 'Male' WHEN 'MALE' THEN 'Male' END AS sex,
               nullif(trim(i.Race), '') AS race, nullif(trim(i.Ethnicity), '') AS ethnicity,
               nullif(trim(i.Means_of_arrival), '') AS arrival_mode,
               g.acuity_level AS triage_level, g.acuity_group AS triage_acuity_group,
               t.triage_pain_score, t.triage_pain_minutes,
               CASE WHEN isfinite(try_cast(i.Triage_HR AS DOUBLE)) THEN try_cast(i.Triage_HR AS DOUBLE) END AS triage_heart_rate,
               h.has_medication_list, d.prior_htn_diagnosis, d.prior_htn_same_day,
               h.prior_htn_medication, h.prior_htn_home_sensitivity_only,
               d.prior_htn_diagnosis OR h.prior_htn_medication AS prior_htn,
               d.prior_htn_diagnosis OR d.prior_htn_same_day OR h.prior_htn_medication AS prior_htn_same_day_sensitivity,
               d.ckd, d.liver_disease, d.hht, i.is_c1 AND NOT (d.prior_htn_diagnosis OR h.prior_htn_medication) AS is_c3,
               h.anticoagulant_use, h.antiplatelet_use, h.anticoagulant_use OR h.antiplatelet_use AS any_antithrombotic_use,
               coalesce(e.ed_acute_aht, false) AS ed_acute_aht, coalesce(e.ed_any_aht, false) AS ed_any_aht,
               coalesce(e.ed_alpha1_sensitivity, false) AS ed_alpha1_sensitivity,
               coalesce(e.ed_anxiolytic, false) AS ed_anxiolytic, coalesce(e.ed_opioid, false) AS ed_opioid,
               e.ed_first_acute_aht_minutes, e.ed_first_any_aht_minutes,
               s.* EXCLUDE (visit_id), l.* EXCLUDE (visit_id)
        FROM p1_mcmed_index i
        JOIN p1_mcmed_pmh_history d ON i.CSN = d.visit_id
        JOIN p1_mcmed_home_meds h ON i.CSN = h.visit_id
        LEFT JOIN p1_mcmed_triage_pain t ON i.CSN = t.visit_id
        LEFT JOIN p1_triage_groups g
          ON coalesce(try_cast(i.Triage_acuity AS DOUBLE),
                      map_extract_value(map($acuity_codes, $acuity_levels), upper(trim(i.Triage_acuity)))) = g.acuity_level
        LEFT JOIN p1_mcmed_ed_meds e ON i.CSN = e.visit_id
        JOIN p1_mcmed_hemostasis s ON i.CSN = s.visit_id
        LEFT JOIN p1_mcmed_lab_values l ON i.CSN = l.visit_id
        """,
        {"acuity_codes": list(cfg.MCMED_TRIAGE_ACUITY_CODES),
         "acuity_levels": list(cfg.MCMED_TRIAGE_ACUITY_CODES.values())},
    )
    con.execute(
        """
        CREATE OR REPLACE TEMP TABLE p1_mimic_covariates AS
        SELECT i.stay_id AS visit_id,
               CASE upper(trim(i.gender)) WHEN 'F' THEN 'Female' WHEN 'FEMALE' THEN 'Female'
                                          WHEN 'M' THEN 'Male' WHEN 'MALE' THEN 'Male' END AS sex,
               nullif(trim(i.race), '') AS race, NULL::VARCHAR AS ethnicity,
               nullif(trim(i.arrival_transport), '') AS arrival_mode,
               g.acuity_level AS triage_level, g.acuity_group AS triage_acuity_group,
               t.triage_pain_score, CASE WHEN t.triage_pain_score IS NOT NULL THEN 0.0::DOUBLE END AS triage_pain_minutes,
               t.triage_heart_rate,
               h.has_medication_list, d.prior_htn_diagnosis, NULL::BOOLEAN AS prior_htn_same_day,
               h.prior_htn_medication, h.prior_htn_home_sensitivity_only,
               d.prior_htn_diagnosis OR h.prior_htn_medication AS prior_htn,
               NULL::BOOLEAN AS prior_htn_same_day_sensitivity,
               d.ckd, d.liver_disease, d.hht, i.is_c1 AND NOT (d.prior_htn_diagnosis OR h.prior_htn_medication) AS is_c3,
               h.anticoagulant_use, h.antiplatelet_use, h.anticoagulant_use OR h.antiplatelet_use AS any_antithrombotic_use,
               coalesce(e.ed_acute_aht, false) AS ed_acute_aht, coalesce(e.ed_any_aht, false) AS ed_any_aht,
               coalesce(e.ed_alpha1_sensitivity, false) AS ed_alpha1_sensitivity,
               coalesce(e.ed_anxiolytic, false) AS ed_anxiolytic, coalesce(e.ed_opioid, false) AS ed_opioid,
               e.ed_first_acute_aht_minutes, e.ed_first_any_aht_minutes,
               NULL::BOOLEAN AS hemostasis_level2, NULL::BOOLEAN AS level2_txa_local, NULL::BOOLEAN AS level2_txa_systemic,
               NULL::BOOLEAN AS level2_txa_route_unknown, NULL::BOOLEAN AS ent_referral, NULL::BOOLEAN AS hemostasis_level3,
               NULL::BOOLEAN AS hemostasis_level3_by_60min, NULL::BOOLEAN AS hemostasis_level3_after_60min,
               NULL::BOOLEAN AS landmark_eligible, NULL::BOOLEAN AS hemostasis_level3_any_admission_diagnostic,
               NULL::BOOLEAN AS escalation_admission_epistaxis_primary, NULL::BOOLEAN AS escalation_admission_any,
               NULL::BOOLEAN AS escalation_reversal_candidate, NULL::DOUBLE AS first_escalation_minutes,
               NULL::DOUBLE AS hemoglobin, NULL::DOUBLE AS platelet_count, NULL::DOUBLE AS inr,
               NULL::BOOLEAN AS hemoglobin_from_poc, NULL::BOOLEAN AS inr_from_poc
        FROM p1_mimic_index i
        JOIN p1_mimic_diagnosis_history d ON i.stay_id = d.visit_id
        JOIN p1_mimic_home_meds h ON i.stay_id = h.visit_id
        LEFT JOIN p1_mimic_triage_values t ON i.stay_id = t.visit_id
        LEFT JOIN p1_triage_groups g ON t.triage_level = g.acuity_level
        LEFT JOIN p1_mimic_ed_meds e ON i.stay_id = e.visit_id
        """
    )


def write_medication_audit(con):
    header = (
        "Medication category", "Matched order rows", "Matched index visits", "Blank First_admin_time rows",
        "Index visits with a blank First_admin_time", "Unparseable nonblank First_admin_time rows",
        "Rows with a parsed administration time outside ED", "Rows with a parsed administration time in ED",
        "Index visits with a parsed administration time in ED",
    )
    rows = con.execute(
        """
        WITH categories AS (
            SELECT * FROM (VALUES ('Acute BP-lowering'), ('Any BP-lowering'),
                ('Alpha-1 sensitivity only'), ('Anxiolytic'), ('Opioid')) AS c(category)
        ), events AS (
            SELECT o.visit_id, o.arrival_ts, o.departure_ts, o.admin_ts, o.admin_blank, c.category
            FROM p1_mcmed_order_events o
            CROSS JOIN LATERAL (VALUES ('Acute BP-lowering', o.aht_acute), ('Any BP-lowering', o.aht_any),
                ('Alpha-1 sensitivity only', o.aht_alpha1_sensitivity), ('Anxiolytic', o.anxiolytic),
                ('Opioid', o.opioid)) AS c(category, flagged)
            WHERE o.medication_order AND c.flagged
        )
        SELECT c.category, count(e.visit_id), count(DISTINCT e.visit_id),
               count(*) FILTER (WHERE e.admin_blank), count(DISTINCT e.visit_id) FILTER (WHERE e.admin_blank),
               count(*) FILTER (WHERE e.admin_ts IS NULL AND NOT e.admin_blank),
               count(*) FILTER (WHERE e.admin_ts < e.arrival_ts OR e.admin_ts > e.departure_ts),
               count(*) FILTER (WHERE e.admin_ts BETWEEN e.arrival_ts AND e.departure_ts),
               count(DISTINCT e.visit_id) FILTER (WHERE e.admin_ts BETWEEN e.arrival_ts AND e.departure_ts)
        FROM categories c LEFT JOIN events e USING (category)
        GROUP BY c.category ORDER BY c.category
        """
    ).fetchall()
    output_path = Path("outputs/tables/p1_medication_audit.csv")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8-sig", newline="") as output:
        writer = csv.writer(output)
        writer.writerow(header)
        writer.writerows(rows)


def run(con):
    create_mcmed_history(con)
    create_mimic_history(con)
    create_triage(con)
    create_labs(con)
    create_order_events(con)
    create_mcmed_escalation(con)
    assemble_covariates(con)
    for site in ("mcmed", "mimic"):
        counts = con.execute(
            f"SELECT count(*), count(DISTINCT visit_id), (SELECT count(*) FROM p1_{site}_index) "
            f"FROM p1_{site}_covariates"
        ).fetchone()
        if not counts[0] == counts[1] == counts[2]:
            raise RuntimeError("Covariate rows must match the full index visit set exactly once.")
    write_medication_audit(con)
