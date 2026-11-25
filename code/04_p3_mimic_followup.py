"""Ascertain outpatient blood-pressure follow-up for MIMIC-IV-ED Aim 2."""

import csv
from pathlib import Path

import duckdb

import config as cfg


TABLES = Path(__file__).resolve().parents[1] / "outputs" / "tables"


def prepare_index(con):
    con.execute("""
        CREATE OR REPLACE TEMP TABLE p3_mimic_omr_index AS
        SELECT 'primary' AS index_variant, visit_id, patient_id, departure_ts::DATE AS discharge_day
        FROM p3_mimic_index
        UNION ALL
        SELECT 'no_prior_ed', visit_id, patient_id, departure_ts::DATE
        FROM p3_mimic_index_no_prior_ed
    """)
    invalid = con.execute("""
        SELECT count(*) > 0
        FROM p3_mimic_omr_index
        WHERE visit_id IS NULL OR patient_id IS NULL OR coalesce(NOT isfinite(discharge_day), true)
    """).fetchone()[0]
    duplicates = con.execute("""
        SELECT count(*) > 0
        FROM (
            SELECT index_variant, patient_id
            FROM p3_mimic_omr_index
            GROUP BY index_variant, patient_id HAVING count(*) > 1
            UNION ALL
            SELECT index_variant, visit_id
            FROM p3_mimic_omr_index
            GROUP BY index_variant, visit_id HAVING count(*) > 1
        )
    """).fetchone()[0]
    if invalid or duplicates:
        raise ValueError("MIMIC follow-up index requires one valid visit per patient in each variant.")


def screen_omr(con):
    con.execute("""
        CREATE OR REPLACE TEMP TABLE p3_mimic_omr_candidates AS
        WITH extracted AS (
            SELECT o.subject_id AS patient_id, o.result_name,
                   try_cast(o.chartdate AS DATE) AS chart_day,
                   coalesce(regexp_full_match(o.result_value, $pattern), false) AS matches_bp_pair,
                   try_cast(regexp_extract(o.result_value, $pattern, 1) AS DOUBLE) AS sbp,
                   try_cast(regexp_extract(o.result_value, $pattern, 2) AS DOUBLE) AS dbp
            FROM mimic_hosp_omr AS o
            SEMI JOIN p3_mimic_omr_index AS i ON i.patient_id = o.subject_id
            WHERE o.result_name IN (SELECT unnest($labels))
        ), linked AS (
            SELECT i.index_variant, i.visit_id, i.patient_id, o.result_name,
                   CASE WHEN isfinite(o.chart_day) THEN o.chart_day END AS chart_day,
                   CASE WHEN isfinite(o.chart_day) THEN o.chart_day - i.discharge_day END
                       AS days_after_discharge,
                   o.sbp, o.dbp,
                   coalesce(o.matches_bp_pair AND isfinite(o.sbp) AND isfinite(o.dbp), false) AS parsed,
                   o.sbp NOT BETWEEN $sbp_min AND $sbp_max AS sbp_outside_bounds,
                   o.dbp NOT BETWEEN $dbp_min AND $dbp_max AS dbp_outside_bounds,
                   o.sbp - o.dbp < $pulse_min AS narrow_pulse_pressure,
                   coalesce(o.chart_day BETWEEN i.discharge_day + $day_min::INTEGER
                                            AND i.discharge_day + $day_max::INTEGER, false) AS in_window
            FROM p3_mimic_omr_index AS i
            JOIN extracted AS o USING (patient_id)
        )
        SELECT *, parsed AND NOT (sbp_outside_bounds OR dbp_outside_bounds OR narrow_pulse_pressure)
            AS plausible
        FROM linked
    """, {
        "pattern": cfg.P3_OMR_BP_PATTERN,
        "labels": list(cfg.P3_OMR_RESULT_NAMES),
        "day_min": cfg.P3_MIMIC_FOLLOWUP_FIRST_DAY,
        "day_max": cfg.P3_FOLLOWUP_LAST_DAY,
        "sbp_min": cfg.BP_SBP_MIN,
        "sbp_max": cfg.BP_SBP_MAX,
        "dbp_min": cfg.BP_DBP_MIN,
        "dbp_max": cfg.BP_DBP_MAX,
        "pulse_min": cfg.BP_MIN_PULSE_PRESSURE,
    })
    con.execute("""
        CREATE OR REPLACE TEMP TABLE p3_mimic_omr_acute_bounds AS
        WITH dates AS (
            SELECT 'Inpatient' AS setting, a.subject_id AS patient_id,
                   try_cast(a.admittime AS TIMESTAMP)::DATE AS start_day,
                   try_cast(a.dischtime AS TIMESTAMP)::DATE AS end_day
            FROM mimic_hosp_admissions AS a
            SEMI JOIN p3_mimic_omr_index AS i ON i.patient_id = a.subject_id
            UNION ALL
            SELECT 'ED', e.subject_id,
                   try_cast(e.intime AS TIMESTAMP)::DATE,
                   try_cast(e.outtime AS TIMESTAMP)::DATE
            FROM mimic_edstays AS e
            SEMI JOIN p3_mimic_omr_index AS i ON i.patient_id = e.subject_id
        )
        SELECT setting, patient_id,
               CASE WHEN isfinite(start_day) THEN start_day END AS start_day,
               CASE WHEN isfinite(end_day) THEN end_day END AS end_day
        FROM dates
    """)
    con.execute("""
        CREATE OR REPLACE TEMP TABLE p3_mimic_omr_acute_flags AS
        SELECT c.*,
            EXISTS (
                SELECT 1 FROM p3_mimic_omr_acute_bounds AS a
                WHERE a.patient_id = c.patient_id AND a.setting = 'Inpatient'
                  AND (c.chart_day BETWEEN a.start_day AND a.end_day
                       OR c.chart_day = a.start_day OR c.chart_day = a.end_day)
            ) AS on_inpatient_date,
            EXISTS (
                SELECT 1 FROM p3_mimic_omr_acute_bounds AS a
                WHERE a.patient_id = c.patient_id AND a.setting = 'ED'
                  AND (c.chart_day BETWEEN a.start_day AND a.end_day
                       OR c.chart_day = a.start_day OR c.chart_day = a.end_day)
            ) AS on_ed_date
        FROM p3_mimic_omr_candidates AS c
        WHERE c.in_window AND c.plausible
    """)


def summarize_followup(con):
    # With missing or reversed acute bounds, the equality checks above exclude each known
    # endpoint only. No dates are inferred beyond the available endpoints.
    con.execute("""
        CREATE OR REPLACE TEMP TABLE p3_mimic_followup AS
        WITH retained AS (
            SELECT index_variant, visit_id, patient_id, count(*) AS omr_readings,
                   count(DISTINCT chart_day) AS omr_reading_days,
                   avg(sbp) AS mean_sbp, avg(dbp) AS mean_dbp
            FROM p3_mimic_omr_acute_flags
            WHERE NOT on_inpatient_date AND NOT on_ed_date
            GROUP BY index_variant, visit_id, patient_id
        ), patient_means AS (
            SELECT i.index_variant, i.visit_id, i.patient_id,
                   coalesce(r.omr_readings, 0) AS omr_readings,
                   coalesce(r.omr_reading_days, 0) AS omr_reading_days,
                   r.mean_sbp, r.mean_dbp
            FROM p3_mimic_omr_index AS i
            LEFT JOIN retained AS r USING (index_variant, visit_id, patient_id)
        )
        SELECT *,
               omr_readings >= $min_main AS has_followup,
               CASE WHEN omr_readings >= $min_main
                    THEN mean_sbp >= $sbp_main OR mean_dbp >= $dbp_main END AS htn_evidence,
               omr_readings >= $min_strict AS has_followup_strict,
               CASE WHEN omr_readings >= $min_strict
                    THEN mean_sbp >= $sbp_strict OR mean_dbp >= $dbp_strict END AS htn_evidence_strict
        FROM patient_means
    """, {
        "min_main": cfg.P3_OMR_MIN_READINGS,
        "min_strict": cfg.P3_OMR_STRICT_MIN_READINGS,
        "sbp_main": cfg.BP_ELEVATED_SBP,
        "dbp_main": cfg.BP_ELEVATED_DBP,
        "sbp_strict": cfg.P3_OMR_STRICT_SBP,
        "dbp_strict": cfg.P3_OMR_STRICT_DBP,
    })


def write_aggregate(con, tables, name, query, parameters=None):
    result = con.execute(query, parameters)
    fields = [column[0] for column in result.description]
    records = result.fetchall()
    with (tables / f"p3_mimic_omr_{name}.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(fields)
        writer.writerows(records)


def audit_extract(con, tables):
    write_aggregate(con, tables, "flow", """
        WITH variants AS (
            SELECT * FROM (VALUES ('primary'), ('no_prior_ed')) AS v(index_variant)
        ), candidates AS (
            SELECT v.index_variant,
                count(c.result_name) AS exact_label_readings_all_dates,
                count(DISTINCT c.patient_id) AS patients_with_exact_label,
                count(*) FILTER (WHERE c.result_name IS NOT NULL AND c.chart_day IS NULL)
                    AS invalid_chart_date,
                count(*) FILTER (WHERE c.result_name IS NOT NULL AND NOT c.parsed)
                    AS invalid_bp_parse_all_dates,
                count(*) FILTER (WHERE c.in_window) AS exact_label_readings_in_window,
                count(DISTINCT c.patient_id) FILTER (WHERE c.in_window) AS patients_in_window,
                count(*) FILTER (WHERE c.in_window AND NOT c.parsed) AS invalid_bp_parse_in_window,
                count(*) FILTER (WHERE c.in_window AND c.parsed) AS parsed_readings_in_window,
                count(*) FILTER (WHERE c.in_window AND c.parsed AND c.sbp_outside_bounds)
                    AS failed_sbp_bounds,
                count(*) FILTER (WHERE c.in_window AND c.parsed AND c.dbp_outside_bounds)
                    AS failed_dbp_bounds,
                count(*) FILTER (WHERE c.in_window AND c.parsed AND c.narrow_pulse_pressure)
                    AS failed_pulse_pressure,
                count(*) FILTER (WHERE c.in_window AND c.parsed AND NOT c.plausible)
                    AS failed_any_bp_plausibility,
                count(*) FILTER (WHERE c.in_window AND c.plausible) AS plausible_readings_in_window
            FROM variants AS v
            LEFT JOIN p3_mimic_omr_candidates AS c USING (index_variant)
            GROUP BY v.index_variant
        ), exclusions AS (
            SELECT v.index_variant,
                count(*) FILTER (WHERE a.on_inpatient_date AND NOT a.on_ed_date)
                    AS inpatient_only_readings,
                count(*) FILTER (WHERE a.on_ed_date AND NOT a.on_inpatient_date)
                    AS ed_only_readings,
                count(*) FILTER (WHERE a.on_inpatient_date AND a.on_ed_date)
                    AS both_inpatient_and_ed_readings,
                count(*) FILTER (WHERE NOT a.on_inpatient_date AND NOT a.on_ed_date)
                    AS qualifying_outpatient_readings,
                min(a.days_after_discharge)
                    FILTER (WHERE NOT a.on_inpatient_date AND NOT a.on_ed_date) AS days_after_discharge_min,
                quantile_cont(a.days_after_discharge, 0.25)
                    FILTER (WHERE NOT a.on_inpatient_date AND NOT a.on_ed_date) AS days_after_discharge_q25,
                median(a.days_after_discharge)
                    FILTER (WHERE NOT a.on_inpatient_date AND NOT a.on_ed_date) AS days_after_discharge_median,
                quantile_cont(a.days_after_discharge, 0.75)
                    FILTER (WHERE NOT a.on_inpatient_date AND NOT a.on_ed_date) AS days_after_discharge_q75,
                max(a.days_after_discharge)
                    FILTER (WHERE NOT a.on_inpatient_date AND NOT a.on_ed_date) AS days_after_discharge_max,
                count(*) FILTER (WHERE NOT a.on_inpatient_date AND NOT a.on_ed_date
                    AND a.days_after_discharge = $first_day) AS qualifying_readings_on_first_day,
                count(*) FILTER (WHERE NOT a.on_inpatient_date AND NOT a.on_ed_date
                    AND a.days_after_discharge = $last_day) AS qualifying_readings_on_last_day
            FROM variants AS v
            LEFT JOIN p3_mimic_omr_acute_flags AS a USING (index_variant)
            GROUP BY v.index_variant
        ), followup AS (
            SELECT v.index_variant, count(f.patient_id) AS index_patients,
                count(*) FILTER (WHERE f.has_followup) AS patients_with_followup,
                count(*) FILTER (WHERE f.htn_evidence) AS patients_with_htn_evidence,
                count(*) FILTER (WHERE f.has_followup_strict) AS patients_with_strict_followup,
                count(*) FILTER (WHERE f.has_followup_strict AND f.omr_reading_days = 1)
                    AS patients_with_strict_followup_on_one_day,
                count(*) FILTER (WHERE f.htn_evidence_strict) AS patients_with_strict_htn_evidence
            FROM variants AS v
            LEFT JOIN p3_mimic_followup AS f USING (index_variant)
            GROUP BY v.index_variant
        )
        SELECT c.index_variant, f.index_patients, c.* EXCLUDE (index_variant),
               100.0 * c.invalid_bp_parse_all_dates / nullif(c.exact_label_readings_all_dates, 0)
                   AS invalid_bp_parse_all_dates_pct,
               100.0 * c.invalid_bp_parse_in_window / nullif(c.exact_label_readings_in_window, 0)
                   AS invalid_bp_parse_in_window_pct,
               e.* EXCLUDE (index_variant), f.* EXCLUDE (index_variant, index_patients)
        FROM candidates AS c
        JOIN exclusions AS e USING (index_variant)
        JOIN followup AS f USING (index_variant)
        ORDER BY CASE WHEN c.index_variant = 'primary' THEN 0 ELSE 1 END
    """, {"first_day": cfg.P3_MIMIC_FOLLOWUP_FIRST_DAY, "last_day": cfg.P3_FOLLOWUP_LAST_DAY})
    write_aggregate(con, tables, "labels", """
        SELECT 'full_omr' AS population, result_name, count(*) AS exact_label_readings,
               count(DISTINCT subject_id) AS patients, NULL::BIGINT AS in_window_readings,
               NULL::BIGINT AS in_window_parse_failures
        FROM mimic_hosp_omr
        WHERE result_name IN (SELECT unnest($labels))
        GROUP BY result_name
        UNION ALL
        SELECT index_variant, result_name, count(*), count(DISTINCT patient_id),
               count(*) FILTER (WHERE in_window), count(*) FILTER (WHERE in_window AND NOT parsed)
        FROM p3_mimic_omr_candidates
        GROUP BY index_variant, result_name
        ORDER BY population, result_name
    """, {"labels": list(cfg.P3_OMR_RESULT_NAMES)})
    write_aggregate(con, tables, "acute_intervals", """
        SELECT i.index_variant, b.setting, count(*) AS all_history_encounters,
               count(*) FILTER (WHERE b.start_day IS NULL AND b.end_day IS NULL)
                   AS no_usable_day,
               count(*) FILTER (WHERE (b.start_day IS NULL) <> (b.end_day IS NULL))
                   AS one_usable_day,
               count(*) FILTER (WHERE b.start_day > b.end_day) AS reversed_date_bounds
        FROM p3_mimic_omr_index AS i
        JOIN p3_mimic_omr_acute_bounds AS b USING (patient_id)
        GROUP BY i.index_variant, b.setting
        ORDER BY i.index_variant, b.setting
    """)
    write_aggregate(con, tables, "distributions", """
        WITH groups AS (
            SELECT index_variant, patient_id, 'Parsed in window' AS scope, sbp, dbp
            FROM p3_mimic_omr_candidates
            WHERE in_window AND parsed
            UNION ALL
            SELECT index_variant, patient_id, 'Qualifying outpatient readings', sbp, dbp
            FROM p3_mimic_omr_acute_flags
            WHERE NOT on_inpatient_date AND NOT on_ed_date
            UNION ALL
            SELECT index_variant, patient_id, 'Per-patient outpatient means', mean_sbp, mean_dbp
            FROM p3_mimic_followup
            WHERE has_followup
        ), measures AS (
            SELECT index_variant, patient_id, scope, 'SBP' AS measure, sbp AS mmhg FROM groups
            UNION ALL
            SELECT index_variant, patient_id, scope, 'DBP', dbp FROM groups
        )
        SELECT index_variant, scope, measure, count(*) AS observations,
               count(DISTINCT patient_id) AS patients, quantile_cont(mmhg, 0.05) AS p05,
               quantile_cont(mmhg, 0.25) AS q25, median(mmhg) AS median,
               quantile_cont(mmhg, 0.75) AS q75, quantile_cont(mmhg, 0.95) AS p95
        FROM measures
        GROUP BY index_variant, scope, measure
        ORDER BY index_variant, scope, measure
    """)


def run(con, tables=TABLES):
    """Create the follow-up relation and write only the prespecified aggregate audits."""
    try:
        prepare_index(con)
        screen_omr(con)
        summarize_followup(con)
        tables.mkdir(parents=True, exist_ok=True)
        audit_extract(con, tables)
    except duckdb.Error as exc:
        raise RuntimeError(f"MIMIC OMR extraction failed ({type(exc).__name__}).") from None
