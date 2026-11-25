"""Derive the MC-MED supporting follow-up measure in frozen SAP S6.2."""

from importlib import import_module
from pathlib import Path

import config as cfg


TABLES = Path(__file__).resolve().parents[1] / "outputs" / "tables"


def prepare_inputs(con):
    con.execute("""
        CREATE OR REPLACE TEMP TABLE p3_mcmed_followup_indices AS
        SELECT 'primary' AS index_variant, visit_id, patient_id, arrival_ts::DATE AS index_date
        FROM p3_mcmed_index
        UNION ALL
        SELECT 'same_day_htn', visit_id, patient_id, arrival_ts::DATE
        FROM p3_mcmed_index_same_day_htn
    """)
    count, visits, patients, missing = con.execute("""
        SELECT count(*), count(DISTINCT (index_variant, visit_id)),
            count(DISTINCT (index_variant, patient_id)),
            count(*) FILTER (WHERE visit_id IS NULL OR patient_id IS NULL OR index_date IS NULL)
        FROM p3_mcmed_followup_indices
    """).fetchone()
    if count != visits or count != patients or missing:
        raise ValueError("MC-MED C3 inputs must have one valid patient and visit per index variant.")

    dictionaries = import_module("02_p1_dictionaries")
    dictionaries.map_labels(
        con, "mcmed_meds", ("Med_class", "Med_subclass"), "p3_mcmed_home_med_dictionary",
        dictionaries.HOME_FIELDS, dictionaries.classify_home_acap, None,
    )


def create_counts(con):
    con.execute("""
        CREATE OR REPLACE TEMP TABLE p3_mcmed_followup_counts AS
        WITH windows AS (
            SELECT i.*, start.first_day, start.first_day = $strict_first_day AS strict_window,
                i.index_date + start.first_day AS first_date, i.index_date + $last_day::INTEGER AS last_date
            FROM p3_mcmed_followup_indices i
            CROSS JOIN (VALUES ($first_day::INTEGER), ($strict_first_day::INTEGER)) AS start(first_day)
        ), people AS (
            SELECT DISTINCT patient_id FROM p3_mcmed_followup_indices
        ), pmh AS (
            SELECT raw.MRN AS patient_id, try_cast(raw.Noted_date AS DATE) AS noted_date,
                CASE upper(trim(raw.CodeType))
                    WHEN 'DX10' THEN regexp_matches(trim(raw.Code), $htn10)
                    WHEN 'DX9' THEN regexp_matches(trim(raw.Code), $htn9)
                    ELSE false
                END AS htn_code
            FROM mcmed_pmh raw JOIN people p ON raw.MRN = p.patient_id
        ), meds AS (
            SELECT raw.MRN AS patient_id, try_cast(raw.Entry_date AS DATE) AS entry_date,
                try_cast(raw.Start_date AS DATE) AS start_date,
                nullif(trim(raw.Start_date), '') IS NULL AS start_missing,
                coalesce(d.aht_home, false) AS home_aht
            FROM mcmed_meds raw JOIN people p ON raw.MRN = p.patient_id
            LEFT JOIN p3_mcmed_home_med_dictionary d
                ON d.Med_class IS NOT DISTINCT FROM raw.Med_class
                AND d.Med_subclass IS NOT DISTINCT FROM raw.Med_subclass
        ), visits AS (
            SELECT raw.MRN AS patient_id, raw.CSN AS visit_id,
                try_cast(raw.Arrival_time AS DATE) AS arrival_date
            FROM mcmed_visits raw JOIN people p ON raw.MRN = p.patient_id
        ), pmh_counts AS (
            SELECT i.index_variant, i.visit_id, i.strict_window,
                count(*) FILTER (WHERE p.noted_date BETWEEN i.first_date AND i.last_date) AS pmh_contact_rows,
                count(*) FILTER (WHERE p.htn_code AND p.noted_date BETWEEN i.first_date AND i.last_date)
                    AS htn_code_rows
            FROM windows i LEFT JOIN pmh p USING (patient_id)
            GROUP BY i.index_variant, i.visit_id, i.strict_window
        ), med_counts AS (
            SELECT i.index_variant, i.visit_id, i.strict_window,
                count(*) FILTER (WHERE m.entry_date BETWEEN i.first_date AND i.last_date)
                    AS med_entry_contact_rows,
                count(*) FILTER (WHERE m.home_aht AND m.start_date BETWEEN i.first_date AND i.last_date)
                    AS aht_start_date_rows,
                count(*) FILTER (WHERE m.home_aht AND m.start_missing
                    AND m.entry_date BETWEEN i.first_date AND i.last_date) AS aht_entry_fallback_rows,
                count(*) FILTER (WHERE m.home_aht AND m.start_date < i.index_date
                    AND m.entry_date BETWEEN i.first_date AND i.last_date) AS old_start_late_entry_rows,
                count(*) FILTER (WHERE m.home_aht AND m.start_date = i.index_date
                    AND m.entry_date BETWEEN i.first_date AND i.last_date) AS day0_start_late_entry_rows,
                count(*) FILTER (WHERE m.home_aht AND NOT m.start_missing AND m.start_date IS NULL
                    AND m.entry_date BETWEEN i.first_date AND i.last_date) AS invalid_start_no_fallback_rows
            FROM windows i LEFT JOIN meds m USING (patient_id)
            GROUP BY i.index_variant, i.visit_id, i.strict_window
        ), visit_counts AS (
            SELECT i.index_variant, i.visit_id, i.strict_window,
                count(DISTINCT v.visit_id) FILTER (WHERE v.visit_id <> i.visit_id
                    AND v.arrival_date BETWEEN i.first_date AND i.last_date) AS subsequent_visit_rows
            FROM windows i LEFT JOIN visits v USING (patient_id)
            GROUP BY i.index_variant, i.visit_id, i.strict_window
        )
        SELECT i.index_variant, i.visit_id, i.patient_id, i.strict_window, i.first_day,
            p.pmh_contact_rows, p.htn_code_rows, m.med_entry_contact_rows, m.aht_start_date_rows,
            m.aht_entry_fallback_rows, m.old_start_late_entry_rows, m.day0_start_late_entry_rows,
            m.invalid_start_no_fallback_rows, v.subsequent_visit_rows,
            p.pmh_contact_rows > 0 OR m.med_entry_contact_rows > 0 OR v.subsequent_visit_rows > 0
                AS has_followup,
            p.htn_code_rows > 0 AS htn_pmh_candidate,
            m.aht_start_date_rows + m.aht_entry_fallback_rows > 0 AS htn_med_candidate
        FROM windows i
        JOIN pmh_counts p USING (index_variant, visit_id, strict_window)
        JOIN med_counts m USING (index_variant, visit_id, strict_window)
        JOIN visit_counts v USING (index_variant, visit_id, strict_window)
    """, {
        "first_day": cfg.P3_MCMED_FOLLOWUP_FIRST_DAY,
        "strict_first_day": cfg.P3_MCMED_FOLLOWUP_STRICT_FIRST_DAY,
        "last_day": cfg.P3_FOLLOWUP_LAST_DAY,
        "htn10": cfg.MCMED_HTN_ICD10_PATTERN,
        "htn9": cfg.MCMED_HTN_ICD9_PATTERN,
    })


def create_followup(con):
    con.execute("""
        CREATE OR REPLACE TEMP TABLE p3_mcmed_followup AS
        SELECT p.index_variant, p.visit_id, p.patient_id, p.has_followup,
            CASE WHEN p.has_followup THEN p.htn_pmh_candidate OR p.htn_med_candidate END AS htn_evidence,
            s.has_followup AS has_followup_strict,
            CASE WHEN s.has_followup THEN s.htn_pmh_candidate OR s.htn_med_candidate END AS htn_evidence_strict
        FROM p3_mcmed_followup_counts p
        JOIN p3_mcmed_followup_counts s USING (index_variant, visit_id)
        WHERE NOT p.strict_window AND s.strict_window
    """)
    count, invalid = con.execute("""
        SELECT count(*), count(*) FILTER (WHERE
            (NOT has_followup AND (has_followup_strict OR htn_evidence IS NOT NULL))
            OR (has_followup AND htn_evidence IS NULL)
            OR (has_followup_strict AND htn_evidence_strict IS NULL)
            OR (NOT has_followup_strict AND htn_evidence_strict IS NOT NULL)
            OR (htn_evidence_strict IS TRUE AND htn_evidence IS DISTINCT FROM true)
        )
        FROM p3_mcmed_followup
    """).fetchone()
    expected = con.execute("SELECT count(*) FROM p3_mcmed_followup_indices").fetchone()[0]
    if count != expected or invalid:
        raise ValueError("MC-MED follow-up denominators, nullable evidence or window nesting did not reconcile.")


def write_audit(con, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    result = con.execute("""
        SELECT index_variant, first_day, $last_day::INTEGER AS last_day,
            count(*) AS index_patients,
            count(*) FILTER (WHERE has_followup) AS followed_patients,
            count(*) FILTER (WHERE pmh_contact_rows > 0) AS pmh_contact_patients,
            count(*) FILTER (WHERE med_entry_contact_rows > 0) AS medication_entry_contact_patients,
            count(*) FILTER (WHERE subsequent_visit_rows > 0) AS subsequent_visit_contact_patients,
            count(*) FILTER (WHERE pmh_contact_rows > 0 AND med_entry_contact_rows = 0
                AND subsequent_visit_rows = 0) AS pmh_only_contact_patients,
            count(*) FILTER (WHERE pmh_contact_rows = 0 AND med_entry_contact_rows > 0
                AND subsequent_visit_rows = 0) AS medication_entry_only_contact_patients,
            count(*) FILTER (WHERE pmh_contact_rows = 0 AND med_entry_contact_rows = 0
                AND subsequent_visit_rows > 0) AS subsequent_visit_only_contact_patients,
            count(*) FILTER (WHERE has_followup AND (htn_pmh_candidate OR htn_med_candidate)) AS evidence_patients,
            count(*) FILTER (WHERE has_followup AND htn_pmh_candidate AND NOT htn_med_candidate)
                AS htn_code_only_patients,
            count(*) FILTER (WHERE has_followup AND htn_med_candidate AND NOT htn_pmh_candidate)
                AS new_medication_only_patients,
            count(*) FILTER (WHERE has_followup AND htn_pmh_candidate AND htn_med_candidate)
                AS both_evidence_sources_patients,
            count(*) FILTER (WHERE has_followup AND aht_start_date_rows > 0)
                AS new_medication_start_date_patients,
            count(*) FILTER (WHERE has_followup AND aht_entry_fallback_rows > 0)
                AS new_medication_entry_fallback_patients,
            count(*) FILTER (WHERE NOT has_followup AND htn_med_candidate)
                AS new_medication_candidate_without_contact_patients,
            sum(pmh_contact_rows) AS pmh_contact_records,
            sum(med_entry_contact_rows) AS medication_entry_contact_records,
            sum(subsequent_visit_rows) AS subsequent_visit_contact_records,
            sum(htn_code_rows) AS htn_code_records,
            sum(aht_start_date_rows) AS aht_start_date_records,
            sum(aht_entry_fallback_rows) AS aht_entry_fallback_records,
            sum(old_start_late_entry_rows) AS old_start_late_entry_records_excluded,
            sum(day0_start_late_entry_rows) AS day0_start_late_entry_records_excluded,
            sum(invalid_start_no_fallback_rows) AS invalid_start_no_fallback_records_excluded
        FROM p3_mcmed_followup_counts
        GROUP BY index_variant, first_day
        ORDER BY index_variant, first_day
    """, {"last_day": cfg.P3_FOLLOWUP_LAST_DAY})
    header = [column[0] for column in result.description]
    import_module("01_p0_audit").write_csv(output_dir / "p3_mcmed_followup_audit.csv", header, result.fetchall())


def run(con, output_dir=TABLES):
    """Create the private per-index follow-up lookups and write the aggregate audit."""
    try:
        prepare_inputs(con)
        create_counts(con)
        create_followup(con)
        if output_dir is not None:
            write_audit(con, Path(output_dir))
    except Exception:
        raise RuntimeError("Could not derive or validate the MC-MED S6.2 follow-up measure.") from None
