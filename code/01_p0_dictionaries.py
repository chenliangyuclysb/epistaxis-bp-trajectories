"""Generate the P0.7 dictionary-candidate frequency tables."""

import csv
from pathlib import Path

from config import (
    AHT_DRUG_PATTERN,
    AHT_IV_ROUTE_PATTERN,
    AHT_ORAL_ROUTE_PATTERN,
    HEMOSTASIS_ESCALATION_PATTERN,
    LAB_COMPONENT_PATTERN,
    MCMED_HTN_ICD10_PATTERN,
    MCMED_HTN_ICD9_PATTERN,
    MIMIC_MEDRECON_ACAP_AHT_PATTERN,
)


CANDIDATE_DIR = Path("outputs/dictionary_candidates")


def write_candidate(con, filename, query, parameters=()):
    output = CANDIDATE_DIR / filename
    pending = output.with_suffix(".csv.tmp")
    try:
        result = con.execute(query, parameters)
        columns = [column[0] for column in result.description]
        rows = result.fetchall()
        with pending.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(columns)
            writer.writerows(rows)
        pending.replace(output)
    except Exception:
        pending.unlink(missing_ok=True)
        raise RuntimeError(f"Failed to generate {filename}") from None


def prepare_orders(con):
    """Scan orders once and aggregate matching records by visit and label in memory."""
    try:
        con.execute(
            """
            CREATE OR REPLACE TEMP TABLE p0_dictionary_order_groups AS
            WITH matched AS (
                SELECT CSN, Order_type, Procedure_name, Order_time, First_admin_time,
                    Order_type = 'Medications' AND regexp_matches(Procedure_name, ?) AS is_aht,
                    regexp_matches(Procedure_name, ?) AS is_escalation
                FROM mcmed_orders
            ), candidate AS (
                SELECT CSN, Order_type, Procedure_name, is_aht, is_escalation,
                    TRY_CAST(Order_time AS TIMESTAMP) AS order_ts,
                    TRY_CAST(First_admin_time AS TIMESTAMP) AS administration_ts
                FROM matched
                WHERE is_aht OR is_escalation
            ), visit_times AS (
                SELECT CSN, TRY_CAST(Arrival_time AS TIMESTAMP) AS arrival_ts,
                    TRY_CAST(Departure_time AS TIMESTAMP) AS departure_ts
                FROM mcmed_visits
            )
            SELECT o.CSN, o.Order_type, o.Procedure_name, o.is_aht, o.is_escalation,
                count(*) AS order_count,
                bool_or(o.administration_ts BETWEEN v.arrival_ts AND v.departure_ts) AS ed_admin,
                bool_or(o.order_ts BETWEEN v.arrival_ts AND v.departure_ts) AS ed_order,
                min(epoch(o.order_ts - v.arrival_ts) / 60.0) AS first_order_minutes
            FROM candidate o
            LEFT JOIN visit_times v ON o.CSN = v.CSN
            GROUP BY o.CSN, o.Order_type, o.Procedure_name, o.is_aht, o.is_escalation
            """,
            (AHT_DRUG_PATTERN, HEMOSTASIS_ESCALATION_PATTERN),
        )
    except Exception:
        raise RuntimeError("Failed to aggregate dictionary-candidate orders") from None


def run(con):
    CANDIDATE_DIR.mkdir(parents=True, exist_ok=True)
    prepare_orders(con)
    try:
        write_candidate(
            con,
            "mcmed_aht_ed_candidates.csv",
            """
            WITH route_tags AS (
                SELECT *, regexp_matches(Procedure_name, ?) AS iv_label,
                    regexp_matches(Procedure_name, ?) AS oral_label
                FROM p0_dictionary_order_groups
                WHERE is_aht
            ), labeled AS (
                SELECT *, CASE
                    WHEN iv_label AND NOT oral_label THEN 'Intravenous'
                    WHEN oral_label AND NOT iv_label THEN 'Oral'
                    ELSE 'Unknown'
                END AS inferred_route
                FROM route_tags
            )
            SELECT Procedure_name, inferred_route AS Inferred_route,
                sum(order_count) AS Order_count,
                count(CSN) AS Visit_count,
                count(CSN) FILTER (WHERE ed_admin) AS ED_admin_visit_count
            FROM labeled
            GROUP BY Procedure_name, inferred_route
            ORDER BY Order_count DESC, Procedure_name
            """,
            (AHT_IV_ROUTE_PATTERN, AHT_ORAL_ROUTE_PATTERN),
        )
        write_candidate(
            con,
            "mcmed_hemostasis_escalation_candidates.csv",
            """
            SELECT Order_type, Procedure_name, sum(order_count) AS Order_count,
                count(CSN) AS Visit_count,
                count(CSN) FILTER (WHERE ed_order) AS ED_order_visit_count,
                median(first_order_minutes) FILTER (WHERE CSN IS NOT NULL)
                    AS Median_minutes_arrival_to_first_order
            FROM p0_dictionary_order_groups
            WHERE is_escalation
            GROUP BY Order_type, Procedure_name
            ORDER BY Order_count DESC, Order_type NULLS LAST, Procedure_name
            """,
        )
    finally:
        con.execute("DROP TABLE IF EXISTS p0_dictionary_order_groups")

    write_candidate(
        con,
        "mcmed_meds_class_candidates.csv",
        """
        SELECT Med_class, Med_subclass,
            count(*) AS Medication_record_count, count(DISTINCT MRN) AS Patient_count
        FROM mcmed_meds
        WHERE MRN IN (SELECT MRN FROM p0_mcmed_candidates)
        GROUP BY Med_class, Med_subclass
        ORDER BY Medication_record_count DESC, Med_class NULLS LAST, Med_subclass NULLS LAST
        """,
    )
    write_candidate(
        con,
        "mimic_pyxis_aht_candidates.csv",
        """
        SELECT name AS Medication_name,
            count(*) AS Dispensing_record_count, count(DISTINCT stay_id) AS Stay_count
        FROM mimic_pyxis
        WHERE regexp_matches(name, ?)
        GROUP BY name
        ORDER BY Dispensing_record_count DESC, Medication_name
        """,
        (AHT_DRUG_PATTERN,),
    )
    write_candidate(
        con,
        "mimic_medrecon_class_candidates.csv",
        """
        SELECT etcdescription AS Medication_class,
            count(*) AS Reconciliation_record_count,
            count(DISTINCT stay_id) AS Stay_count,
            count(DISTINCT subject_id) AS Patient_count
        FROM mimic_medrecon
        WHERE regexp_matches(etcdescription, ?)
        GROUP BY etcdescription
        ORDER BY Reconciliation_record_count DESC, Medication_class
        """,
        (MIMIC_MEDRECON_ACAP_AHT_PATTERN,),
    )
    write_candidate(
        con,
        "mcmed_pmh_htn_timing.csv",
        """
        WITH htn AS (
            SELECT p.MRN,
                CAST(TRY_CAST(p.Noted_date AS TIMESTAMP) AS DATE) AS noted_day,
                CAST(i.arrival_ts AS DATE) AS index_day
            FROM mcmed_pmh p
            INNER JOIN p0_mcmed_index i ON p.MRN = i.MRN
            WHERE (p.CodeType = 'Dx10' AND regexp_matches(trim(p.Code), ?))
               OR (p.CodeType = 'Dx9' AND regexp_matches(trim(p.Code), ?))
        ), timed AS (
            SELECT MRN, CASE
                WHEN noted_day IS NULL OR index_day IS NULL THEN 'Missing'
                WHEN noted_day < index_day THEN 'Before'
                WHEN noted_day = index_day THEN 'Same_day'
                ELSE 'After'
            END AS timing
            FROM htn
        ), counts AS (
            SELECT timing, count(*) AS HTN_code_record_count,
                count(DISTINCT MRN) AS Patient_count
            FROM timed
            GROUP BY timing
        ), groups(timing, sort_order) AS (
            VALUES ('Before', 1), ('Same_day', 2), ('After', 3), ('Missing', 4)
        )
        SELECT g.timing AS Timing, coalesce(c.HTN_code_record_count, 0) AS HTN_code_record_count,
            coalesce(c.Patient_count, 0) AS Patient_count
        FROM groups g
        LEFT JOIN counts c ON g.timing = c.timing
        ORDER BY g.sort_order
        """,
        (MCMED_HTN_ICD10_PATTERN, MCMED_HTN_ICD9_PATTERN),
    )
    write_candidate(
        con,
        "mcmed_labs_candidates.csv",
        """
        SELECT l.Component_name,
            count(*) AS All_period_record_count,
            count(*) FILTER (WHERE TRY_CAST(l.Order_time AS TIMESTAMP)
                BETWEEN v.arrival_ts AND v.departure_ts) AS ED_order_record_count
        FROM mcmed_labs l
        INNER JOIN p0_mcmed_candidates v ON l.CSN = v.CSN
        WHERE regexp_matches(l.Component_name, ?)
        GROUP BY l.Component_name
        ORDER BY All_period_record_count DESC, l.Component_name
        """,
        (LAB_COMPONENT_PATTERN,),
    )
