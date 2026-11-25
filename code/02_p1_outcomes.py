"""Ascertain the P1 outcomes independently of covariates and BP."""

import config as cfg


def add_returns(con, site, visit_id, patient_id):
    con.execute(
        f"""
        CREATE OR REPLACE TEMP TABLE p1_{site}_returns AS
        WITH next_visit AS (
            SELECT i.{visit_id} AS visit_id, i.incomplete_visit, i.departure_ts,
                min(r.arrival_ts) AS return_ts
            FROM p1_{site}_index i
            LEFT JOIN p1_{site}_adult_episodes r ON
                i.{patient_id} = r.{patient_id}
                AND i.{visit_id} <> r.{visit_id}
                AND r.arrival_ts > i.departure_ts
                AND r.arrival_ts <= i.departure_ts + $long * INTERVAL '1 hour'
            GROUP BY i.{visit_id}, i.incomplete_visit, i.departure_ts
        )
        SELECT visit_id, NOT incomplete_visit AS return_outcome_eligible,
            CASE WHEN NOT incomplete_visit
                THEN coalesce(return_ts <= departure_ts + $short * INTERVAL '1 hour', false)
            END AS return_7d_event,
            CASE WHEN NOT incomplete_visit THEN return_ts IS NOT NULL END AS return_30d_event,
            CASE WHEN NOT incomplete_visit THEN epoch(return_ts - departure_ts) / 60.0
            END AS return_minutes_within_30d
        FROM next_visit
        """,
        {"short": cfg.RETURN_7D_HOURS, "long": cfg.RETURN_30D_HOURS},
    )


def add_inpatient_ent_hemostasis(con):
    con.execute(
        """
        CREATE OR REPLACE TEMP TABLE p1_mimic_procedure_labels AS
        SELECT icd_code, icd_version, long_title
        FROM mimic_hosp_d_icd_procedures
        WHERE regexp_matches(long_title, ?)
        """,
        [cfg.INPATIENT_ENT_HEMOSTASIS_TITLE_PATTERN],
    )
    con.execute(
        """
        CREATE OR REPLACE TEMP TABLE p1_mimic_inpatient_outcome AS
        WITH admissions AS (
            SELECT i.stay_id, i.subject_id, a.hadm_id,
                try_cast(a.admittime AS DATE) AS admission_day,
                try_cast(a.dischtime AS DATE) AS discharge_day
            FROM p1_mimic_index i
            JOIN mimic_hosp_admissions a ON a.hadm_id = i.hadm_id AND a.subject_id = i.subject_id
        ), procedure_hits AS (
            SELECT a.stay_id, TRUE AS event
            FROM mimic_hosp_procedures_icd p
            JOIN admissions a ON a.hadm_id = p.hadm_id AND a.subject_id = p.subject_id
            JOIN p1_mimic_procedure_labels l USING (icd_code, icd_version)
            WHERE try_cast(p.chartdate AS DATE) BETWEEN a.admission_day AND a.discharge_day
            GROUP BY a.stay_id
        )
        SELECT i.stay_id AS visit_id,
            coalesce(a.admission_day <= a.discharge_day, false) AS inpatient_ent_outcome_eligible,
            CASE WHEN a.admission_day <= a.discharge_day THEN coalesce(h.event, false) END
                AS inpatient_ent_hemostasis_event
        FROM p1_mimic_index i
        LEFT JOIN admissions a USING (stay_id)
        LEFT JOIN procedure_hits h USING (stay_id)
        """
    )


def run(con):
    add_returns(con, "mcmed", "CSN", "MRN")
    add_returns(con, "mimic", "stay_id", "subject_id")
    add_inpatient_ent_hemostasis(con)
