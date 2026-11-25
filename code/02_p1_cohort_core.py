"""Construct the diagnosis-based P1 index cohorts and ED blood-pressure readings."""

import csv
from pathlib import Path

import config as cfg


TABLES = Path(__file__).resolve().parents[1] / "outputs" / "tables"
SITES = (
    ("mcmed", "MC-MED", "CSN", "MRN", "ED_dispo", "DIRECT_HOME_MCMED"),
    ("mimic", "MIMIC-IV-ED", "stay_id", "subject_id", "disposition", "DIRECT_HOME_MIMIC"),
)


def write_aggregate(con, name, query):
    result = con.execute(query)
    columns = [column[0] for column in result.description]
    rows = result.fetchall()
    TABLES.mkdir(parents=True, exist_ok=True)
    with (TABLES / f"{name}.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(columns)
        writer.writerows(rows)


def prepare_candidates(con):
    con.execute(
        """
        CREATE OR REPLACE TEMP TABLE p1_mcmed_candidates AS
        WITH flagged AS (
            SELECT v.*,
                try_cast(v.Arrival_time AS TIMESTAMP) AS arrival_ts,
                try_cast(v.Departure_time AS TIMESTAMP) AS departure_ts,
                try_cast(v.Age AS DOUBLE) AS age_years,
                coalesce(regexp_matches(trim(v.Dx_ICD10), $dx10), false) AS dx_match,
                coalesce(regexp_matches(v.CC, $cc), false) AS cc_match,
                coalesce(regexp_matches(trim(v.Dx_ICD10), $trauma10)
                      OR regexp_matches(trim(v.Dx_ICD9), $trauma9), false) AS trauma_match,
                coalesce(regexp_matches(v.CC, $postop_cc)
                      OR regexp_matches(v.Dx_name, $postop_dx), false) AS postoperative_match
            FROM mcmed_visits AS v
        )
        SELECT * FROM flagged WHERE dx_match OR cc_match
        """,
        {
            "dx10": cfg.MCMED_DX_PATTERN,
            "cc": cfg.EPIS_CC_PATTERN,
            "trauma10": cfg.TRAUMA_DX_PATTERN,
            "trauma9": cfg.TRAUMA_DX9_PATTERN,
            "postop_cc": cfg.POST_OP_CC_PATTERN,
            "postop_dx": cfg.POST_OP_DX_PATTERN,
        },
    )
    con.execute(
        """
        CREATE OR REPLACE TEMP TABLE p1_mimic_candidates AS
        WITH coded AS (
            SELECT stay_id,
                bool_or(coalesce(try_cast(seq_num AS INTEGER) = $primary AND (
                    (icd_version = '10' AND regexp_matches(trim(icd_code), $dx10))
                    OR (icd_version = '9' AND trim(icd_code) = $dx9)
                ), false)) AS dx_match,
                bool_or(coalesce(
                    (icd_version = '10' AND regexp_matches(trim(icd_code), $trauma10))
                    OR (icd_version = '9' AND regexp_matches(trim(icd_code), $trauma9))
                , false)) AS trauma_match,
                bool_or(coalesce(regexp_matches(icd_title, $postop_dx), false)) AS postoperative_match
            FROM mimic_diagnosis
            GROUP BY stay_id
        ), complaints AS (
            SELECT stay_id,
                bool_or(coalesce(regexp_matches(chiefcomplaint, $cc), false)) AS cc_match,
                bool_or(coalesce(regexp_matches(chiefcomplaint, $postop_cc), false))
                    AS postoperative_match
            FROM mimic_triage
            GROUP BY stay_id
        ), flagged AS (
            SELECT e.*,
                try_cast(e.intime AS TIMESTAMP) AS arrival_ts,
                try_cast(e.outtime AS TIMESTAMP) AS departure_ts,
                coalesce(d.dx_match, false) AS dx_match,
                coalesce(t.cc_match, false) AS cc_match,
                coalesce(d.trauma_match, false) AS trauma_match,
                coalesce(d.postoperative_match, false) OR coalesce(t.postoperative_match, false)
                    AS postoperative_match
            FROM mimic_edstays AS e
            LEFT JOIN coded AS d USING (stay_id)
            LEFT JOIN complaints AS t USING (stay_id)
            WHERE coalesce(d.dx_match, false) OR coalesce(t.cc_match, false)
        )
        SELECT f.*,
            try_cast(p.anchor_age AS DOUBLE) + (year(f.arrival_ts) - try_cast(p.anchor_year AS INTEGER))
                AS age_years,
            p.subject_id IS NOT NULL AS matched_hosp_patient
        FROM flagged AS f
        LEFT JOIN mimic_hosp_patients AS p USING (subject_id)
        """,
        {
            "primary": cfg.MIMIC_PRIMARY_DX_SEQ,
            "dx10": cfg.MIMIC_DX10_PATTERN,
            "dx9": cfg.MIMIC_DX9_CODE,
            "cc": cfg.EPIS_CC_PATTERN,
            "trauma10": cfg.TRAUMA_DX_PATTERN,
            "trauma9": cfg.TRAUMA_DX9_PATTERN,
            "postop_cc": cfg.POST_OP_CC_PATTERN,
            "postop_dx": cfg.POST_OP_DX_PATTERN,
        },
    )


def prepare_episode_cohorts(con):
    for short, _, visit_id, patient_id, _, _ in SITES:
        con.execute(
            f"""
            CREATE OR REPLACE TEMP TABLE p1_{short}_screened AS
            SELECT *,
                CASE
                    WHEN age_years IS NULL OR NOT isfinite(age_years) THEN 'Age unavailable'
                    WHEN age_years < $adult THEN 'Age below adult cutoff'
                    WHEN nullif(trim({visit_id}), '') IS NULL
                      OR nullif(trim({patient_id}), '') IS NULL THEN 'Identifier unavailable'
                    WHEN trauma_match THEN 'Trauma diagnosis'
                    WHEN arrival_ts IS NULL OR departure_ts IS NULL OR departure_ts <= arrival_ts
                        THEN 'Invalid ED times'
                    WHEN postoperative_match THEN 'Postoperative or postprocedural'
                END AS exclusion_reason
            FROM p1_{short}_candidates
            """,
            {"adult": cfg.MIN_ADULT_AGE},
        )
        con.execute(
            f"""
            CREATE OR REPLACE TEMP TABLE p1_{short}_adult_episodes AS
            SELECT *
            FROM p1_{short}_screened
            WHERE isfinite(age_years) AND age_years >= $adult
            """,
            {"adult": cfg.MIN_ADULT_AGE},
        )
        con.execute(
            f"""
            CREATE OR REPLACE TEMP TABLE p1_{short}_eligible_visits AS
            WITH flags AS (
                SELECT e.*,
                    p.last_departure IS NOT NULL AS prior_epistaxis_365d,
                    p.last_departure IS NULL
                      OR p.last_departure < e.arrival_ts - $return_hours * INTERVAL '1 hour'
                        AS is_new_episode
                FROM p1_{short}_adult_episodes AS e
                LEFT JOIN LATERAL (
                    SELECT max(prior.departure_ts) AS last_departure
                    FROM p1_{short}_adult_episodes AS prior
                    WHERE prior.{patient_id} = e.{patient_id}
                      AND prior.departure_ts > prior.arrival_ts
                      AND prior.departure_ts < e.arrival_ts
                      AND prior.departure_ts >= e.arrival_ts - $prior_days * INTERVAL '1 day'
                ) AS p ON true
                WHERE e.exclusion_reason IS NULL
            )
            SELECT *,
                row_number() OVER (
                    PARTITION BY {patient_id} ORDER BY arrival_ts, {visit_id}
                ) = 1 AS is_first_eligible_visit
            FROM flags
            """,
            {
                "return_hours": cfg.RETURN_7D_HOURS,
                "prior_days": cfg.PRIOR_EPISODES_LOOKBACK_DAYS,
            },
        )
        con.execute(
            f"""
            CREATE OR REPLACE TEMP TABLE p1_{short}_index AS
            SELECT * FROM p1_{short}_eligible_visits WHERE is_new_episode
            """
        )


def prepare_bp_slots(con):
    con.execute(
        """
        CREATE OR REPLACE TEMP TABLE p1_bp_slots AS
        WITH mc_components AS (
            SELECT n.CSN AS visit_id,
                CASE upper(trim(n.Source)) WHEN 'CHART' THEN 'Chart' ELSE 'Monitor' END AS source,
                upper(trim(n.Measure)) AS measure,
                try_cast(n."Time" AS TIMESTAMP) AS reading_ts, n."Time" AS raw_time,
                CASE WHEN isfinite(try_cast(n.Value AS DOUBLE)) THEN try_cast(n.Value AS DOUBLE)
                END AS value
            FROM mcmed_numerics AS n
            JOIN p1_mcmed_index AS i USING (CSN)
            WHERE upper(trim(n.Source)) IN ('CHART', 'MONITOR')
              AND upper(trim(n.Measure)) IN ('SBP', 'DBP')
        ), mc_paired_slots AS (
            SELECT 'MC-MED' AS site, visit_id, source, reading_ts,
                min(value) FILTER (WHERE measure = 'SBP') AS sbp,
                min(value) FILTER (WHERE measure = 'DBP') AS dbp,
                count(DISTINCT value) FILTER (WHERE measure = 'SBP') > 1
                    OR count(DISTINCT value) FILTER (WHERE measure = 'DBP') > 1 AS conflict
            FROM mc_components
            GROUP BY visit_id, source, reading_ts, CASE WHEN reading_ts IS NULL THEN raw_time END
        ), raw_slots AS (
            SELECT * FROM mc_paired_slots
            UNION ALL
            SELECT 'MC-MED', CSN, 'Triage', arrival_ts,
                try_cast(Triage_SBP AS DOUBLE), try_cast(Triage_DBP AS DOUBLE), false
            FROM p1_mcmed_index
            UNION ALL
            SELECT 'MIMIC-IV-ED', i.stay_id, 'Triage', i.arrival_ts,
                try_cast(t.sbp AS DOUBLE), try_cast(t.dbp AS DOUBLE), false
            FROM mimic_triage AS t
            JOIN p1_mimic_index AS i USING (stay_id)
            UNION ALL
            SELECT 'MIMIC-IV-ED', i.stay_id, 'Vitalsign', try_cast(v.charttime AS TIMESTAMP),
                try_cast(v.sbp AS DOUBLE), try_cast(v.dbp AS DOUBLE), false
            FROM mimic_vitalsign AS v
            JOIN p1_mimic_index AS i USING (stay_id)
        ), index_times AS (
            SELECT 'MC-MED' AS site, CSN AS visit_id, arrival_ts, departure_ts
            FROM p1_mcmed_index
            UNION ALL
            SELECT 'MIMIC-IV-ED', stay_id, arrival_ts, departure_ts
            FROM p1_mimic_index
        ), flags AS (
            SELECT r.*, i.arrival_ts,
                coalesce(isfinite(r.sbp) AND isfinite(r.dbp), false) AND NOT r.conflict AS paired,
                coalesce(r.reading_ts BETWEEN i.arrival_ts AND i.departure_ts
                    AND r.reading_ts <= i.arrival_ts + $max_hours * INTERVAL '1 hour', false)
                    AS within_ed,
                r.sbp NOT BETWEEN $sbp_min AND $sbp_max AS sbp_limit_failure,
                r.dbp NOT BETWEEN $dbp_min AND $dbp_max AS dbp_limit_failure,
                r.sbp - r.dbp < $pulse_pressure_min AS pulse_pressure_failure
            FROM raw_slots AS r
            JOIN index_times AS i USING (site, visit_id)
        )
        SELECT *,
            NOT (sbp_limit_failure OR dbp_limit_failure OR pulse_pressure_failure) AS plausible
        FROM flags
        """,
        {
            "max_hours": cfg.ED_MAX_HOURS,
            "sbp_min": cfg.BP_SBP_MIN,
            "sbp_max": cfg.BP_SBP_MAX,
            "dbp_min": cfg.BP_DBP_MIN,
            "dbp_max": cfg.BP_DBP_MAX,
            "pulse_pressure_min": cfg.BP_MIN_PULSE_PRESSURE,
        },
    )


def deduplicate_bp(con, near_minutes=cfg.G0_BP_DUP_MINUTES):
    con.execute(
        """
        CREATE OR REPLACE TEMP TABLE p1_bp_deduplication AS
        WITH usable AS MATERIALIZED (
            SELECT site, visit_id, source, reading_ts, arrival_ts, sbp, dbp
            FROM p1_bp_slots
            WHERE paired AND within_ed AND plausible
        ), rule_a AS (
            SELECT b.*,
                b.site = 'MC-MED' AND b.source = 'Chart'
                AND b.reading_ts <= b.arrival_ts + $exact_minutes * INTERVAL '1 minute'
                AND EXISTS (
                    SELECT 1 FROM usable AS t
                    WHERE t.site = b.site AND t.visit_id = b.visit_id AND t.source = 'Triage'
                      AND t.sbp = b.sbp AND t.dbp = b.dbp
                ) AS removed_a
            FROM usable AS b
        ), rule_b AS (
            SELECT a.*,
                NOT a.removed_a AND a.site = 'MC-MED' AND a.source = 'Chart'
                AND EXISTS (
                    SELECT 1 FROM usable AS m
                    WHERE m.site = a.site AND m.visit_id = a.visit_id AND m.source = 'Monitor'
                      AND m.reading_ts BETWEEN a.reading_ts - $near_minutes * INTERVAL '1 minute'
                                           AND a.reading_ts + $near_minutes * INTERVAL '1 minute'
                      AND abs(a.sbp - m.sbp) <= $near_sbp
                ) AS removed_b
            FROM rule_a AS a
        )
        SELECT b.*,
            NOT b.removed_a AND NOT b.removed_b AND b.source <> 'Triage'
            AND EXISTS (
                SELECT 1 FROM usable AS t
                WHERE t.site = b.site AND t.visit_id = b.visit_id AND t.source = 'Triage'
                  AND t.reading_ts BETWEEN b.reading_ts - $near_minutes * INTERVAL '1 minute'
                                       AND b.reading_ts + $near_minutes * INTERVAL '1 minute'
                  AND abs(b.sbp - t.sbp) <= $near_sbp
            ) AS removed_c
        FROM rule_b AS b
        """,
        {
            "exact_minutes": cfg.BP_CHART_TRIAGE_DUP_MINUTES,
            "near_minutes": near_minutes,
            "near_sbp": cfg.BP_CROSS_SOURCE_MAX_SBP_DIFF,
        },
    )
    for short, site, visit_id, _, disposition_field, home_config in SITES:
        con.execute(
            f"""
            CREATE OR REPLACE TEMP TABLE p1_bp_{short} AS
            SELECT visit_id, source, reading_ts,
                date_diff('microsecond', arrival_ts, reading_ts) / 60000000.0 AS minutes_from_arrival,
                sbp, dbp
            FROM p1_bp_deduplication
            WHERE site = $site AND NOT (removed_a OR removed_b OR removed_c)
            """,
            {"site": site},
        )
        con.execute(
            f"""
            CREATE OR REPLACE TEMP TABLE p1_{short}_index AS
            WITH counts AS (
                SELECT visit_id, count(*) AS bp_readings,
                    date_diff('microsecond', min(reading_ts), max(reading_ts)) / 60000000.0
                        AS bp_span_minutes
                FROM p1_bp_{short}
                GROUP BY visit_id
            ), flags AS (
                SELECT i.*, coalesce(c.bp_readings, 0) AS bp_readings, c.bp_span_minutes,
                    coalesce(c.bp_readings >= $min_readings AND c.bp_span_minutes >= $min_span, false)
                        AS is_c1,
                    coalesce(i.{disposition_field} = $home, false)
                        AS direct_home,
                    coalesce(upper(trim(i.{disposition_field})) IN (SELECT unnest($incomplete)), false)
                        AS incomplete_visit
                FROM p1_{short}_index AS i
                LEFT JOIN counts AS c ON c.visit_id = i.{visit_id}
            )
            SELECT *, is_c1 AND direct_home AND NOT incomplete_visit AS is_c2
            FROM flags
            """,
            {
                "min_readings": cfg.C1_MIN_BP_READINGS,
                "min_span": cfg.C1_MIN_SPAN_MINUTES,
                "home": getattr(cfg, home_config),
                "incomplete": [value.upper() for value in cfg.INCOMPLETE_DISPOSITIONS],
            },
        )


def audit_flow(con):
    selections = []
    for short, site, _, patient_id, _, _ in SITES:
        stages = [
            (1, "All ED visits", "mcmed_visits" if short == "mcmed" else "mimic_edstays", "true"),
            (2, "Primary diagnosis or chief complaint", f"p1_{short}_screened", "true"),
            (3, "Adult age available and eligible", f"p1_{short}_adult_episodes", "true"),
            (4, "Available visit and patient identifiers", f"p1_{short}_adult_episodes",
             "exclusion_reason IS DISTINCT FROM 'Identifier unavailable'"),
            (5, "After trauma exclusion", f"p1_{short}_adult_episodes",
             "exclusion_reason IS NULL OR exclusion_reason IN "
             "('Invalid ED times', 'Postoperative or postprocedural')"),
            (6, "After invalid ED times exclusion", f"p1_{short}_adult_episodes",
             "exclusion_reason IS NULL OR exclusion_reason = 'Postoperative or postprocedural'"),
            (7, "After postoperative or postprocedural exclusion", f"p1_{short}_eligible_visits", "true"),
            (8, "New episodes (index)", f"p1_{short}_index", "true"),
            (9, "C1: eligible blood-pressure data", f"p1_{short}_index", "is_c1"),
            (10, "C2: C1 with direct discharge home", f"p1_{short}_index", "is_c2"),
            (11, "Sensitivity: first eligible visit", f"p1_{short}_eligible_visits", "is_first_eligible_visit"),
        ]
        selections.extend(
            f"""SELECT '{site}' AS site, {number} AS step, '{name}' AS stage,
                   count(*) AS visits, count(DISTINCT nullif(trim({patient_id}), '')) AS patients
                   FROM {table} WHERE {predicate}"""
            for number, name, table, predicate in stages
        )
    union = " UNION ALL ".join(selections)
    write_aggregate(
        con,
        "p1_flow_counts",
        f"""
        WITH flow AS ({union})
        SELECT *,
            CASE WHEN step BETWEEN 2 AND 10 THEN lag(visits) OVER w - visits END AS visits_removed,
            CASE WHEN step BETWEEN 2 AND 10 THEN lag(patients) OVER w - patients END AS patients_removed
        FROM flow
        WINDOW w AS (PARTITION BY site ORDER BY step)
        ORDER BY site, step
        """,
    )


def audit_bp(con):
    write_aggregate(
        con,
        "p1_bp_quality",
        """
        WITH sources AS (
            SELECT * FROM (VALUES
                ('MC-MED', 'Triage'), ('MC-MED', 'Chart'), ('MC-MED', 'Monitor'),
                ('MIMIC-IV-ED', 'Triage'), ('MIMIC-IV-ED', 'Vitalsign')
            ) AS source_list(site, source)
        ), quality AS (
            SELECT s.site, s.source,
                count(r.visit_id) AS raw_slots,
                count(*) FILTER (WHERE r.conflict) AS conflicting_slots,
                count(*) FILTER (WHERE NOT r.paired AND NOT r.conflict) AS missing_or_nonnumeric_slots,
                count(*) FILTER (WHERE r.paired) AS paired_slots,
                count(*) FILTER (WHERE r.paired AND r.reading_ts IS NULL) AS missing_timestamp_slots,
                count(*) FILTER (WHERE r.paired AND r.reading_ts IS NOT NULL AND NOT r.within_ed)
                    AS outside_ed_slots,
                count(*) FILTER (WHERE r.paired AND r.within_ed) AS paired_ed_slots,
                count(*) FILTER (WHERE r.paired AND r.within_ed AND r.sbp_limit_failure)
                    AS sbp_limit_failures,
                count(*) FILTER (WHERE r.paired AND r.within_ed AND r.dbp_limit_failure)
                    AS dbp_limit_failures,
                count(*) FILTER (WHERE r.paired AND r.within_ed AND r.pulse_pressure_failure)
                    AS pulse_pressure_failures,
                count(*) FILTER (WHERE r.paired AND r.within_ed AND NOT r.plausible)
                    AS plausibility_failure_slots,
                count(*) FILTER (WHERE r.paired AND r.within_ed AND r.plausible) AS valid_before_dedup
            FROM sources AS s
            LEFT JOIN p1_bp_slots AS r USING (site, source)
            GROUP BY s.site, s.source
        ), dedup AS (
            SELECT s.site, s.source,
                count(*) FILTER (WHERE d.removed_a) AS rule_a_removed,
                count(*) FILTER (WHERE d.removed_b) AS rule_b_removed,
                count(*) FILTER (WHERE d.removed_c) AS rule_c_removed,
                count(d.visit_id) FILTER (WHERE NOT (d.removed_a OR d.removed_b OR d.removed_c))
                    AS retained_readings
            FROM sources AS s
            LEFT JOIN p1_bp_deduplication AS d USING (site, source)
            GROUP BY s.site, s.source
        )
        SELECT q.*,
            round(100.0 * sbp_limit_failures / nullif(paired_ed_slots, 0), 4) AS sbp_limit_failure_pct,
            round(100.0 * dbp_limit_failures / nullif(paired_ed_slots, 0), 4) AS dbp_limit_failure_pct,
            round(100.0 * pulse_pressure_failures / nullif(paired_ed_slots, 0), 4)
                AS pulse_pressure_failure_pct,
            round(100.0 * plausibility_failure_slots / nullif(paired_ed_slots, 0), 4)
                AS plausibility_failure_pct,
            d.* EXCLUDE (site, source)
        FROM quality AS q
        JOIN dedup AS d USING (site, source)
        ORDER BY q.site, q.source
        """,
    )


def run(con):
    prepare_candidates(con)
    prepare_episode_cohorts(con)
    prepare_bp_slots(con)
    deduplicate_bp(con)
    audit_flow(con)
    audit_bp(con)
