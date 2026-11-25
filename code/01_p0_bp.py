"""Audit candidate-visit blood pressure data without exporting any readings."""

import csv
from pathlib import Path

import config


TABLES_DIR = Path(__file__).resolve().parents[1] / "outputs" / "tables"


def write_table(con, name, query):
    result = con.execute(query)
    columns = [column[0] for column in result.description]
    data = result.fetchall()
    with (TABLES_DIR / f"{name}.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(columns)
        writer.writerows(data)


def prepare_readings(con):
    con.execute(
        """
        CREATE OR REPLACE TEMP TABLE p0_bp_settings AS
        SELECT ? AS sbp_min, ? AS sbp_max, ? AS dbp_min, ? AS dbp_max,
               ? AS min_pulse_pressure, ? AS exact_dup_minutes, ? AS dup_sbp_difference,
               ? AS ed_max_hours, ? AS c1_min_readings, ? AS c1_min_span_minutes
        """,
        [
            config.BP_SBP_MIN,
            config.BP_SBP_MAX,
            config.BP_DBP_MIN,
            config.BP_DBP_MAX,
            config.BP_MIN_PULSE_PRESSURE,
            config.BP_CHART_TRIAGE_DUP_MINUTES,
            config.BP_CROSS_SOURCE_MAX_SBP_DIFF,
            config.ED_MAX_HOURS,
            config.C1_MIN_BP_READINGS,
            config.C1_MIN_SPAN_MINUTES,
        ],
    )
    con.execute("CREATE OR REPLACE TEMP TABLE p0_bp_windows (window_minutes INTEGER)")
    con.executemany(
        "INSERT INTO p0_bp_windows VALUES (?)",
        [(window,) for window in config.BP_CROSS_SOURCE_WINDOWS_MINUTES],
    )
    con.execute(
        """
        CREATE OR REPLACE TEMP TABLE p0_mc_bp_components AS
        SELECT
            v.CSN AS visit_id,
            CASE upper(trim(n.Source)) WHEN 'CHART' THEN 'Chart' ELSE 'Monitor' END AS source,
            upper(trim(n.Measure)) AS measure,
            try_cast(n."Time" AS TIMESTAMP) AS measurement_time,
            n."Time" AS raw_measurement_time,
            CASE WHEN isfinite(try_cast(n.Value AS DOUBLE)) THEN try_cast(n.Value AS DOUBLE)
            END AS value
        FROM mcmed_numerics AS n
        JOIN p0_mcmed_candidates AS v ON n.CSN = v.CSN
        WHERE upper(trim(n.Measure)) IN ('SBP', 'DBP')
          AND upper(trim(n.Source)) IN ('CHART', 'MONITOR')
        """
    )
    con.execute(
        """
        CREATE OR REPLACE TEMP TABLE p0_mimic_vital_rows AS
        SELECT
            v.stay_id AS visit_id,
            try_cast(b.charttime AS TIMESTAMP) AS measurement_time,
            b.charttime AS raw_measurement_time,
            CASE WHEN isfinite(try_cast(b.sbp AS DOUBLE)) THEN try_cast(b.sbp AS DOUBLE) END AS sbp,
            CASE WHEN isfinite(try_cast(b.dbp AS DOUBLE)) THEN try_cast(b.dbp AS DOUBLE) END AS dbp
        FROM mimic_vitalsign AS b
        JOIN p0_mimic_candidates AS v ON b.stay_id = v.stay_id
        """
    )
    con.execute(
        """
        CREATE OR REPLACE TEMP TABLE p0_bp_stays AS
        SELECT 'MC-MED' AS site, CSN AS visit_id, MRN AS person_id, arrival_ts, departure_ts
        FROM p0_mcmed_candidates
        UNION ALL
        SELECT 'MIMIC-IV-ED', stay_id, subject_id, arrival_ts, departure_ts
        FROM p0_mimic_candidates
        """
    )
    con.execute(
        """
        CREATE OR REPLACE TEMP TABLE p0_bp_index_stays AS
        SELECT 'MC-MED' AS site, CSN AS visit_id, MRN AS person_id
        FROM p0_mcmed_index
        UNION ALL
        SELECT 'MIMIC-IV-ED', stay_id, subject_id
        FROM p0_mimic_index
        """
    )
    con.execute(
        """
        CREATE OR REPLACE TEMP TABLE p0_bp_sources AS
        SELECT * FROM (VALUES
            ('MC-MED', 'Triage'),
            ('MC-MED', 'Chart'),
            ('MC-MED', 'Monitor'),
            ('MIMIC-IV-ED', 'Triage'),
            ('MIMIC-IV-ED', 'Vitalsign')
        ) AS sources(site, source)
        """
    )
    con.execute(
        """
        CREATE OR REPLACE TEMP TABLE p0_bp_slots AS
        WITH mc_slots AS (
            SELECT
                visit_id, source, measurement_time,
                min(value) FILTER (WHERE measure = 'SBP') AS sbp,
                min(value) FILTER (WHERE measure = 'DBP') AS dbp,
                count(DISTINCT value) FILTER (WHERE measure = 'SBP') > 1
                    OR count(DISTINCT value) FILTER (WHERE measure = 'DBP') > 1 AS conflict,
                count(*) FILTER (WHERE measure = 'SBP') AS sbp_records,
                count(*) FILTER (WHERE measure = 'DBP') AS dbp_records
            FROM p0_mc_bp_components
            GROUP BY visit_id, source, measurement_time,
                     CASE WHEN measurement_time IS NULL THEN raw_measurement_time END
        ), raw_slots AS (
            SELECT 'MC-MED' AS site, * FROM mc_slots
            UNION ALL
            SELECT 'MIMIC-IV-ED', visit_id, 'Vitalsign', measurement_time, sbp, dbp,
                   false, 1, 1
            FROM p0_mimic_vital_rows
            UNION ALL
            SELECT 'MC-MED', CSN, 'Triage', arrival_ts,
                   try_cast(Triage_SBP AS DOUBLE), try_cast(Triage_DBP AS DOUBLE), false, 1, 1
            FROM p0_mcmed_candidates
            UNION ALL
            SELECT 'MIMIC-IV-ED', v.stay_id, 'Triage', v.arrival_ts,
                   try_cast(t.sbp AS DOUBLE), try_cast(t.dbp AS DOUBLE), false, 1, 1
            FROM p0_mimic_candidates AS v
            JOIN mimic_triage AS t ON v.stay_id = t.stay_id
        ), flags AS (
            SELECT
                s.site, s.visit_id, s.source, s.measurement_time, s.sbp, s.dbp, s.conflict,
                s.sbp_records, s.dbp_records, v.arrival_ts, v.departure_ts,
                coalesce(isfinite(s.sbp), false) AS sbp_available,
                coalesce(isfinite(s.dbp), false) AS dbp_available,
                coalesce(s.measurement_time BETWEEN v.arrival_ts AND v.departure_ts
                    AND s.measurement_time <= v.arrival_ts + p.ed_max_hours * INTERVAL '1 hour',
                    false) AS within_ed,
                s.sbp NOT BETWEEN p.sbp_min AND p.sbp_max AS sbp_outside_range,
                s.dbp NOT BETWEEN p.dbp_min AND p.dbp_max AS dbp_outside_range,
                s.sbp - s.dbp < p.min_pulse_pressure AS low_pulse_pressure
            FROM raw_slots AS s
            JOIN p0_bp_stays AS v USING (site, visit_id)
            CROSS JOIN p0_bp_settings AS p
        )
        SELECT *,
               sbp_available AND dbp_available AND NOT conflict AS paired,
               NOT (sbp_outside_range OR dbp_outside_range OR low_pulse_pressure) AS plausible
        FROM flags
        """
    )
    con.execute(
        """
        CREATE OR REPLACE TEMP TABLE p0_bp_reads AS
        SELECT site, visit_id, source, measurement_time, sbp, dbp, arrival_ts
        FROM p0_bp_slots
        WHERE paired AND plausible AND within_ed
        """
    )


def prepare_duplicates(con):
    con.execute(
        """
        CREATE OR REPLACE TEMP TABLE p0_bp_mc_exact AS
        SELECT
            c.*,
            c.measurement_time <= c.arrival_ts + p.exact_dup_minutes * INTERVAL '1 minute'
                AS within_exact_window,
            EXISTS (
                SELECT 1 FROM p0_bp_reads AS t
                WHERE t.site = 'MC-MED' AND t.source = 'Triage' AND t.visit_id = c.visit_id
                  AND t.measurement_time = c.arrival_ts AND t.sbp = c.sbp AND t.dbp = c.dbp
            ) AS identical_triage
        FROM p0_bp_reads AS c
        CROSS JOIN p0_bp_settings AS p
        WHERE c.site = 'MC-MED' AND c.source = 'Chart'
        """
    )
    con.execute(
        """
        CREATE OR REPLACE TEMP TABLE p0_bp_cross_matches AS
        WITH lower_priority AS (
            SELECT site, visit_id, source, measurement_time, sbp, 'Monitor' AS reference_source,
                   within_exact_window AND identical_triage AS eligible_for_exact_removal
            FROM p0_bp_mc_exact
            UNION ALL
            SELECT site, visit_id, source, measurement_time, sbp, 'Triage',
                   within_exact_window AND identical_triage
            FROM p0_bp_mc_exact
            UNION ALL
            SELECT site, visit_id, source, measurement_time, sbp, 'Monitor', false
            FROM p0_bp_reads WHERE site = 'MC-MED' AND source = 'Triage'
            UNION ALL
            SELECT site, visit_id, source, measurement_time, sbp, 'Triage', false
            FROM p0_bp_reads WHERE site = 'MIMIC-IV-ED' AND source = 'Vitalsign'
        )
        SELECT
            c.*, w.window_minutes,
            EXISTS (
                SELECT 1 FROM p0_bp_reads AS r
                WHERE r.site = c.site AND r.visit_id = c.visit_id
                  AND r.source = c.reference_source
                  AND r.measurement_time BETWEEN
                      c.measurement_time - w.window_minutes * INTERVAL '1 minute'
                      AND c.measurement_time + w.window_minutes * INTERVAL '1 minute'
                  AND abs(c.sbp - r.sbp) <= p.dup_sbp_difference
            ) AS has_near_match
        FROM lower_priority AS c
        CROSS JOIN p0_bp_windows AS w
        CROSS JOIN p0_bp_settings AS p
        """
    )


def audit_sources(con):
    write_table(
        con,
        "p0_bp_sources",
        """
        WITH reading_counts AS (
            SELECT site, visit_id, source, count(*) AS n_readings
            FROM p0_bp_reads
            GROUP BY site, visit_id, source
        ), visit_counts AS (
            SELECT
                s.site, s.source, count(*) AS candidate_stays,
                count(c.n_readings) AS stays_with_valid_ed_bp,
                sum(coalesce(c.n_readings, 0)) AS valid_ed_bp,
                quantile_cont(coalesce(c.n_readings, 0), 0.25) AS per_stay_all_q25,
                quantile_cont(coalesce(c.n_readings, 0), 0.50) AS per_stay_all_q50,
                quantile_cont(coalesce(c.n_readings, 0), 0.75) AS per_stay_all_q75,
                quantile_cont(c.n_readings, 0.25) AS per_stay_observed_q25,
                quantile_cont(c.n_readings, 0.50) AS per_stay_observed_q50,
                quantile_cont(c.n_readings, 0.75) AS per_stay_observed_q75
            FROM p0_bp_sources AS s
            JOIN p0_bp_stays AS v USING (site)
            LEFT JOIN reading_counts AS c USING (site, visit_id, source)
            GROUP BY s.site, s.source
        ), adjacent AS (
            SELECT site, source,
                   date_diff('microsecond', lag(measurement_time) OVER (
                       PARTITION BY site, visit_id, source ORDER BY measurement_time
                   ), measurement_time) / 60000000.0 AS gap_minutes
            FROM p0_bp_reads
        ), gaps AS (
            SELECT site, source, count(gap_minutes) AS n_adjacent_gaps,
                   count(*) FILTER (WHERE gap_minutes = 0) AS zero_minute_gaps,
                   min(gap_minutes) AS gap_min,
                   quantile_cont(gap_minutes, 0.25) AS gap_q25,
                   quantile_cont(gap_minutes, 0.50) AS gap_q50,
                   quantile_cont(gap_minutes, 0.75) AS gap_q75,
                   max(gap_minutes) AS gap_max
            FROM adjacent
            GROUP BY site, source
        )
        SELECT v.site, v.source, v.* EXCLUDE (site, source),
               coalesce(g.n_adjacent_gaps, 0) AS n_adjacent_gaps,
               coalesce(g.zero_minute_gaps, 0) AS zero_minute_gaps,
               g.* EXCLUDE (site, source, n_adjacent_gaps, zero_minute_gaps)
        FROM visit_counts AS v
        LEFT JOIN gaps AS g USING (site, source)
        ORDER BY v.site, v.source
        """,
    )
    write_table(
        con,
        "p0_bp_pairing_quality",
        """
        WITH counts AS (
            SELECT
                s.site, s.source, count(b.visit_id) AS candidate_slots,
                count(*) FILTER (WHERE b.sbp_available) AS slots_with_numeric_sbp,
                count(*) FILTER (WHERE b.dbp_available) AS slots_with_numeric_dbp,
                count(*) FILTER (WHERE b.sbp_records > 1 OR b.dbp_records > 1)
                    AS slots_with_repeated_components,
                count(*) FILTER (WHERE b.conflict) AS slots_with_conflicting_components,
                count(*) FILTER (WHERE b.paired) AS paired_slots,
                count(*) FILTER (WHERE b.paired AND b.measurement_time IS NULL)
                    AS paired_slots_without_parsed_time,
                count(*) FILTER (
                    WHERE b.paired AND b.measurement_time IS NOT NULL
                      AND (b.arrival_ts IS NULL OR b.departure_ts IS NULL)
                ) AS paired_slots_without_ed_bounds,
                count(*) FILTER (
                    WHERE b.paired AND b.measurement_time IS NOT NULL
                      AND b.arrival_ts IS NOT NULL AND b.departure_ts IS NOT NULL
                      AND NOT b.within_ed
                ) AS paired_slots_outside_ed,
                count(*) FILTER (WHERE b.paired AND b.within_ed) AS paired_ed_slots,
                count(*) FILTER (WHERE b.paired AND b.within_ed AND b.sbp_outside_range)
                    AS sbp_range_failures,
                count(*) FILTER (WHERE b.paired AND b.within_ed AND b.dbp_outside_range)
                    AS dbp_range_failures,
                count(*) FILTER (WHERE b.paired AND b.within_ed AND b.low_pulse_pressure)
                    AS pulse_pressure_failures,
                count(*) FILTER (WHERE b.paired AND b.within_ed AND NOT b.plausible)
                    AS any_plausibility_failures,
                count(*) FILTER (WHERE b.paired AND b.within_ed AND b.plausible) AS valid_ed_bp
            FROM p0_bp_sources AS s
            LEFT JOIN p0_bp_slots AS b USING (site, source)
            GROUP BY s.site, s.source
        )
        SELECT *,
               round(100.0 * paired_slots / nullif(candidate_slots, 0), 4) AS paired_pct,
               round(100.0 * paired_slots / nullif(slots_with_numeric_sbp, 0), 4)
                   AS paired_pct_of_numeric_sbp,
               round(100.0 * sbp_range_failures / nullif(paired_ed_slots, 0), 4)
                   AS sbp_range_failure_pct,
               round(100.0 * dbp_range_failures / nullif(paired_ed_slots, 0), 4)
                   AS dbp_range_failure_pct,
               round(100.0 * pulse_pressure_failures / nullif(paired_ed_slots, 0), 4)
                   AS pulse_pressure_failure_pct,
               round(100.0 * any_plausibility_failures / nullif(paired_ed_slots, 0), 4)
                   AS any_plausibility_failure_pct
        FROM counts
        ORDER BY site, source
        """,
    )


def audit_duplicates(con):
    write_table(
        con,
        "p0_bp_exact_duplicates",
        """
        WITH counts AS (
            SELECT count(*) AS chart_reads,
                   count(*) FILTER (WHERE c.within_exact_window) AS chart_reads_in_window,
                   count(*) FILTER (WHERE c.within_exact_window AND c.identical_triage)
                       AS chart_reads_equal_to_triage_in_window,
                   count(DISTINCT visit_id) FILTER (WHERE c.within_exact_window AND c.identical_triage)
                       AS stays_with_exact_match
            FROM p0_bp_mc_exact AS c
        )
        SELECT
            p.exact_dup_minutes AS window_minutes, counts.*,
            round(100.0 * chart_reads_equal_to_triage_in_window / nullif(chart_reads_in_window, 0), 4)
                AS exact_match_pct_in_window,
            round(100.0 * chart_reads_equal_to_triage_in_window / nullif(chart_reads, 0), 4)
                AS exact_match_pct_all_chart
        FROM counts
        CROSS JOIN p0_bp_settings AS p
        """,
    )
    write_table(
        con,
        "p0_bp_cross_source_duplicates",
        """
        WITH comparisons AS (
            SELECT * FROM (VALUES
                ('MC-MED', 'Chart', 'Monitor', 'All valid ED readings'),
                ('MC-MED', 'Chart', 'Monitor', 'After exact Chart/Triage removal'),
                ('MC-MED', 'Chart', 'Triage', 'All valid ED readings'),
                ('MC-MED', 'Chart', 'Triage', 'After exact Chart/Triage removal'),
                ('MC-MED', 'Triage', 'Monitor', 'All valid ED readings'),
                ('MIMIC-IV-ED', 'Vitalsign', 'Triage', 'All valid ED readings')
            ) AS comparisons(site, source, reference_source, denominator)
        ), counts AS (
            SELECT
                c.site, c.source, c.reference_source, c.denominator, w.window_minutes,
                count(m.visit_id) AS readings,
                count(*) FILTER (WHERE m.has_near_match) AS matched_readings,
                count(DISTINCT m.visit_id) AS stays,
                count(DISTINCT m.visit_id) FILTER (WHERE m.has_near_match) AS matched_stays
            FROM comparisons AS c
            CROSS JOIN p0_bp_windows AS w
            LEFT JOIN p0_bp_cross_matches AS m
              ON (c.site, c.source, c.reference_source, w.window_minutes) =
                 (m.site, m.source, m.reference_source, m.window_minutes)
             AND (c.denominator <> 'After exact Chart/Triage removal'
                  OR NOT m.eligible_for_exact_removal)
            GROUP BY c.site, c.source, c.reference_source, c.denominator, w.window_minutes
        )
        SELECT *,
               round(100.0 * matched_readings / nullif(readings, 0), 4) AS matched_readings_pct,
               round(100.0 * matched_stays / nullif(stays, 0), 4) AS matched_stays_pct
        FROM counts
        ORDER BY site, source, reference_source, denominator, window_minutes
        """,
    )


def audit_c1(con):
    write_table(
        con,
        "p0_bp_c1",
        """
        WITH rules AS (
            SELECT 'MC-MED' AS site, 'No deduplication' AS rule, NULL::INTEGER AS window_minutes
            UNION ALL
            SELECT 'MC-MED', 'Exact Chart/Triage removal', NULL::INTEGER
            UNION ALL
            SELECT 'MC-MED', 'Exact Chart/Triage and near Chart/Monitor removal', window_minutes
            FROM p0_bp_windows
            UNION ALL
            SELECT 'MC-MED', 'Also remove Triage near Monitor', window_minutes
            FROM p0_bp_windows
            UNION ALL
            SELECT 'MIMIC-IV-ED', 'No deduplication', NULL::INTEGER
            UNION ALL
            SELECT 'MIMIC-IV-ED', 'Remove Vitalsign near Triage', window_minutes
            FROM p0_bp_windows
        ), visits AS (
            SELECT site, visit_id, person_id, 'Diagnosis or chief complaint candidates' AS population
            FROM p0_bp_stays
            UNION ALL
            SELECT site, visit_id, person_id, 'Provisional eligible index stays'
            FROM p0_bp_index_stays
        ), retained AS (
            SELECT
                r.site, r.rule, r.window_minutes, b.visit_id, b.measurement_time
            FROM rules AS r
            JOIN p0_bp_reads AS b USING (site)
            WHERE (
                r.rule = 'No deduplication'
                OR b.site <> 'MC-MED'
                OR b.source <> 'Chart'
                OR NOT EXISTS (
                    SELECT 1 FROM p0_bp_mc_exact AS e
                    WHERE e.visit_id = b.visit_id AND e.measurement_time = b.measurement_time
                      AND e.within_exact_window AND e.identical_triage
                )
            )
            AND NOT EXISTS (
                SELECT 1 FROM p0_bp_cross_matches AS m
                WHERE m.site = b.site AND m.source = b.source AND m.visit_id = b.visit_id
                  AND m.measurement_time = b.measurement_time AND m.sbp = b.sbp
                  AND m.window_minutes = r.window_minutes AND m.has_near_match
                  AND (m.site <> 'MC-MED'
                       OR m.reference_source = 'Monitor' AND (
                           m.source = 'Chart' OR r.rule = 'Also remove Triage near Monitor'
                       ))
            )
        ), per_stay AS (
            SELECT site, rule, window_minutes, visit_id, count(*) AS reading_count,
                   date_diff('microsecond', min(measurement_time), max(measurement_time))
                       / 60000000.0 AS span_minutes
            FROM retained
            GROUP BY site, rule, window_minutes, visit_id
        ), qualifies AS (
            SELECT *
            FROM per_stay
            CROSS JOIN p0_bp_settings
            WHERE reading_count >= c1_min_readings AND span_minutes >= c1_min_span_minutes
        ), populations AS (
            SELECT * FROM (VALUES
                ('Diagnosis or chief complaint candidates'),
                ('Provisional eligible index stays')
            ) AS populations(population)
        )
        SELECT
            r.site, p.population, r.rule, r.window_minutes,
            CASE WHEN r.site = 'MIMIC-IV-ED' THEN 'Adult eligibility unverified'
                 WHEN p.population = 'Provisional eligible index stays' THEN 'Adult eligibility verified'
                 ELSE 'No age filter on candidate pool'
            END AS age_scope,
            count(v.visit_id) AS stays,
            count(DISTINCT v.person_id) AS patients,
            count(q.visit_id) AS c1_stays,
            count(DISTINCT v.person_id) FILTER (WHERE q.visit_id IS NOT NULL) AS c1_patients,
            round(100.0 * count(q.visit_id) / nullif(count(v.visit_id), 0), 4) AS c1_stays_pct
        FROM rules AS r
        CROSS JOIN populations AS p
        LEFT JOIN visits AS v ON r.site = v.site AND p.population = v.population
        LEFT JOIN qualifies AS q
          ON (r.site, r.rule, v.visit_id) = (q.site, q.rule, q.visit_id)
         AND r.window_minutes IS NOT DISTINCT FROM q.window_minutes
        GROUP BY r.site, p.population, r.rule, r.window_minutes
        ORDER BY r.site, p.population, r.rule, r.window_minutes
        """,
    )


def run(con):
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    prepare_readings(con)
    prepare_duplicates(con)
    audit_sources(con)
    audit_duplicates(con)
    audit_c1(con)
