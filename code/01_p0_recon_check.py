"""Recheck the aggregate reconnaissance figures in docs/03 sections 2 and 3."""

import csv
from datetime import timedelta
from pathlib import Path

import config


def number(value, digits=2):
    return "NA" if value is None else f"{value:.{digits}f}"


def proportion(numerator, denominator):
    return 100 * numerator / denominator if denominator else None


def as_comparison(reference, observed, digits=None):
    comparable = round(observed, digits) if digits is not None and observed is not None else observed
    return "Yes" if comparable == reference else "No"


def mcmed_visit_rows(con, add):
    visits, patients, dx, cc, min_age, max_age, adults, los, missing_sbp, missing_dbp, high, severe = con.execute(
        """
        SELECT count(*), count(DISTINCT MRN), count(*) FILTER (WHERE dx_match),
               count(*) FILTER (WHERE cc_match), min(age_years), max(age_years),
               count(*) FILTER (WHERE age_years >= ?), quantile_cont(try_cast(ED_LOS AS DOUBLE), [0.25, 0.5, 0.75]),
               count(*) FILTER (WHERE NOT coalesce(isfinite(try_cast(Triage_SBP AS DOUBLE)), false)),
               count(*) FILTER (WHERE NOT coalesce(isfinite(try_cast(Triage_DBP AS DOUBLE)), false)),
               count(*) FILTER (WHERE isfinite(try_cast(Triage_SBP AS DOUBLE)) AND try_cast(Triage_SBP AS DOUBLE) >= ?),
               count(*) FILTER (WHERE isfinite(try_cast(Triage_SBP AS DOUBLE)) AND try_cast(Triage_SBP AS DOUBLE) >= ?)
        FROM p0_mcmed_candidates
        """,
        [config.MIN_ADULT_AGE, *config.RECON_SBP_CUTOFFS],
    ).fetchone()
    add("MC-MED", "Candidate visits", 447, visits, as_comparison(447, visits), "S3.1 diagnosis or chief complaint before exclusions")
    add("MC-MED", "Candidate patients", 378, patients, as_comparison(378, patients), "Distinct nonmissing MRN among candidates")
    add("MC-MED", "Primary diagnosis matches", 363, dx, as_comparison(363, dx), "Candidate visits whose primary ICD-10 matches")
    add("MC-MED", "Chief complaint matches", 418, cc, as_comparison(418, cc), "Candidate visits whose chief complaint matches")
    add("MC-MED", "Adult candidates and age range", "447 | 18–90", f"{adults} | {number(min_age, 0)}–{number(max_age, 0)}",
        as_comparison((447, 18, 90), (adults, min_age, max_age)), "Age range uses the numeric source age field")

    dispositions = dict(con.execute(
        "SELECT ED_dispo, count(*) FROM p0_mcmed_candidates GROUP BY ED_dispo"
    ).fetchall())
    for label, reference in (("Discharge", 323), ("Observation", 64), ("Inpatient", 58), ("ICU", 2)):
        observed = dispositions.get(label, 0)
        add("MC-MED", f"ED disposition: {label}", reference, observed, as_comparison(reference, observed), "Candidate visits")

    for index, label, reference in ((1, "median", 4.25), (0, "Q1", 2.85), (2, "Q3", 6.43)):
        observed = None if los is None else los[index]
        add("MC-MED", f"ED LOS {label}, h", reference, number(observed), as_comparison(reference, observed, 2), "Source ED_LOS field among candidates")
    add("MC-MED", "Missing triage SBP and DBP", "0 SBP | 0 DBP", f"{missing_sbp} SBP | {missing_dbp} DBP",
        as_comparison((0, 0), (missing_sbp, missing_dbp)), "Missing includes unparseable and nonfinite components")
    for cutoff, n, reference in zip(config.RECON_SBP_CUTOFFS, (high, severe), (55.5, 9.6)):
        observed = proportion(n, visits)
        add("MC-MED", f"Triage SBP ≥{cutoff}, %", reference, f"{number(observed, 1)} ({n}/{visits})",
            as_comparison(reference, observed, 1), "SBP-only reconnaissance threshold among all candidate visits")


def revisit_rows(con, add):
    short_days, long_days = config.RECON_RETURN_DAYS
    short_interval, long_interval = (timedelta(days=days) for days in config.RECON_RETURN_DAYS)
    for site, candidates, visit_id, patient_id, dispo, home_label, references in (
        ("MC-MED", "p0_mcmed_candidates", "CSN", "MRN", "ED_dispo", "Discharge", (37, 34, 51)),
        ("MIMIC", "p0_mimic_candidates", "stay_id", "subject_id", "disposition", "HOME", (None, 134, None)),
    ):
        short_all, short_home, long_all = con.execute(
            f"""
            WITH first_candidate_return AS (
                SELECT i.{visit_id}, i.{dispo}, i.departure_ts, min(r.arrival_ts) AS next_arrival
                FROM {candidates} AS i
                LEFT JOIN {candidates} AS r
                  ON i.{patient_id} = r.{patient_id} AND i.{visit_id} <> r.{visit_id}
                 AND r.arrival_ts > i.departure_ts AND r.arrival_ts <= i.departure_ts + ?
                GROUP BY i.{visit_id}, i.{dispo}, i.departure_ts
            )
            SELECT count(*) FILTER (WHERE next_arrival <= departure_ts + ?),
                   count(*) FILTER (WHERE next_arrival <= departure_ts + ? AND {dispo} = ?),
                   count(next_arrival)
            FROM first_candidate_return
            """,
            [long_interval if references[2] is not None else short_interval, short_interval, short_interval, home_label],
        ).fetchone()
        summaries = [
            (f"{short_days}-day epistaxis return, all candidate visits", references[0], short_all),
            (f"{short_days}-day epistaxis return, home disposition", references[1], short_home),
        ]
        if references[2] is not None:
            summaries.append((f"{long_days}-day epistaxis return, all candidate visits", references[2], long_all))
        for label, reference, observed in summaries:
            if site == "MC-MED" and label.startswith(f"{long_days}-day"):
                immediate = con.execute(
                    """
                    WITH next_ed AS (
                        SELECT CSN, lead(CSN) OVER (
                            PARTITION BY MRN ORDER BY try_cast(Arrival_time AS TIMESTAMP), CSN
                        ) AS next_csn
                        FROM mcmed_visits
                    )
                    SELECT count(*)
                    FROM p0_mcmed_candidates i JOIN next_ed n ON n.CSN = i.CSN
                    JOIN p0_mcmed_candidates r ON r.CSN = n.next_csn
                    WHERE r.arrival_ts > i.departure_ts AND r.arrival_ts <= i.departure_ts + ?
                    """,
                    [long_interval],
                ).fetchone()[0]
                add(site, label, reference, f"Any later {observed}; immediate next ED visit only {immediate}",
                    "Yes, immediate-next definition" if immediate == reference else as_comparison(reference, observed),
                    "Same-patient S3.1 return after departure; the primary count allows intervening visits without epistaxis")
                continue
            add(site, label, reference if reference is not None else "—", observed,
                as_comparison(reference, observed) if reference is not None else "Context",
                "Distinct triggering stays with a later same-patient S3.1 candidate arriving after departure within the stated window")

    paired, from_departure, from_arrival = con.execute(
        """
        WITH next_arrival AS (
            SELECT c.CSN, c.arrival_ts, c.departure_ts, try_cast(c.Hours_to_next_visit AS DOUBLE) AS recorded_hours,
                   min(try_cast(v.Arrival_time AS TIMESTAMP)) AS next_arrival_ts
            FROM p0_mcmed_candidates c
            LEFT JOIN mcmed_visits v ON v.MRN = c.MRN AND v.CSN <> c.CSN
                                   AND try_cast(v.Arrival_time AS TIMESTAMP) > c.arrival_ts
            GROUP BY c.CSN, c.arrival_ts, c.departure_ts, c.Hours_to_next_visit
        )
        SELECT count(*) FILTER (WHERE isfinite(recorded_hours) AND next_arrival_ts IS NOT NULL AND departure_ts IS NOT NULL),
               quantile_cont(abs(recorded_hours - date_diff('second', departure_ts, next_arrival_ts) / 3600.0), [0.5, 0.95])
                   FILTER (WHERE isfinite(recorded_hours) AND departure_ts IS NOT NULL),
               quantile_cont(abs(recorded_hours - date_diff('second', arrival_ts, next_arrival_ts) / 3600.0), [0.5, 0.95])
                   FILTER (WHERE isfinite(recorded_hours) AND departure_ts IS NOT NULL)
        FROM next_arrival
        """
    ).fetchone()
    departure_errors = " / ".join(number(v, 4) for v in (from_departure or (None, None)))
    arrival_errors = " / ".join(number(v, 4) for v in (from_arrival or (None, None)))
    add("MC-MED", "Hours_to_next_visit field consistency", "Pearson r=0.999",
        f"Paired n={paired}; absolute error P50/P95 h: departure {departure_errors}; arrival {arrival_errors}",
        "Correlation not recalculated", "Direct field reconciliation against the next observed same-patient ED visit of any diagnosis")


def mcmed_bp_rows(con, add):
    by_source = con.execute(
        """
        WITH source_sbp AS (
            SELECT visit_id, source, measurement_time,
                   date_diff('second', lag(measurement_time) OVER (
                       PARTITION BY visit_id, source ORDER BY measurement_time NULLS LAST), measurement_time) / 60.0 AS gap_min
            FROM p0_mc_bp_components
            WHERE measure = 'SBP' AND isfinite(value)
        ), raw_counts AS (
            SELECT source, count(*) AS raw_readings
            FROM p0_mc_bp_components WHERE measure = 'SBP' GROUP BY source
        ), source_counts AS (
            SELECT source, visit_id, count(*) AS n FROM source_sbp GROUP BY source, visit_id
        ), count_summary AS (
            SELECT source, count(*) AS visits, sum(n) AS readings,
                   quantile_cont(n, [0.25, 0.5, 0.75]) AS q_counts
            FROM source_counts GROUP BY source
        ), gap_summary AS (
            SELECT source, quantile_cont(gap_min, [0.25, 0.5, 0.75]) AS q_gaps
            FROM source_sbp GROUP BY source
        )
        SELECT source, visits, readings, q_counts, q_gaps, raw_readings
        FROM count_summary JOIN gap_summary USING (source) JOIN raw_counts USING (source)
        WHERE source IN ('Monitor', 'Chart')
        """
    ).fetchall()
    references = {"Monitor": (302, 1971, (3, 5, 9), (8, 15, 29)), "Chart": (372, 1848, (2, 3, 6), (None, 37, None))}
    source_rows = {source: (visits, readings, counts, gaps, raw) for source, visits, readings, counts, gaps, raw in by_source}
    for source, (ref_visits, ref_readings, ref_counts, ref_gaps) in references.items():
        visits, readings, counts, gaps, raw = source_rows.get(source, (0, 0, None, None, 0))
        add("MC-MED", f"{source} visits with SBP", ref_visits, visits, as_comparison(ref_visits, visits),
            "Visits with at least one finite numeric SBP")
        add("MC-MED", f"{source} SBP readings", ref_readings, f"Raw {raw}; finite {readings}",
            as_comparison(ref_readings, raw), "The reference includes raw SBP rows without a usable numeric value")
        for name, units, values, ref_values in (("SBP count", "readings", counts, ref_counts), ("adjacent SBP gap", "min", gaps, ref_gaps)):
            for index, quantile in enumerate(("Q1", "median", "Q3")):
                observed = values[index] if values is not None else None
                reference = ref_values[index]
                if reference is not None:
                    add("MC-MED", f"{source} {name} {quantile}, {units}", reference, number(observed),
                        as_comparison(reference, observed, 0), "Positive counts per visit; elapsed minutes between consecutive timestamped readings within visit/source for gaps")

    same_times, first_chart_equal, first_chart_total, sbp_rows, out_of_ed, unknown_ed = con.execute(
        """
        WITH raw_sbp AS (
            SELECT visit_id, source, measurement_time, value AS sbp
            FROM p0_mc_bp_components
            WHERE measure = 'SBP' AND isfinite(value) AND source IN ('Chart', 'Monitor')
        ), exact_match AS (
            SELECT visit_id, measurement_time FROM raw_sbp WHERE measurement_time IS NOT NULL
            GROUP BY visit_id, measurement_time HAVING count(DISTINCT source) = 2
        ), chart_first AS (
            SELECT visit_id, sbp FROM raw_sbp WHERE source = 'Chart'
            QUALIFY row_number() OVER (PARTITION BY visit_id ORDER BY measurement_time NULLS LAST, sbp) = 1
        ), chart_summary AS (
            SELECT count(*) FILTER (WHERE c.sbp = try_cast(v.Triage_SBP AS DOUBLE)) AS exact_value, count(*) AS comparable
            FROM chart_first c JOIN p0_mcmed_candidates v ON v.CSN = c.visit_id
            WHERE isfinite(try_cast(v.Triage_SBP AS DOUBLE))
        ), bounds AS (
            SELECT count(*) AS n, count(*) FILTER (WHERE measurement_time < arrival_ts OR measurement_time > departure_ts) AS outside,
                   count(*) FILTER (WHERE measurement_time IS NULL OR arrival_ts IS NULL OR departure_ts IS NULL) AS unknown
            FROM raw_sbp b JOIN p0_mcmed_candidates v ON v.CSN = b.visit_id
        )
        SELECT (SELECT count(*) FROM exact_match), exact_value, comparable, n, outside, unknown
        FROM chart_summary CROSS JOIN bounds
        """
    ).fetchone()
    add("MC-MED", "Exact Chart/Monitor timestamp collisions", 0, same_times, as_comparison(0, same_times), "Visit/timestamp pairs containing both raw SBP sources")
    first_match = proportion(first_chart_equal, first_chart_total)
    add("MC-MED", "First Chart SBP equals triage SBP, %", 96,
        f"{number(first_match)} ({first_chart_equal}/{first_chart_total})", as_comparison(96, first_match, 0),
        "Full Chart history; first timestamp, ties ordered by SBP; SBP-only equality without a 10-minute restriction")
    add("MC-MED", "Raw SBP outside arrival–departure", 0, f"{out_of_ed}/{sbp_rows}; unassessable={unknown_ed}",
        as_comparison((0, 0), (out_of_ed, unknown_ed)), "Both source types, inclusive observed arrival and departure, no added ED time cap")

    without_min_span, without_dense, with_min_span, with_dense = con.execute(
        """
        WITH numerics AS (
            SELECT b.visit_id, b.measurement_time, false AS is_triage
            FROM p0_mc_bp_components b JOIN p0_mcmed_candidates v ON v.CSN = b.visit_id
            WHERE b.measure = 'SBP' AND isfinite(b.value) AND b.source IN ('Chart', 'Monitor')
              AND b.measurement_time BETWEEN v.arrival_ts AND v.departure_ts
        ), readings AS (
            SELECT * FROM numerics
            UNION ALL
            SELECT CSN, arrival_ts, true FROM p0_mcmed_candidates
            WHERE isfinite(try_cast(Triage_SBP AS DOUBLE)) AND arrival_ts <= departure_ts
        ), visits AS (
            SELECT visit_id, count(*) AS all_count, count(*) FILTER (WHERE NOT is_triage) AS numeric_count,
                   date_diff('second', min(measurement_time), max(measurement_time)) / 60.0 AS all_span,
                   date_diff('second', min(measurement_time) FILTER (WHERE NOT is_triage),
                                       max(measurement_time) FILTER (WHERE NOT is_triage)) / 60.0 AS numeric_span
            FROM readings GROUP BY visit_id
        )
        SELECT count(*) FILTER (WHERE numeric_count >= ? AND numeric_span >= ?), count(*) FILTER (WHERE numeric_count >= ?),
               count(*) FILTER (WHERE all_count >= ? AND all_span >= ?), count(*) FILTER (WHERE all_count >= ?)
        FROM visits
        """,
        [config.C1_MIN_BP_READINGS, config.C1_MIN_SPAN_MINUTES, config.RECON_MCMED_MIN_DENSE_BP_READINGS] * 2,
    ).fetchone()
    for item, reference, numerics, with_triage in (
        (f"Raw SBP visits with ≥{config.C1_MIN_BP_READINGS} readings and span ≥{config.C1_MIN_SPAN_MINUTES} min", 352, without_min_span, with_min_span),
        (f"Raw SBP visits with ≥{config.RECON_MCMED_MIN_DENSE_BP_READINGS} readings", 287, without_dense, with_dense),
    ):
        status = "Yes" if reference in (numerics, with_triage) else "No"
        add("MC-MED", item, reference, f"Numerics {numerics}; numerics + triage {with_triage}", status,
            "Two source interpretations of raw SBP in the observed ED window, before pairing or deduplication")


def mimic_rows(con, add):
    any_dx, primary_dx = con.execute(
        "SELECT count(*) FILTER (WHERE dx_any_match), count(*) FILTER (WHERE dx_match) FROM p0_mimic_dx_flags"
    ).fetchone()
    stays, patients, cc, union_any, home, los = con.execute(
        """
        WITH ed_any AS (
            SELECT stay_id FROM p0_mimic_candidates
            UNION
            SELECT stay_id FROM p0_mimic_dx_flags WHERE dx_any_match
        )
        SELECT count(*), count(DISTINCT v.subject_id), count(*) FILTER (WHERE v.cc_match),
               (SELECT count(*) FROM ed_any), count(*) FILTER (WHERE v.disposition = 'HOME'),
               median(date_diff('second', v.arrival_ts, v.departure_ts) / 3600.0)
        FROM p0_mimic_candidates AS v
        """
    ).fetchone()
    for item, reference, observed, note in (
        ("Any-position epistaxis diagnosis", 1371, any_dx, "Distinct stays with matching ICD-9 or ICD-10 in any sequence"),
        ("Primary epistaxis diagnosis", 1161, primary_dx, "Distinct stays with matching ICD-9 or ICD-10 in the primary sequence"),
        ("Chief complaint matches", 1318, cc, "Distinct candidates with matching chief complaint"),
        ("Any-position diagnosis or chief complaint union", 1526, union_any, "Distinct stays after union of any diagnosis position and complaint"),
        ("Candidate stays", 1413, stays, "S3.1 primary diagnosis or chief complaint before exclusions"),
        ("Candidate patients", 1127, patients, "Distinct nonmissing subject_id among candidates"),
        ("HOME disposition", 1129, home, "Candidate stays with disposition HOME"),
    ):
        add("MIMIC", item, reference, observed, as_comparison(reference, observed), note)
    add("MIMIC", "ED LOS median, h", 3.8, number(los), as_comparison(3.8, los, 1),
        "Arrival-to-departure hours; truncation to one decimal would yield 3.8")

    counts, counts_including_zero, observed_stays, gap, at_least_first, at_least_second, ed_first, ed_second = con.execute(
        """
        WITH sbp AS (
            SELECT visit_id, measurement_time,
                   date_diff('second', lag(measurement_time) OVER (
                       PARTITION BY visit_id ORDER BY measurement_time NULLS LAST), measurement_time) / 60.0 AS gap_min
            FROM p0_mimic_vital_rows WHERE isfinite(sbp)
        ), counts AS (
            SELECT v.stay_id, count(s.visit_id) AS n,
                   count(*) FILTER (WHERE s.measurement_time BETWEEN v.arrival_ts AND v.departure_ts) AS n_in_ed
            FROM p0_mimic_candidates v LEFT JOIN sbp s ON s.visit_id = v.stay_id
            GROUP BY v.stay_id
        )
        SELECT quantile_cont(n, [0.25, 0.5, 0.75]) FILTER (WHERE n > 0), quantile_cont(n, [0.25, 0.5, 0.75]),
               count(*) FILTER (WHERE n > 0), (SELECT median(gap_min) FROM sbp),
               count(*) FILTER (WHERE n >= ?), count(*) FILTER (WHERE n >= ?),
               count(*) FILTER (WHERE n_in_ed >= ?), count(*) FILTER (WHERE n_in_ed >= ?)
        FROM counts
        """,
        list(config.RECON_MIMIC_BP_READING_COUNTS) * 2,
    ).fetchone()
    for index, label, reference in ((1, "median", 2), (0, "Q1", 1), (2, "Q3", 3)):
        observed = counts[index] if counts is not None else None
        with_zero = counts_including_zero[index] if counts_including_zero is not None else None
        add("MIMIC", f"Vitalsign SBP count {label}", reference,
            f"{number(observed)}; including zero-read stays {number(with_zero)}", as_comparison(reference, observed, 0),
            f"Finite SBP only in vitalsign, {observed_stays} stays with ≥1 record")
    add("MIMIC", "Vitalsign adjacent SBP gap median, min", 103, number(gap), as_comparison(103, gap, 0),
        "Weighted by valid adjacent within-stay finite SBP observations")
    for minimum, reference, observed, in_ed in zip(
        config.RECON_MIMIC_BP_READING_COUNTS, (470, 131), (at_least_first, at_least_second), (ed_first, ed_second)
    ):
        add("MIMIC", f"Stays with ≥{minimum} vitalsign SBP readings", reference, f"Any time {observed}; observed ED window {in_ed}",
            as_comparison(reference, observed), "Finite vitalsign SBP only, excludes the separate triage record")


def drug_pattern(name):
    return config.RECON_DRUG_PATTERNS[name]


def medication_rows(con, add):
    con.execute(
        """
        CREATE OR REPLACE TEMP TABLE p0_recon_mcmed_orders AS
        SELECT o.CSN, o.Order_type, o.Procedure_name, o.First_admin_time, v.arrival_ts, v.departure_ts
        FROM mcmed_orders o JOIN p0_mcmed_candidates v ON v.CSN = o.CSN
        """
    )
    all_total, all_documented = con.execute(
        """
        SELECT count(*), count(*) FILTER (WHERE nullif(trim(First_admin_time), '') IS NOT NULL)
        FROM mcmed_orders WHERE Order_type = 'Medications'
        """
    ).fetchone()
    candidate_total, candidate_documented = con.execute(
        """
        SELECT count(*), count(*) FILTER (WHERE nullif(trim(First_admin_time), '') IS NOT NULL)
        FROM p0_recon_mcmed_orders WHERE Order_type = 'Medications'
        """
    ).fetchone()
    all_rate = proportion(all_documented, all_total)
    candidate_rate = proportion(candidate_documented, candidate_total)
    status = "Yes" if "Yes" in (as_comparison(65.6, all_rate, 1), as_comparison(65.6, candidate_rate, 1)) else "No"
    add("MC-MED", "Medications First_admin_time populated, %", 65.6,
        f"All {number(all_rate, 1)} ({all_documented}/{all_total}); candidates {number(candidate_rate, 1)} ({candidate_documented}/{candidate_total})",
        status, "Order rows with a nonblank administration timestamp; original denominator unspecified")

    for key, reference, meds_only in (
        ("oxymetazoline", 264, True), ("tranexamic_acid", 104, True), ("silver_nitrate", 45, True),
        ("mcmed_epistaxis_supply", 73, False), ("mcmed_ent", None, False), ("mcmed_transfusion", 33, False),
        ("labetalol", 12, True), ("hydralazine", 17, True), ("metoprolol", 10, True),
        ("amlodipine", 7, True), ("acei_arb", 13, True), ("clonidine", 5, True), ("nicardipine", 0, True),
    ):
        patterns = [drug_pattern(name) for name in (("acei", "arb") if key == "acei_arb" else (key,))]
        signals = " OR ".join(["regexp_matches(Procedure_name, ?, 'i')"] * len(patterns))
        observed = con.execute(
            f"""
            SELECT count(DISTINCT CSN) FROM p0_recon_mcmed_orders
            WHERE ({signals}) AND (NOT ? OR Order_type = 'Medications')
            """,
            [*patterns, meds_only],
        ).fetchone()[0]
        agreement = "Context" if reference is None else as_comparison(reference, observed)
        if key == "mcmed_epistaxis_supply" and observed != reference:
            agreement = "Approximate agreement"
        add("MC-MED", f"Candidate visits with order-name match: {key}",
            "Frequent; count unspecified" if reference is None else ("About 73" if key == "mcmed_epistaxis_supply" else reference),
            observed, agreement,
            f"{'Medication orders' if meds_only else 'All order types'} linked to candidates, any order time; original matching dictionary unspecified")

    acute_names = ("labetalol", "hydralazine", "nicardipine", "enalaprilat", "metoprolol", "esmolol")
    signals = " OR ".join(["regexp_matches(Procedure_name, ?, 'i')"] * len(acute_names))
    candidate_count, in_ed_count = con.execute(
        f"""
        SELECT count(DISTINCT CSN), count(DISTINCT CSN) FILTER (
            WHERE try_cast(First_admin_time AS TIMESTAMP) BETWEEN arrival_ts AND departure_ts
        )
        FROM p0_recon_mcmed_orders WHERE Order_type = 'Medications' AND ({signals})
        """,
        [drug_pattern(name) for name in acute_names],
    ).fetchone()
    denominator = con.execute("SELECT count(*) FROM p0_mcmed_candidates").fetchone()[0]
    add("MC-MED", "Acute agent name candidates with in-ED documented administration, %", "About 5%; definition unspecified",
        f"{number(proportion(in_ed_count, denominator), 1)} ({in_ed_count}/{denominator}); {candidate_count} with an order at any time",
        "Definition unspecified", "S8.1 acute-agent names; route unconfirmed and denominator is all candidates")

    con.execute(
        """
        CREATE OR REPLACE TEMP TABLE p0_recon_mimic_pyxis AS
        SELECT p.stay_id, p.name, try_cast(p.charttime AS TIMESTAMP) AS pyxis_time,
               c.arrival_ts, c.departure_ts
        FROM mimic_pyxis p JOIN p0_mimic_candidates c USING (stay_id)
        """
    )
    for key, reference in (
        ("labetalol", 19), ("hydralazine", 14), ("metoprolol", 38), ("amlodipine", 16),
        ("acei", 28), ("nicardipine", 1), ("clonidine", 1), ("oxymetazoline", 350), ("tranexamic_acid", 24),
    ):
        observed, in_ed = con.execute(
            """
            SELECT count(DISTINCT stay_id), count(DISTINCT stay_id) FILTER (WHERE pyxis_time BETWEEN arrival_ts AND departure_ts)
            FROM p0_recon_mimic_pyxis
            WHERE regexp_matches(name, ?, 'i')
            """,
            [drug_pattern(key)],
        ).fetchone()
        add("MIMIC", f"Candidate stays with pyxis name match: {key}", reference, f"Any time {observed}; observed ED window {in_ed}",
            as_comparison(reference, observed), "Counts use generic and brand label candidates; original matching dictionary and denominator unspecified")


def pmh_rows(con, add):
    matching, patients = con.execute(
        """
        WITH patients AS (SELECT DISTINCT MRN FROM p0_mcmed_candidates WHERE MRN IS NOT NULL)
        SELECT count(DISTINCT p.MRN), (SELECT count(*) FROM patients)
        FROM mcmed_pmh p JOIN patients USING (MRN)
        WHERE regexp_matches(trim(p.Code), ?, 'i')
        """,
        [config.MCMED_HTN_ICD10_PATTERN],
    ).fetchone()
    observed = proportion(matching, patients)
    add("MC-MED", "Any-time PMH I10–I16 patients, %", 32.4, f"{number(observed, 1)} ({matching}/{patients})",
        as_comparison(32.4, observed, 1), "ICD-10 pattern without date restrictions; 32.4% cannot result from rounding an integer count out of 378")


def run(con):
    rows = []

    def add(site, metric, reference, observed, agreement, definition):
        rows.append((site, metric, reference, observed, agreement, definition))

    mcmed_visit_rows(con, add)
    revisit_rows(con, add)
    mcmed_bp_rows(con, add)
    mimic_rows(con, add)
    medication_rows(con, add)
    pmh_rows(con, add)

    path = Path("outputs/tables/p0_recon_check.csv")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as output:
        writer = csv.writer(output)
        writer.writerow(("Site", "Metric", "docs/03 value", "Recalculated value", "Agreement", "Audit definition"))
        writer.writerows(row for site in ("MC-MED", "MIMIC") for row in rows if row[0] == site)
