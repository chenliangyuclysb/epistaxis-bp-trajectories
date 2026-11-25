"""Export aggregate cohort totals, descriptive Table 1 and missingness."""

from math import sqrt
from pathlib import Path

import config as cfg


TABLES = Path("outputs/tables")
POPULATIONS = {"Index": "true", "C1": "is_c1", "C2": "is_c2", "C3": "is_c3"}
CONTINUOUS = (
    "age_years", "triage_pain_score", "triage_heart_rate",
    "hemoglobin", "platelet_count", "inr",
)
CATEGORICAL = ("sex", "race", "ethnicity", "arrival_mode", "triage_acuity_group")
BINARY = (
    "prior_htn", "ckd", "liver_disease", "hht", "prior_epistaxis_365d",
    "anticoagulant_use", "antiplatelet_use", "any_antithrombotic_use",
    "has_medication_list", "hemostasis_level2",
)
STRUCTURALLY_ABSENT_MIMIC = {"ethnicity", "hemoglobin", "platelet_count", "inr", "hemostasis_level2"}


def percent(n, total):
    return round(100 * n / total, 4) if total else None


def wilson(n, total):
    if not total:
        return None, None
    z = cfg.WILSON_95_Z
    p = n / total
    adjustment = 1 + z * z / total
    center = (p + z * z / (2 * total)) / adjustment
    half_width = z * sqrt((p * (1 - p) + z * z / (4 * total)) / total) / adjustment
    return round(100 * max(0, center - half_width), 4), round(100 * min(1, center + half_width), 4)


def cohort_counts(con, write_csv):
    data = []
    for site, source in (("MC-MED", "p1_mcmed_cohort"), ("MIMIC-IV-ED", "p1_mimic_cohort")):
        for population, condition in POPULATIONS.items():
            visits, patients, eligible, e7, e30 = con.execute(f"""
                SELECT count(*), count(DISTINCT patient_id),
                    count(*) FILTER (WHERE return_outcome_eligible),
                    count(*) FILTER (WHERE return_7d_event),
                    count(*) FILTER (WHERE return_30d_event)
                FROM {source} WHERE {condition}
            """).fetchone()
            estimates_7d = (percent(e7, eligible), *wilson(e7, eligible))
            estimates_30d = (percent(e30, eligible), *wilson(e30, eligible))
            if population != "C2":
                estimates_7d = estimates_30d = (None, None, None)
            data.append((site, population, visits, patients, eligible,
                         e7, *estimates_7d, e30, *estimates_30d))
    write_csv(
        TABLES / "p1_cohort_counts.csv",
        ("site", "cohort", "visits", "patients", "return_outcome_eligible_visits",
         "return_7d_events", "return_7d_pct", "return_7d_wilson95_lo_pct", "return_7d_wilson95_hi_pct",
         "return_30d_events", "return_30d_pct", "return_30d_wilson95_lo_pct", "return_30d_wilson95_hi_pct"),
        data,
    )


def extra_outcomes(con, write_csv):
    # Summary denominators are outcome-specific eligibility, with no stratification.
    definitions = (
        ("MC-MED", "p1_mcmed_cohort", "Landmark hemostasis escalation", "landmark_eligible",
         "hemostasis_level3_after_60min"),
        ("MIMIC-IV-ED", "p1_mimic_cohort", "Inpatient ENT hemostasis procedure",
         "inpatient_ent_outcome_eligible", "inpatient_ent_hemostasis_event"),
    )
    rows = []
    for site, relation, label, eligible, event in definitions:
        n, count = con.execute(f"""
            SELECT count(*) FILTER (WHERE {eligible}),
                count(*) FILTER (WHERE {eligible} AND {event})
            FROM {relation}
        """).fetchone()
        rows.append((site, label, n, count, percent(count, n), *wilson(count, n)))
    write_csv(
        TABLES / "p1_secondary_outcomes.csv",
        ("site", "outcome", "eligible_visits", "events", "event_pct", "wilson95_lo_pct", "wilson95_hi_pct"),
        rows,
    )


def subcohort_exclusions(con, write_csv):
    rows = []
    rules = (
        ("Index to C1", "Fewer than two valid BP readings", "bp_readings < ?"),
        ("Index to C1", "BP span below minimum", "NOT is_c1 AND bp_readings >= ?"),
        ("C1 to C2", "Incomplete ED visit", "is_c1 AND incomplete_visit"),
        ("C1 to C2", "Other than direct discharge home", "is_c1 AND NOT is_c2 AND NOT incomplete_visit"),
        ("C1 to C3", "Documented prior hypertension", "is_c1 AND prior_htn"),
    )
    for site, relation in (("MC-MED", "p1_mcmed_cohort"), ("MIMIC-IV-ED", "p1_mimic_cohort")):
        for transition, reason, condition in rules:
            values = [cfg.C1_MIN_BP_READINGS] if transition == "Index to C1" else []
            visits, patients = con.execute(f"""
                SELECT count(*), count(DISTINCT patient_id)
                FROM {relation} WHERE {condition}
            """, values).fetchone()
            rows.append((site, transition, reason, visits, patients))
    write_csv(
        TABLES / "p1_subcohort_exclusions.csv",
        ("site", "transition", "reason", "excluded_visits", "affected_patients"), rows,
    )


def covariate_summaries(con, write_csv):
    rows, table1 = [], []
    fields = CONTINUOUS + CATEGORICAL + BINARY
    for site, relation in (("MC-MED", "p1_mcmed_cohort"), ("MIMIC-IV-ED", "p1_mimic_cohort")):
        field_counts = ", ".join(f"count({field})" for field in fields)
        for population in ("Index", "C1", "C2"):
            n, *observed = con.execute(f"""
                SELECT count(*), {field_counts}
                FROM {relation} WHERE {POPULATIONS[population]}
            """).fetchone()
            for field, count in zip(fields, observed):
                availability = ("structurally_unavailable"
                                if site == "MIMIC-IV-ED" and field in STRUCTURALLY_ABSENT_MIMIC
                                else "available_in_source")
                rows.append((site, population, field, n, count, n - count,
                             percent(n - count, n), availability))

        n, patients = con.execute(f"""
            SELECT count(*), count(DISTINCT patient_id) FROM {relation} WHERE is_c1
        """).fetchone()
        table1.extend([
            (site, "visits", "", n, n, None, None, None, None, None, None),
            (site, "patients", "", n, patients, None, None, None, None, None, None),
        ])
        for field in CONTINUOUS:
            observed, median, q1, q3 = con.execute(f"""
                SELECT count({field}), median({field}),
                    quantile_cont({field}, 0.25), quantile_cont({field}, 0.75)
                FROM {relation} WHERE is_c1
            """).fetchone()
            table1.append((site, field, "", observed, None, None, None, None, median, q1, q3))
        for field in CATEGORICAL:
            categories = con.execute(f"""
                SELECT {field}, count(*) FROM {relation} WHERE is_c1 AND {field} IS NOT NULL
                GROUP BY {field} ORDER BY {field}
            """).fetchall()
            observed = sum(count for _, count in categories)
            if not categories:
                table1.append((site, field, "", 0, None, None, None, None, None, None, None))
            for label, count in categories:
                table1.append((site, field, label, observed, count, percent(count, observed),
                               *wilson(count, observed), None, None, None))
        for field in BINARY:
            observed, yes = con.execute(f"""
                SELECT count({field}), count(*) FILTER (WHERE {field})
                FROM {relation} WHERE is_c1
            """).fetchone()
            table1.append((site, field, "Yes", observed, yes if observed else None,
                           percent(yes, observed), *wilson(yes, observed), None, None, None))
    write_csv(
        TABLES / "p1_covariate_missingness.csv",
        ("site", "cohort", "covariate", "visits", "observed", "missing", "missing_pct", "source_availability"),
        rows,
    )
    write_csv(
        TABLES / "p1_table1.csv",
        ("site", "characteristic", "level", "denominator", "count", "pct", "wilson95_lo_pct",
         "wilson95_hi_pct", "median", "q1", "q3"),
        table1,
    )


def run(con, write_csv):
    cohort_counts(con, write_csv)
    extra_outcomes(con, write_csv)
    subcohort_exclusions(con, write_csv)
    covariate_summaries(con, write_csv)
