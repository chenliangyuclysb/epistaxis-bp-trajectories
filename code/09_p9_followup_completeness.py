"""Audit aggregate ED date ranges and nominal administrative censoring for task 9.2."""

from datetime import date, timedelta
from importlib import import_module
from pathlib import Path
import sys

import duckdb

import config as cfg


OUTPUT = Path("outputs/tables/p9_followup_completeness.csv")
SOURCES = (
    ("MC-MED", "mcmed", Path("mc-med/visits.csv"), "Arrival_time", "Departure_time"),
    (
        "MIMIC-IV-ED", "mimic", Path("data_raw/mimic-iv-ed/ed/edstays.csv.gz"),
        "intime", "outtime",
    ),
)
YEAR_METRICS = (
    "visits", "arrival_parsed", "arrival_distinct_shifted_years",
    "arrival_min_shifted_year", "arrival_max_shifted_year",
    "departure_parsed", "departure_distinct_shifted_years",
    "departure_min_shifted_year", "departure_max_shifted_year",
)
HEADER = (
    "site", "ed_source", "locked_return_cohort_source",
    "official_documentation", "date_shift_scope", "common_shifted_cutoff_applicable",
    "seven_day_sensitivity_refit_status", "documented_coverage",
    "coverage_assumption", "nominal_coverage_days",
    *(f"all_ed_{metric}" for metric in YEAR_METRICS),
    *(f"return_cohort_{metric}" for metric in YEAR_METRICS),
    "return_window_days", "nominal_truncation_fraction_if_uniform",
    "nominal_truncation_percent_if_uniform", "scenario_interpretation",
)


def year_counts(con, source, arrival, departure):
    return con.execute(f"""
        WITH shifted_years AS (
            SELECT year(try_cast({arrival} AS TIMESTAMP)) AS arrival_year,
                year(try_cast({departure} AS TIMESTAMP)) AS departure_year
            FROM {source}
        )
        SELECT count(*), count(arrival_year), count(DISTINCT arrival_year),
            min(arrival_year), max(arrival_year),
            count(departure_year), count(DISTINCT departure_year),
            min(departure_year), max(departure_year)
        FROM shifted_years
    """).fetchone()


def source_reader(p0, path):
    columns = ", ".join(f"{p0.literal(name)}: 'VARCHAR'" for name in p0.csv_columns(path))
    return (
        f"read_csv({p0.literal(path.as_posix())}, columns = {{{columns}}}, "
        "auto_detect = false, header = true, delim = ',', strict_mode = true)"
    )


def audit(con, p0):
    rows = []
    for site, short, source_file, arrival, departure in SOURCES:
        cohort_file = Path(f"derived/cohort_{short}.parquet")
        source = source_reader(p0, source_file)
        cohort = f"(SELECT arrival_ts, departure_ts FROM read_parquet({p0.literal(cohort_file.as_posix())}) WHERE is_c2)"
        all_ed = year_counts(con, source, p0.ident(arrival), p0.ident(departure))
        return_cohort = year_counts(con, cohort, "arrival_ts", "departure_ts")
        if not (0 < return_cohort[0] <= all_ed[0]
                and return_cohort[0] == return_cohort[1] == return_cohort[5]):
            raise ValueError("Unexpected empty or incomplete locked return cohort.")

        start, end, coverage, assumption, documentation = cfg.P9_ED_COVERAGE[site]
        duration = (date(*end) - date(*start)).days
        for hours in (cfg.RETURN_7D_HOURS, cfg.RETURN_30D_HOURS):
            window = timedelta(hours=hours).days
            fraction = window / duration
            rows.append((
                site, source_file.as_posix(), cohort_file.as_posix(),
                documentation, "Patient-specific", False,
                "Not applicable with patient-specific time shifts", coverage,
                assumption, duration, *all_ed, *return_cohort, window,
                round(fraction, 10), round(100 * fraction, 8),
                "Nominal end-of-coverage scenario under uniform index arrivals "
                "and negligible arrival-to-departure lag. This is a scenario, "
                "without estimating censored visits or bounding their true fraction.",
            ))
    return rows


def main():
    try:
        import_module("02_p1_setup").verify_project_runtime()
        p0 = import_module("01_p0_audit")
        with duckdb.connect(":memory:") as con:
            rows = audit(con, p0)
        p0.write_csv(OUTPUT, HEADER, rows)
    except Exception as error:
        # CSV scanner exceptions may contain individual records.
        print(f"9.2 stopped during the date-shift audit ({type(error).__name__}).", file=sys.stderr)
        return 1
    print("9.2 aggregate date-shift and nominal follow-up audit complete.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
