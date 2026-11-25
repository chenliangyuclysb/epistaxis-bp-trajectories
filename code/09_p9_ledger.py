"""Append numeric P9 aggregate results without changing existing ledger bytes."""

import argparse
from importlib import import_module


p8_ledger = import_module("08_p8_ledger")
SOURCES = {
    "outputs/tables/p9_risk_difference.csv": ("analysis", "contrast"),
    "outputs/tables/p9_followup_completeness.csv": ("site", "return_window_days"),
}


def interval_reference(path, line, column, row):
    if column != "risk_difference_pp":
        return None
    low, high = "ci95_low_pp", "ci95_high_pp"
    if p8_ledger.decimal(row[low]) is None or p8_ledger.decimal(row[high]) is None:
        return None
    return {
        "ci95_low": row[low], "ci95_high": row[high],
        "ci95_low_reference": f"{path}:{line}:{low}",
        "ci95_high_reference": f"{path}:{line}:{high}",
    }


def append_results(check_only=False):
    prefix = p8_ledger.LEDGER.read_bytes()
    ledger_path = p8_ledger.LEDGER.relative_to(p8_ledger.ROOT).as_posix()
    fields, previous = p8_ledger.csv_rows(ledger_path)
    ledger = p8_ledger.Ledger(fields, previous.values())
    if len(ledger.items) != len(previous):
        raise ValueError("The existing number ledger contains duplicate keys.")
    existing = {}
    for row in previous.values():
        if row["source_file"] in SOURCES:
            locator = row["source_file"], int(row["source_line"]), row["source_column"]
            if locator in existing:
                raise ValueError("The ledger already contains duplicate P9 result locators.")
            existing[locator] = row
    located = set()
    for path, context_fields in SOURCES.items():
        _, rows = p8_ledger.csv_rows(path)
        if not rows:
            raise ValueError("A P9 aggregate result table is empty.")
        for line, row in rows.items():
            context = ", ".join(row[field] for field in context_fields)
            for column, value in row.items():
                if p8_ledger.decimal(value) is None:
                    continue
                locator = path, line, column
                located.add(locator)
                extra = interval_reference(path, line, column, row)
                if locator in existing:
                    previous_entry = existing[locator]
                    if (previous_entry["source_cell_value"] != value or previous_entry["value"] != value
                            or any(previous_entry[k] != v for k, v in (extra or {}).items())):
                        raise ValueError("A P9 result no longer matches its existing ledger entry.")
                    continue
                label = f"P9 post hoc {context}: {column.replace('_', ' ')}"
                ledger.add(label, value, path, line, column, value, extra)
    if existing.keys() - located:
        raise ValueError("A previously catalogued P9 result is absent from the source tables.")
    if not check_only:
        p8_ledger.write_ledger(prefix, ledger)
        if not p8_ledger.LEDGER.read_bytes().startswith(prefix):
            raise RuntimeError("The original number-ledger bytes did not survive the append.")
    return len(ledger.pending), len(located)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    import_module("02_p1_setup").verify_project_runtime()
    count, checked = append_results(args.check_only)
    action = "pending" if args.check_only else "appended"
    print(f"P9 validated {checked} numeric source cells; {count} entries {action}.")


if __name__ == "__main__":
    main()
