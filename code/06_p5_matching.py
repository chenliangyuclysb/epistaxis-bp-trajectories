"""G4 5.2 risk-set matching and aggregate SBP response estimates."""

from collections import Counter
from dataclasses import dataclass
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm

import config as cfg


TABLES = Path("outputs/tables")
SITES = ("MC-MED", "MIMIC-IV-ED")
SCENARIOS = ("primary", "exclude_b_controls")
MICROSECONDS_PER_MINUTE = 60_000_000


@dataclass(frozen=True)
class Visit:
    site: str
    key: str
    person: str
    acute: bool
    any_aht: bool
    dose: int | None
    times: np.ndarray
    sbp: np.ndarray
    dbp: np.ndarray


@dataclass(frozen=True)
class Anchor:
    visit: Visit
    position: int
    reference: int

    @property
    def sbp(self):
        return self.visit.sbp[self.position]

    @property
    def dbp(self):
        return self.visit.dbp[self.position]


@dataclass(frozen=True)
class Match:
    treated: Anchor
    controls: tuple[Anchor, ...]


@dataclass(frozen=True)
class Observation:
    anchor: Anchor
    position: int

    @property
    def sbp(self):
        return self.anchor.visit.sbp[self.position]

    @property
    def dbp(self):
        return self.anchor.visit.dbp[self.position]

    @property
    def change(self):
        return self.sbp - self.anchor.sbp


def to_microseconds(minutes):
    return np.rint(np.asarray(minutes) * MICROSECONDS_PER_MINUTE).astype(np.int64)


def make_visits(frame, readings):
    c1 = frame.loc[frame["is_c1"]].sort_values(["site", "visit_id"])
    required = (
        "site", "visit_id", "patient_id", "ed_acute_aht", "ed_any_aht", "ed_duration_minutes",
    )
    if c1[list(required)].isna().any().any() or c1.duplicated(["site", "visit_id"]).any():
        raise ValueError("Unexpected null or duplicate keys in matching cohort.")
    if not c1["site"].isin(SITES).all():
        raise ValueError("Unexpected center in matching cohort.")
    visits = []
    for row in c1.itertuples(index=False):
        records = np.asarray(readings[row.site][row.visit_id], dtype=float)
        if records.ndim != 2 or records.shape[1] != 3 or not np.isfinite(records).all():
            raise ValueError("C1 matching requires finite, paired and ordered BP records.")
        if not np.isfinite(row.ed_duration_minutes) or row.ed_duration_minutes <= 0:
            raise ValueError("C1 matching requires a valid ED observation interval.")
        times = to_microseconds(records[:, 0])
        if np.any(np.diff(times) < 0):
            raise ValueError("C1 matching BP sequence is not time ordered.")
        duration = int(to_microseconds(row.ed_duration_minutes))
        within = (times >= 0) & (times <= duration)
        dose = None
        if pd.notna(row.ed_first_acute_aht_minutes):
            if not np.isfinite(row.ed_first_acute_aht_minutes):
                raise ValueError("Nonfinite class A dose time.")
            dose = int(to_microseconds(row.ed_first_acute_aht_minutes))
            if not row.ed_acute_aht or not 0 <= dose <= duration:
                raise ValueError("Class A dose time is outside the documented ED interval.")
        if row.ed_acute_aht and not row.ed_any_aht:
            raise ValueError("Class B medication flag must include every class A visit.")
        visits.append(
            Visit(row.site, row.visit_id, row.patient_id, row.ed_acute_aht, row.ed_any_aht, dose,
                  times[within], records[within, 1], records[within, 2])
        )
    return visits


def final_anchor(visit):
    if visit.dose is None or not len(visit.times):
        return None
    last = np.searchsorted(visit.times, visit.dose, side="right") - 1
    lookback = int(to_microseconds(cfg.P5_MATCH_ANCHOR_LOOKBACK_MINUTES))
    if last < 0 or visit.times[last] < visit.dose - lookback:
        return None
    # Preserve the same-time priority already assigned by P2.
    position = int(np.searchsorted(visit.times, visit.times[last], side="left"))
    return Anchor(visit, position, visit.dose)


def has_post_reading(visit, reference):
    j = np.searchsorted(visit.times, reference, side="right")
    return j < len(visit.times) and visit.times[j] <= reference + to_microseconds(cfg.P5_MATCH_FOLLOWUP_MINUTES)


def eligible_control_positions(visit):
    times = visit.times
    next_distinct = np.searchsorted(times, times, side="right")
    has_later = next_distinct < len(times)
    next_positions = next_distinct[has_later]
    gap = times[next_positions] - times[has_later]
    positions = np.flatnonzero(has_later)
    return positions[gap <= to_microseconds(cfg.P5_MATCH_FOLLOWUP_MINUTES)]


def best_control_anchor(treated, control, positions):
    times = control.times[positions]
    time_gap = np.abs(times - treated.reference)
    sbp_gap = np.abs(control.sbp[positions] - treated.sbp)
    eligible = (time_gap <= to_microseconds(cfg.P5_MATCH_TIME_CALIPER_MINUTES)) & (
        sbp_gap <= cfg.P5_MATCH_SBP_CALIPER_MMHG
    )
    candidates = np.flatnonzero(eligible)
    if not len(candidates):
        return None
    # For tied readings in a control visit, favor the earlier reading, then the P2 order.
    order = np.lexsort((positions[candidates], times[candidates], time_gap[candidates], sbp_gap[candidates]))
    j = candidates[order[0]]
    anchor = Anchor(control, int(positions[j]), int(times[j]))
    return sbp_gap[j], time_gap[j], anchor


def person_count(anchors):
    return len({(a.visit.site, a.visit.person) for a in anchors})


def match_visits(visits):
    by_scenario = {label: [] for label in SCENARIOS}
    attrition = []
    rng = np.random.default_rng(cfg.RANDOM_SEED)
    for site in SITES:
        pool = [visit for visit in visits if visit.site == site]
        treated = [visit for visit in pool if visit.acute]
        timed = sum(visit.dose is not None for visit in treated)
        untreated = [(visit, eligible_control_positions(visit)) for visit in pool if not visit.acute]
        baseline = [anchor for visit in treated if (anchor := final_anchor(visit)) is not None]
        ready = [anchor for anchor in baseline if has_post_reading(anchor.visit, anchor.reference)]
        eligible_controls = [(visit, positions) for visit, positions in untreated if len(positions)]
        sets = {label: [] for label in SCENARIOS}
        for anchor in ready:
            available = []
            for control, positions in eligible_controls:
                candidate = best_control_anchor(anchor, control, positions)
                if candidate is not None:
                    sbp_gap, time_gap, control_anchor = candidate
                    available.append((sbp_gap, time_gap, rng.random(), control_anchor))
            ranked = [item[-1] for item in sorted(available, key=lambda x: x[:3])]
            control_sets = {
                "primary": ranked,
                "exclude_b_controls": [c for c in ranked if not c.visit.any_aht],
            }
            for scenario, selected in control_sets.items():
                selection = tuple(selected[: cfg.P5_MATCH_MAX_CONTROLS])
                if selection:
                    sets[scenario].append(Match(anchor, selection))
        for label, matched in sets.items():
            by_scenario[label].extend(matched)
            controls = [c for match in matched for c in match.controls]
            frequencies = Counter(c.visit.key for c in controls)
            allow_distances = len(matched) >= 2 and person_count(controls) >= 2
            time_gaps = [
                abs(c.reference - match.treated.reference) / MICROSECONDS_PER_MINUTE
                for match in matched for c in match.controls
            ]
            sbp_gaps = [abs(c.sbp - match.treated.sbp) for match in matched for c in match.controls]
            include = lambda visit: label == "primary" or not visit.any_aht
            attrition.append({
                "scenario": label, "site": site, "c1_visits": len(pool),
                "class_a_treated_visits": len(treated),
                "treated_with_timed_dose": timed,
                "treated_without_timed_dose": len(treated) - timed,
                "treated_with_anchor": len(baseline),
                "timed_treated_without_anchor": timed - len(baseline),
                "treated_with_anchor_and_followup": len(ready),
                "anchored_treated_without_post_reading": len(baseline) - len(ready),
                "treated_without_eligible_match": len(ready) - len(matched),
                "class_a_treated_without_matched_set": len(treated) - len(matched),
                "eligible_control_pool_visits": sum(include(visit) for visit, _ in untreated),
                "control_pool_visits_with_candidate_anchor": sum(include(visit) for visit, _ in eligible_controls),
                "matched_sets": len(matched), "control_slots": len(controls),
                "unique_control_visits": len(frequencies),
                "control_visits_used_in_multiple_sets": sum(n > 1 for n in frequencies.values()),
                "sets_with_one_control": sum(len(m.controls) == 1 for m in matched),
                "sets_with_two_controls": sum(len(m.controls) == 2 for m in matched),
                "sets_with_three_controls": sum(len(m.controls) == 3 for m in matched),
                "median_abs_time_gap_minutes": np.median(time_gaps) if allow_distances else None,
                "median_abs_sbp_gap_mmhg": np.median(sbp_gaps) if allow_distances else None,
            })
    return by_scenario, attrition


def measure_at(anchor, offset, lower, upper):
    if offset == 0:
        return Observation(anchor, anchor.position)
    time = anchor.visit.times
    start = anchor.reference + int(to_microseconds(lower))
    stop = anchor.reference + int(to_microseconds(upper))
    target = anchor.reference + int(to_microseconds(offset))
    lo = np.searchsorted(time, start, side="left")
    hi = np.searchsorted(time, stop, side="right")
    if hi == lo:
        return None
    # The nearest target takes priority. Ties favor the earlier time and then P2 order.
    position = lo + np.argmin(np.abs(time[lo:hi] - target))
    return Observation(anchor, int(position))


def measured_sets(matches, window):
    observed = []
    for match in matches:
        treated = measure_at(match.treated, *window)
        controls = tuple(obs for anchor in match.controls if (obs := measure_at(anchor, *window)) is not None)
        if treated is not None and controls:
            observed.append((treated, controls))
    return observed


def percentile_bootstrap(matrix, rng):
    if len(matrix) < 2:
        return np.full((3, matrix.shape[1]), np.nan)
    draws = rng.integers(0, len(matrix), size=(cfg.P5_MATCH_BOOTSTRAPS, len(matrix)))
    resampled = matrix[draws].mean(axis=1)
    confidence = np.quantile(resampled, cfg.P5_MATCH_CI_QUANTILES, axis=0)
    return np.vstack([matrix.mean(axis=0), confidence])


def mixed_check(observed, key, allow_estimate):
    result = {
        **key, "matched_sets": len(observed), "measurements": sum(1 + len(c) for _, c in observed),
        "treated_minus_control_change": None, "ci95_low": None, "ci95_high": None,
        "estimated_matched_set_variance": None, "residual_variance": None,
        "optimizer": None,
        "status": "Insufficient contributing matched sets or distinct patients",
    }
    if not allow_estimate:
        return result
    y, treatment, groups, sites = [], [], [], []
    for group, (treated, controls) in enumerate(observed):
        for member, active in [(treated, 1), *((c, 0) for c in controls)]:
            y.append(member.change)
            treatment.append(active)
            groups.append(group)
            sites.append(member.anchor.visit.site)
    columns = [np.ones(len(y)), treatment]
    if len(set(sites)) > 1:
        columns.append([site == "MIMIC-IV-ED" for site in sites])
    x = np.column_stack(columns)
    for optimizer in ("powell", "lbfgs", "nm"):
        try:
            with warnings.catch_warnings(record=True) as recorded:
                warnings.simplefilter("always")
                model = sm.MixedLM(y, x, groups=np.asarray(groups)).fit(reml=True, method=optimizer)
                interval = np.asarray(model.conf_int(alpha=2 * cfg.P5_MATCH_CI_QUANTILES[0]))[1]
                estimate = model.fe_params[1]
        except (np.linalg.LinAlgError, ValueError, RuntimeError, ZeroDivisionError) as exc:
            # Statsmodels messages may contain values. Report the error class only.
            result["status"] = f"Mixed model unavailable ({type(exc).__name__})"
            continue
        if not (model.converged and np.isfinite(estimate) and np.isfinite(interval).all()):
            result["status"] = "Mixed model failed to converge or interval unavailable"
            continue
        result.update({
            "treated_minus_control_change": estimate, "ci95_low": interval[0], "ci95_high": interval[1],
            "estimated_matched_set_variance": np.asarray(model.cov_re)[0, 0],
            "residual_variance": model.scale, "optimizer": optimizer,
            "status": "Converged with numerical warnings" if recorded else "Converged",
        })
        break
    return result


def bp_normal(sbp, dbp):
    return sbp < cfg.BP_ELEVATED_SBP and dbp < cfg.BP_ELEVATED_DBP


def recovery_summary(observed, key, rng):
    for choice, get_group in (("Treated", lambda t, c: (t,)), ("Matched controls", lambda t, c: c)):
        for elevated_only in (False, True):
            groups = []
            counts = []
            people = []
            for treated, controls in observed:
                members = get_group(treated, controls)
                if elevated_only:
                    members = tuple(
                        m for m in members if not bp_normal(m.anchor.sbp, m.anchor.dbp)
                    )
                if not members:
                    continue
                normal = sum(bp_normal(m.sbp, m.dbp) for m in members)
                counts.append((len(members), normal))
                people.extend(m.anchor for m in members)
                groups.append(100.0 * normal / len(members))
            measured = sum(a for a, _ in counts)
            allowed = len(groups) >= 2 and person_count(people) >= 2
            interval = (
                percentile_bootstrap(np.asarray(groups)[:, None], rng)[:, 0]
                if allowed else (None, None, None)
            )
            yield {
                **key, "group": choice,
                "population": "Anchor BP elevated" if elevated_only else "All observed",
                "matched_sets": len(groups), "observations": measured,
                "normal_observations": sum(b for _, b in counts) if allowed else None,
                "set_weighted_percent_normal": interval[0], "ci95_low_pct": interval[1],
                "ci95_high_pct": interval[2],
            }


def summarize_matches(scenarios):
    means, changes, recovery, mixed = [], [], [], []
    rng = np.random.default_rng(cfg.RANDOM_SEED)
    for scenario in SCENARIOS:
        matches = scenarios[scenario]
        for site in (*SITES, "Pooled"):
            selected = [match for match in matches if site == "Pooled" or match.treated.visit.site == site]
            for window in ((0, 0, 0), *cfg.P5_MATCH_RESPONSE_WINDOWS):
                minute = window[0]
                key = {"scenario": scenario, "site": site, "minutes": minute}
                measured = measured_sets(selected, window)
                t_anchors = [t.anchor for t, _ in measured]
                c_anchors = [c.anchor for _, controls in measured for c in controls]
                t_count, c_count = len(t_anchors), len(c_anchors)
                distinct_t, distinct_c = person_count(t_anchors), person_count(c_anchors)
                allowed_t = len(measured) >= 2 and distinct_t >= 2
                allowed_c = len(measured) >= 2 and distinct_c >= 2
                allowed_effect = allowed_t and allowed_c
                rows = np.array([
                    (t.sbp, np.mean([c.sbp for c in controls]), t.change,
                     np.mean([c.change for c in controls]))
                    for t, controls in measured
                ]).reshape(-1, 4)
                matrix = np.column_stack([rows, rows[:, 2] - rows[:, 3]])
                estimates = percentile_bootstrap(matrix, rng)
                for j, (group, count, allowed) in enumerate((
                    ("Treated", t_count, allowed_t), ("Matched controls", c_count, allowed_c),
                )):
                    values = estimates[:, j] if allowed else (None, None, None)
                    means.append({
                        **key, "group": group, "matched_sets": len(measured), "readings": count,
                        "mean_sbp": values[0], "ci95_low": values[1], "ci95_high": values[2],
                    })
                if minute:
                    difference = estimates[:, 4] if allowed_effect else (None, None, None)
                    changes.append({
                        **key, "matched_sets": len(measured), "treated_visits": t_count,
                        "control_slots": c_count, "distinct_treated_patients": distinct_t,
                        "distinct_control_patients": distinct_c,
                        "unique_control_visits": len({(a.visit.site, a.visit.key) for a in c_anchors}),
                        "treated_mean_sbp_change": estimates[0, 2] if allowed_t else None,
                        "control_mean_sbp_change": estimates[0, 3] if allowed_c else None,
                        "treated_minus_control_change": difference[0],
                        "ci95_low": difference[1], "ci95_high": difference[2],
                        "bootstrap_reps": cfg.P5_MATCH_BOOTSTRAPS if allowed_effect else 0,
                        "status": "Estimated" if allowed_effect else "Insufficient contributing sets or distinct patients",
                    })
                    mixed.append(mixed_check(measured, key, allowed_effect))
                if minute == cfg.P5_MATCH_RESPONSE_WINDOWS[-1][0]:
                    recovery.extend(recovery_summary(measured, key, rng))
    return {
        "p5_matching_group_means": means, "p5_matched_response": changes,
        "p5_matching_recovery": recovery, "p5_matching_mixedlm": mixed,
    }


def run(frame, readings, tables=TABLES):
    """Produce only aggregates from a C1-capable index frame and ordered in-memory BP arrays."""
    matched, flow = match_visits(make_visits(frame, readings))
    aggregates = summarize_matches(matched)
    aggregates["p5_matching_flow"] = flow
    tables.mkdir(parents=True, exist_ok=True)
    for name, rows in aggregates.items():
        path = tables / f"{name}.csv"
        staged = path.with_suffix(".csv.tmp")
        try:
            pd.DataFrame(rows).to_csv(staged, index=False, float_format="%.8f")
            staged.replace(path)
        finally:
            staged.unlink(missing_ok=True)
