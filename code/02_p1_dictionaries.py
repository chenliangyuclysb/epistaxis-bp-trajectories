"""Compile the G0 dictionaries into reviewable label mappings and SQL lookups."""

import csv
import re
from pathlib import Path

import pyarrow as pa

import config as cfg


DICTIONARY_DIR = Path("outputs/dictionary_final")
ED_FLAGS = ("aht_acute", "aht_any", "aht_alpha1_sensitivity", "anxiolytic", "opioid")
HOME_FLAGS = ("aht_home", "aht_home_sensitivity_only", "anticoagulant", "antiplatelet")
ED_FIELDS = tuple((name, pa.bool_()) for name in ED_FLAGS) + (
    ("aht_ed_category", pa.string()), ("aht_ed_exclusion_reason", pa.string()),
)
HOME_FIELDS = tuple((name, pa.bool_()) for name in HOME_FLAGS) + (
    ("final_categories", pa.string()), ("aht_home_exclusion_reason", pa.string()),
    ("acap_exclusion_reason", pa.string()),
)
ESCALATION_FIELDS = (
    ("hemostasis_level", pa.int32()), ("hemostasis_category", pa.string()),
    ("txa_route", pa.string()), ("ent_referral", pa.bool_()),
    ("reversal_candidate", pa.bool_()), ("reversal_family", pa.string()),
    ("hemostasis_exclusion_reason", pa.string()),
)
LAB_FIELDS = (
    ("lab_measure", pa.string()), ("lab_priority", pa.int32()),
    ("lab_exclusion_reason", pa.string()),
)


def word_pattern(labels):
    return re.compile(r"\b(?:" + "|".join(re.escape(label) for label in labels) + r")\b")


ED_ACUTE = word_pattern(cfg.G0_AHT_ED_ACUTE_DRUGS)
ED_EXTRA = word_pattern(cfg.G0_AHT_ED_ANY_EXTRA_DRUGS)
ED_ORAL_BETA = word_pattern(cfg.G0_AHT_ED_ORAL_BETA_DRUGS)
ED_ALPHA1 = word_pattern(cfg.G0_AHT_ED_ALPHA1_DRUGS)
ED_NIFEDIPINE = word_pattern(cfg.G0_AHT_ED_NIFEDIPINE)
ED_METOPROLOL = word_pattern(cfg.G0_AHT_ED_METOPROLOL)
ED_ORAL_ONLY = word_pattern(cfg.G0_AHT_ED_ORAL_ONLY)
ED_ANXIOLYTICS = word_pattern(cfg.G0_AHT_ED_ANXIOLYTICS)
ED_OPIOIDS = word_pattern(cfg.G0_AHT_ED_OPIOIDS)


def normalized(label):
    return " ".join(label.upper().split()) if label is not None else ""


def first_match(rules, label):
    return next((category for category, pattern in rules.items() if re.search(pattern, label)), None)


def classify_aht_ed(name, order_type=None, mimic=False):
    label = normalized(name)
    reason = first_match(cfg.G0_AHT_ED_EXCLUSIONS, label)
    flags = dict.fromkeys(ED_FLAGS, False)
    category = "excluded"
    if not label:
        reason = "Missing medication label"
    elif not mimic and normalized(order_type) != "MEDICATIONS":
        reason = "Order type is not Medications"
    elif not reason:
        iv_pattern = cfg.G0_AHT_ED_MIMIC_IV_PATTERN if mimic else cfg.G0_AHT_ED_MCMED_IV_PATTERN
        iv = bool(re.search(iv_pattern, label))
        oral_label = bool(re.search(cfg.G0_AHT_ED_ORAL_PATTERN, label))
        extended = bool(re.search(cfg.G0_AHT_ED_EXTENDED_PATTERN, label))
        local = bool(re.search(cfg.G0_AHT_ED_LOCAL_PATTERN, label))
        flags["anxiolytic"] = bool(ED_ANXIOLYTICS.search(label))
        flags["opioid"] = bool(ED_OPIOIDS.search(label))
        if ED_METOPROLOL.search(label):
            metoprolol_iv = bool(re.search(cfg.G0_AHT_ED_MIMIC_METOPROLOL_IV_PATTERN, label)) if mimic else iv
            if metoprolol_iv and not oral_label:
                category = "acute"
            elif not iv and (oral_label or extended or re.search(cfg.G0_AHT_ED_METOPROLOL_ORAL_FORM_PATTERN, label)):
                category = "any_only"
            else:
                reason = "Metoprolol route not resolved by the approved formulation labels"
        elif ED_NIFEDIPINE.search(label):
            if local or iv:
                reason = "Nifedipine preparation is not an eligible systemic release form"
            else:
                category = "any_only" if extended else "acute"
        elif ED_ORAL_ONLY.search(label):
            if not iv and not local and (oral_label or extended):
                category = "any_only"
            else:
                reason = "Only oral minoxidil, diltiazem or verapamil qualifies"
        elif ED_ACUTE.search(label):
            category = "acute"
        elif ED_EXTRA.search(label) or (ED_ORAL_BETA.search(label) and not iv and not local):
            category = "any_only"
        elif ED_ALPHA1.search(label):
            category = "alpha1_sensitivity_only"
        elif ED_ORAL_BETA.search(label):
            reason = "Only oral beta blocker preparations qualify for class B"
        else:
            reason = "Outside G0 antihypertensive inclusion list"
    flags["aht_acute"] = category == "acute"
    flags["aht_any"] = category in ("acute", "any_only")
    flags["aht_alpha1_sensitivity"] = category == "alpha1_sensitivity_only"
    return flags | {"aht_ed_category": category, "aht_ed_exclusion_reason": reason}


def classify_home_acap(med_class, subclass=None, mimic=False):
    label = normalized(med_class)
    sublabel = normalized(subclass)
    combined = " ".join((label, sublabel))
    home_rules = cfg.G0_MIMIC_HOME_AHT_EXCLUSIONS if mimic else cfg.G0_MCMED_HOME_AHT_EXCLUSIONS
    acap_rules = cfg.G0_MIMIC_ACAP_EXCLUSIONS if mimic else cfg.G0_MCMED_ACAP_EXCLUSIONS
    home_reason = first_match(home_rules, combined.strip())
    acap_reason = first_match(acap_rules, combined.strip())
    if mimic:
        home = label in cfg.G0_MIMIC_HOME_AHT_CLASSES
        sensitivity = bool(re.search(cfg.G0_MIMIC_HOME_AHT_SENSITIVITY_PATTERN, label))
        anticoagulant = label in cfg.G0_MIMIC_ACAP_ANTICOAGULANT_CLASSES
        antiplatelet = label in cfg.G0_MIMIC_ACAP_ANTIPLATELET_CLASSES
    else:
        home = label in cfg.G0_MCMED_HOME_AHT_CLASSES or bool(re.search(cfg.G0_MCMED_HOME_AHT_COMBINATION_PATTERN, label))
        sensitivity = bool(re.search(cfg.G0_MCMED_HOME_AHT_SENSITIVITY_PATTERN, label)) or bool(
            re.search(cfg.G0_MCMED_HOME_AHT_SENSITIVITY_PATTERN, sublabel)
        )
        anticoagulant = label in cfg.G0_MCMED_ACAP_ANTICOAGULANT_CLASSES or bool(
            re.search(cfg.G0_MCMED_ACAP_THROMBIN_PATTERN, label)
        )
        antiplatelet = label in cfg.G0_MCMED_ACAP_ANTIPLATELET_CLASSES
        if label == cfg.G0_MCMED_ACAP_SALICYLATE_CLASS:
            antiplatelet = sublabel in cfg.G0_MCMED_ACAP_SALICYLATE_SUBCLASSES
            if not antiplatelet:
                acap_reason = "Only verified single-agent or buffered salicylate subclasses qualify"
    flags = {
        "aht_home": home and not sensitivity and not home_reason,
        "aht_home_sensitivity_only": sensitivity and not home_reason,
        "anticoagulant": anticoagulant and not acap_reason,
        "antiplatelet": antiplatelet and not acap_reason,
    }
    missing = "Missing medication class label"
    if not flags["aht_home"] and not flags["aht_home_sensitivity_only"]:
        home_reason = home_reason or ("Outside G0 home antihypertensive inclusion list" if label else missing)
    if not flags["anticoagulant"] and not flags["antiplatelet"]:
        acap_reason = acap_reason or ("Outside G0 antithrombotic inclusion list" if label else missing)
    return {name: bool(value) for name, value in flags.items()} | {
        "final_categories": "|".join(name for name, value in flags.items() if value) or "excluded",
        "aht_home_exclusion_reason": home_reason,
        "acap_exclusion_reason": acap_reason,
    }


def classify_escalation(order_type, name):
    label = normalized(name)
    category = None
    level = None
    route = None
    referral = False
    reason = first_match(cfg.G0_ESC_EXCLUSIONS, label)
    reversal_family = first_match(cfg.G0_ESC_REVERSAL_PATTERNS, label)
    if not label:
        reason = "Missing procedure label"
    elif normalized(order_type) in cfg.G0_ESC_EXCLUDED_ORDER_TYPES:
        reason = "Laboratory or imaging order excluded"
    elif not reason:
        if label in cfg.G0_ESC_LOCAL_ORDERS:
            level, category = 2, "epistaxis_local_care"
        elif re.search(cfg.G0_ESC_SILVER_PATTERN, label):
            level, category = 2, "silver_nitrate"
        elif re.search(cfg.G0_ESC_TXA_PATTERN, label):
            level, category = 2, "tranexamic_acid"
            if re.search(cfg.G0_ESC_TXA_LOCAL_PATTERN, label):
                route = "topical_or_nebulized"
            elif re.search(cfg.G0_ESC_TXA_SYSTEMIC_PATTERN, label):
                route = "intravenous_or_oral"
            else:
                route = "unknown"
        elif label in cfg.G0_ESC_ENT_ORDERS:
            level, category = 3, "ent_consult"
        elif label in cfg.G0_ESC_ANESTHESIA_ORDERS:
            level, category = 3, "anesthesia_consult"
        elif normalized(order_type) in cfg.G0_ESC_TRANSFUSION_ORDER_TYPES and re.search(cfg.G0_ESC_TRANSFUSION_PATTERN, label):
            level, category = 3, "transfusion"
        elif label in cfg.G0_ESC_ENT_REFERRALS:
            category, referral = "outpatient_ent_referral", True
        elif reversal_family:
            category = "reversal_candidate_for_g1"
            reason = "Candidate awaiting G1 confirmation"
        else:
            reason = "Outside G0 hemostasis inclusion list"
    return {
        "hemostasis_level": level, "hemostasis_category": category, "txa_route": route,
        "ent_referral": referral, "reversal_candidate": reversal_family is not None,
        "reversal_family": reversal_family, "hemostasis_exclusion_reason": reason,
    }


def classify_g1_reversal(order_type, name, reversal_family):
    """Apply only the reviewer-approved reversal additions to the G0 order dictionary."""
    kind, label = normalized(order_type), normalized(name)
    if label == cfg.G1_ESC_VITAMIN_K_LAB_LABEL or kind in cfg.G0_ESC_EXCLUDED_ORDER_TYPES:
        reason = "Laboratory or imaging order excluded"
    elif reversal_family not in cfg.G1_ESC_REVERSAL_FAMILIES:
        reason = "Reversal family not confirmed by G1"
    elif reversal_family in ("phytonadione", "desmopressin") and kind != cfg.G1_ESC_REVERSAL_MEDICATION_ORDER:
        reason = "G1 requires a Medications order for this reversal drug"
    elif reversal_family == "desmopressin" and (
        not re.search(cfg.G1_ESC_DESMOPRESSIN_ROUTE_PATTERN, label)
        or re.search(cfg.G1_ESC_DESMOPRESSIN_EXCLUDED_ROUTE_PATTERN, label)
    ):
        reason = "G1 requires IVPB/INJ desmopressin and excludes oral and nasal forms"
    else:
        reason = None
    return {"reversal_confirmed": reason is None, "g1_hemostasis_level": 3 if reason is None else None,
            "g1_reversal_exclusion_reason": reason}


def classify_lab(name):
    label = normalized(name)
    measure, priority = cfg.G0_LAB_LABELS.get(label, (None, None))
    reason = None
    if measure is None:
        if not label:
            reason = "Missing lab component label"
        elif re.search(cfg.G0_LAB_EXCLUSION_PATTERN, label):
            reason = "Excluded lab component"
        else:
            reason = "Outside G0 lab inclusion list"
    return {"lab_measure": measure, "lab_priority": priority, "lab_exclusion_reason": reason}


def classify_disposition(name):
    admission_labels = tuple(normalized(value) for value in cfg.G0_ESC_ADMISSION_DISPOSITIONS)
    admission = normalized(name) in admission_labels
    return {
        "hemostasis_level": 3 if admission else None,
        "hemostasis_category": "admission" if admission else None,
        "hemostasis_exclusion_reason": None if admission else "Outside G0 admission categories",
    }


def write_table(output, fields, records):
    temporary = output.with_suffix(".csv.tmp")
    try:
        with temporary.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fields)
            writer.writeheader()
            writer.writerows(records)
        temporary.replace(output)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise RuntimeError(f"Could not write {output.name}") from None


def map_labels(con, source, keys, table_name, fields, classify, output_dir):
    key_columns = ", ".join(keys)
    try:
        labels = con.execute(f"""
            SELECT {key_columns}, count(*) AS record_count
            FROM {source}
            GROUP BY {key_columns}
            ORDER BY {key_columns}
        """).fetchall()
        records = [dict(zip(keys, row[:-1])) | {"record_count": row[-1]} | classify(*row[:-1]) for row in labels]
        schema = pa.schema([(key, pa.string()) for key in keys] + [("record_count", pa.int64()), *fields])
        lookup = pa.Table.from_pylist(records, schema=schema)
        con.register("p1_dictionary_mapping", lookup)
        try:
            con.execute(f"CREATE OR REPLACE TEMP TABLE {table_name} AS SELECT * FROM p1_dictionary_mapping")
        finally:
            con.unregister("p1_dictionary_mapping")
        if output_dir is not None:
            write_table(output_dir / f"{table_name.removeprefix('p1_dict_')}_mapping.csv", schema.names, records)
    except Exception:
        raise RuntimeError(f"Could not generate G0 mapping for {source}") from None


def reversal_frequencies(con, output_dir):
    try:
        result = con.execute("""
            WITH visits AS (
                SELECT CSN, TRY_CAST(Arrival_time AS TIMESTAMP) AS arrival_ts,
                    TRY_CAST(Departure_time AS TIMESTAMP) AS departure_ts
                FROM mcmed_visits
            ), candidate AS (
                SELECT o.CSN, d.Order_type, d.Procedure_name, d.reversal_family,
                    TRY_CAST(o.Order_time AS TIMESTAMP) AS order_ts
                FROM mcmed_orders o
                JOIN p1_dict_mcmed_orders d
                    ON o.Order_type IS NOT DISTINCT FROM d.Order_type AND o.Procedure_name = d.Procedure_name
                WHERE d.reversal_candidate
            )
            SELECT c.Order_type, c.Procedure_name, c.reversal_family,
                count(*) AS order_count, count(DISTINCT c.CSN) AS visit_count,
                count(DISTINCT c.CSN) FILTER (WHERE c.order_ts BETWEEN v.arrival_ts AND v.departure_ts)
                    AS ed_order_visit_count
            FROM candidate c
            LEFT JOIN visits v ON c.CSN = v.CSN
            GROUP BY c.Order_type, c.Procedure_name, c.reversal_family
            ORDER BY order_count DESC, c.Order_type, c.Procedure_name
        """)
        fields = [field[0] for field in result.description]
        rows = [dict(zip(fields, values)) for values in result.fetchall()]
        write_table(output_dir / "mcmed_reversal_candidates.csv", fields, rows)
    except Exception:
        raise RuntimeError("Could not generate anticoagulant reversal candidate frequencies") from None


def run(con, output_dir=DICTIONARY_DIR):
    if output_dir is not None:
        output_dir.mkdir(parents=True, exist_ok=True)
    map_labels(
        con, "mcmed_orders", ("Order_type", "Procedure_name"), "p1_dict_mcmed_orders", ED_FIELDS + ESCALATION_FIELDS,
        lambda kind, name: classify_aht_ed(name, kind) | classify_escalation(kind, name), output_dir,
    )
    map_labels(
        con, "mimic_pyxis", ("name",), "p1_dict_mimic_pyxis", ED_FIELDS,
        lambda name: classify_aht_ed(name, mimic=True), output_dir,
    )
    map_labels(
        con, "mcmed_meds", ("Med_class", "Med_subclass"), "p1_dict_mcmed_meds", HOME_FIELDS,
        classify_home_acap, output_dir,
    )
    map_labels(
        con, "mimic_medrecon", ("etcdescription",), "p1_dict_mimic_medrecon", HOME_FIELDS,
        lambda name: classify_home_acap(name, mimic=True), output_dir,
    )
    map_labels(
        con, "mcmed_labs", ("Component_name",), "p1_dict_mcmed_labs", LAB_FIELDS,
        classify_lab, output_dir,
    )
    map_labels(
        con, "mcmed_visits", ("ED_dispo",), "p1_dict_mcmed_disposition", (
            ("hemostasis_level", pa.int32()), ("hemostasis_category", pa.string()),
            ("hemostasis_exclusion_reason", pa.string()),
        ), classify_disposition, output_dir,
    )
    if output_dir is not None:
        reversal_frequencies(con, output_dir)
