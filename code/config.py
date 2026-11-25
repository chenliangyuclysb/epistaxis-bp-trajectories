"""Parameters used by the P0 environment and data audit."""

RANDOM_SEED = 20261125

# S3.1: candidate diagnoses and chief complaints.
EPIS_CC_PATTERN = r"(?i)EPISTAXIS|NOSE ?BLEED"
MCMED_DX_PATTERN = r"(?i)^R04\.0"
MIMIC_DX10_PATTERN = r"(?i)^R040"
MIMIC_DX9_CODE = "7847"
MIMIC_PRIMARY_DX_SEQ = 1
MIN_ADULT_AGE = 18

# S3.2: trauma and postoperative exclusions for the provisional index sample.
TRAUMA_DX_PATTERN = r"(?i)^(?:S02\.?2|S01\.?2|S00\.?3)"
POST_OP_CC_PATTERN = r"(?i)POST-OP|POSTOPERATIVE"
POST_OP_DX_PATTERN = r"(?i)postsurgical|postoperative|postprocedural"

# S3.4: blood-pressure eligibility.
C1_MIN_BP_READINGS = 2
C1_MIN_SPAN_MINUTES = 30

# S4.3: exact chart duplicates and near-duplicate audit windows.
BP_CHART_TRIAGE_DUP_MINUTES = 10
BP_CROSS_SOURCE_WINDOWS_MINUTES = (1, 3, 5)
BP_CROSS_SOURCE_MAX_SBP_DIFF = 3

# S4.4 and S4.5: reading validity and observation limit.
BP_SBP_MIN, BP_SBP_MAX = 50, 300
BP_DBP_MIN, BP_DBP_MAX = 20, 200
BP_MIN_PULSE_PRESSURE = 10
ED_MAX_HOURS = 24

# S5.1, S5.4, S5.5: historical reconnaissance checks from docs/03.
RECON_SBP_CUTOFFS = (140, 180)
RECON_MCMED_MIN_DENSE_BP_READINGS = 4
RECON_MIMIC_BP_READING_COUNTS = (3, 6)

# S6.1: coded hypertension in MC-MED PMH.
MCMED_HTN_ICD10_PATTERN = r"(?i)^I1[0-6](?:\.?[0-9A-Z]*)$"
MCMED_HTN_ICD9_PATTERN = r"^40[1-5](?:\.?[0-9]*)$"

# S7.1 and S7.2: unstratified historical return-count checks.
RECON_RETURN_DAYS = (7, 30)

# S7.2: broad hemostatic escalation label candidates for G0 review.
HEMOSTASIS_ESCALATION_PATTERN = (
    r"(?i)\bENT\b|OTOLARYN|EPISTAXIS|PACK|RHINO|TAMPON|BALLOON|CAUTER|"
    r"SILVER|TRANEXAMIC|TRANSFUS|RED BLOOD|PLATELET|PLASMA|EMBOLI|"
    r"ANGIOGRAM|ANGIOGRAPHY|SURGERY|ANESTHES"
)

# S7.2 and S8.1: named-drug and procedure checks used by docs/03.
RECON_DRUG_PATTERNS = {
    "labetalol": r"(?i)\b(?:labetalol|trandate|normodyne)\b",
    "hydralazine": r"(?i)\b(?:hydralazine|apresoline)\b",
    "metoprolol": r"(?i)\b(?:metoprolol|lopressor|toprol)\b",
    "amlodipine": r"(?i)\b(?:amlodipine|norvasc)\b",
    "clonidine": r"(?i)\b(?:clonidine|catapres)\b",
    "nicardipine": r"(?i)\b(?:nicardipine|cardene)\b",
    "enalaprilat": r"(?i)\b(?:enalaprilat|vasotec iv)\b",
    "esmolol": r"(?i)\b(?:esmolol|brevibloc)\b",
    "oxymetazoline": r"(?i)\b(?:oxymetazoline|afrin)\b",
    "tranexamic_acid": r"(?i)\b(?:tranexamic(?: acid)?|cyklokapron|lysteda)\b",
    "silver_nitrate": r"(?i)\bsilver nitrate\b",
    "mcmed_epistaxis_supply": r"(?i)\b(?:EPISTAXIS MANAGEMENT|BEDSIDE SUPPLIES FOR EPISTAXIS)\b",
    "mcmed_ent": r"(?i)\b(?:CONSULT(?:ATION)?(?: TO)? ENT|OTOLARYNGOLOGY)\b",
    "mcmed_transfusion": r"(?i)\b(?:TRANSFUS(?:E|ED|ION)|RED BLOOD CELLS)\b",
    "acei": (
        r"(?i)\b(?:benazepril|captopril|enalapril(?:at)?|fosinopril|lisinopril|"
        r"moexipril|perindopril|quinapril|ramipril|trandolapril|vasotec|prinivil|zestril)\b"
    ),
    "arb": (
        r"(?i)\b(?:azilsartan|candesartan|eprosartan|irbesartan|losartan|"
        r"olmesartan|telmisartan|valsartan|cozaar|diovan)\b"
    ),
}

# S8.1 and S6.1: broad list for G0 review, including drugs with other indications.
AHT_DRUG_PATTERN = (
    r"(?i)\b(?:labetalol|hydralazine|nicardipine|enalapril(?:at)?|metoprolol|esmolol|"
    r"amlodipine|benazepril|captopril|fosinopril|lisinopril|moexipril|perindopril|"
    r"quinapril|ramipril|trandolapril|azilsartan|candesartan|eprosartan|irbesartan|"
    r"losartan|olmesartan|telmisartan|valsartan|aliskiren|clonidine|guanfacine|"
    r"methyldopa|clevidipine|diltiazem|felodipine|isradipine|nimodipine|nifedipine|"
    r"nisoldipine|verapamil|atenolol|bisoprolol|carvedilol|nebivolol|nadolol|"
    r"propranolol|acebutolol|betaxolol|pindolol|sotalol|timolol|minoxidil|doxazosin|"
    r"prazosin|terazosin|nitroprusside|chlorthalidone|chlorothiazide|hydrochlorothiazide|"
    r"indapamide|metolazone|triamterene|amiloride|spironolactone|eplerenone|furosemide|"
    r"bumetanide|torsemide|ethacrynic acid|norvasc|lopressor|toprol|trandate|normodyne|"
    r"apresoline|cardene|vasotec|brevibloc|catapres|prinivil|zestril|cozaar|diovan|coreg|lasix)\b"
)
AHT_IV_ROUTE_PATTERN = r"(?i)\b(?:iv|intravenous|ivpb|ivp|iv push)\b"
AHT_ORAL_ROUTE_PATTERN = r"(?i)\b(?:po|oral|by mouth|tabs?|tablets?|caps?|capsules?|chewable)\b"

# S6.1 and S9: medication reconciliation terminology candidates.
MIMIC_MEDRECON_ACAP_AHT_PATTERN = (
    r"(?i)antihypertens|(?:angiotensin|renin)\b|\bACE\b|calcium.channel|"
    r"beta.{0,30}(?:block|antagon)|alpha.{0,30}(?:block|antagon)|"
    r"(?:loop|thiazide|potassium.sparing|aldosterone).{0,30}(?:diuretic|antagon)|"
    r"diuretic|vasodilator|adrenergic.{0,30}(?:inhibitor|block|antagon)|"
    r"anticoagula|antiplatelet|platelet.{0,30}(?:inhib|antagon)|heparin|"
    r"factor.{0,30}(?:inhibitor|antagon)|vitamin k antagonist|thrombin inhibitor|"
    r"central alpha[- ]2 receptor agonist|antiarrhythmic\s*-\s*class\s+(?:ii|iv)\b|"
    r"\bsalicylate analgesic"
)

# S9: hemoglobin, platelet count, INR and prothrombin time candidates.
LAB_COMPONENT_PATTERN = r"(?i)\b(?:hgb|hb|hemoglobin|haemoglobin|plt|platelets?|inr|pt|prothrombin)\b"

# P0 G0: minimum parse rate among populated timestamp fields.
TIMESTAMP_MIN_PARSE_PERCENT = 99.9

# SAP v0.3 S5.6, S10: MC-MED spline models and numerical optimization.
P2_MODEL_HOURS = (0, 6)
P2_MODEL_SPLINE_DF = 3
P2_MODEL_CURVE_GRID_POINTS = 61
P2_MODEL_CONFIDENCE_LEVEL = 0.95
P2_LCMM_CLASS_COUNTS = (1, 2, 3, 4)
P2_LCMM_MIN_CLASS_PROPORTION = 0.05
P2_LCMM_MIN_MEAN_POSTERIOR = 0.70
P2_LCMM_GRID_STARTS = 6
P2_LCMM_GRID_ITERATIONS = 12
P2_LCMM_FINAL_ITERATIONS = 100
P2_R_TIMEOUT_SECONDS = 600

# S8.1 and S11.15: finalized G0 D-AHT-ED generics and product-name aliases.
G0_AHT_ED_ACUTE_DRUGS = (
    "LABETALOL", "TRANDATE", "NORMODYNE", "HYDRALAZINE", "APRESOLINE",
    "NICARDIPINE", "CARDENE", "CLEVIDIPINE", "CLEVIPREX", "ENALAPRILAT",
    "ESMOLOL", "BREVIBLOC", "NITROPRUSSIDE", "NIPRIDE", "NITROPRESS", "CLONIDINE",
    "CATAPRES", "CAPTOPRIL", "CAPOTEN",
)
G0_AHT_ED_ANY_EXTRA_DRUGS = (
    "AMLODIPINE", "NORVASC", "FELODIPINE", "PLENDIL", "LISINOPRIL", "PRINIVIL",
    "ZESTRIL", "ENALAPRIL", "VASOTEC", "BENAZEPRIL", "LOTENSIN", "RAMIPRIL",
    "ALTACE", "QUINAPRIL", "ACCUPRIL", "FOSINOPRIL", "MONOPRIL", "LOSARTAN",
    "COZAAR", "VALSARTAN", "DIOVAN", "IRBESARTAN", "AVAPRO", "OLMESARTAN",
    "BENICAR", "TELMISARTAN", "MICARDIS", "CANDESARTAN", "ATACAND",
    "HYDROCHLOROTHIAZIDE", "HCTZ", "CHLORTHALIDONE", "INDAPAMIDE",
    "HYZAAR", "EXFORGE", "LOTREL", "AVALIDE",
)
G0_AHT_ED_ORAL_BETA_DRUGS = (
    "ATENOLOL", "TENORMIN", "CARVEDILOL", "COREG", "BISOPROLOL", "ZEBETA",
    "NEBIVOLOL", "BYSTOLIC", "NADOLOL", "CORGARD",
)
G0_AHT_ED_ALPHA1_DRUGS = ("DOXAZOSIN", "PRAZOSIN", "TERAZOSIN", "CARDURA", "MINIPRESS", "HYTRIN")
G0_AHT_ED_NIFEDIPINE = ("NIFEDIPINE", "PROCARDIA", "ADALAT")
G0_AHT_ED_METOPROLOL = ("METOPROLOL", "LOPRESSOR", "TOPROL")
G0_AHT_ED_ORAL_ONLY = ("MINOXIDIL", "LONITEN", "DILTIAZEM", "CARDIZEM", "VERAPAMIL", "CALAN", "ISOPTIN")
G0_AHT_ED_ANXIOLYTICS = ("LORAZEPAM", "ATIVAN", "DIAZEPAM", "VALIUM", "ALPRAZOLAM", "XANAX", "MIDAZOLAM", "VERSED")
G0_AHT_ED_OPIOIDS = (
    "MORPHINE", "HYDROMORPHONE", "DILAUDID", "FENTANYL", "SUFENTANIL",
    "REMIFENTANIL", "ALFENTANIL", "OXYCODONE", "OXYCONTIN", "PERCOCET",
    "HYDROCODONE", "NORCO", "VICODIN", "CODEINE", "TRAMADOL", "ULTRAM",
    "MEPERIDINE", "DEMEROL", "METHADONE", "BUPRENORPHINE", "BUTRANS",
    "OXYMORPHONE", "TAPENTADOL", "NUCYNTA", "BUTORPHANOL", "NALBUPHINE",
    "PENTAZOCINE", "LEVORPHANOL",
)
G0_AHT_ED_EXCLUSIONS = {
    "Loop or other excluded diuretic": r"\b(?:FUROSEMIDE|LASIX|BUMETANIDE|TORSEMIDE|ETHACRYNIC ACID|METOLAZONE|CHLOROTHIAZIDE)\b",
    "Aldosterone antagonist excluded from ED exposure": r"\b(?:SPIRONOLACTONE|EPLERENONE)\b",
    "Excluded alternative-indication drug": r"\b(?:NITROGLYCERIN|NIMODIPINE|SOTALOL|PROPRANOLOL|GUANFACINE|PHENTOLAMINE)\b",
    "Sacubitril-valsartan excluded": r"\b(?:SACUBITRIL|ENTRESTO)\b",
    "Ophthalmic preparation excluded": r"\b(?:OPHTH(?:ALMIC)?|EYE DROPS?|TIMOLOL|BETAXOLOL)\b",
    "Investigational preparation excluded": r"^INV",
}
G0_AHT_ED_EXTENDED_PATTERN = r"\b(?:ER|CR|XL|XR|SR|LA|EXTENDED|SUSTAINED|CONTROLLED|TR24|TB24|TBSR|CP12|CP24|CS24|C24P|CDCR)\b"
G0_AHT_ED_MCMED_IV_PATTERN = r"\b(?:IV|IVPB|IVP|INTRAVENOUS|INJ|INJECTION)\b"
G0_AHT_ED_MIMIC_IV_PATTERN = r"\b(?:IV|IVPB|IVP|INTRAVENOUS|INJ|INJECTION|VIAL|FOR DILUTION)\b|MG\s*/\s*(?:[0-9]+(?:\.[0-9]+)?\s*)?ML\b"
G0_AHT_ED_MIMIC_METOPROLOL_IV_PATTERN = r"\bVIAL\b|MG\s*/\s*(?:[0-9]+(?:\.[0-9]+)?\s*)?ML\b"
G0_AHT_ED_ORAL_PATTERN = r"\b(?:PO|ORAL|BY MOUTH|TABS?|TABLETS?|CAPS?|CAPSULES?|CHEWABLE|SPLIT TABLET)\b"
G0_AHT_ED_LOCAL_PATTERN = r"\b(?:TOPICAL|TP|OINTMENT|CREAM|GEL|OPHTH(?:ALMIC)?|EYE DROPS?)\b"
G0_AHT_ED_METOPROLOL_ORAL_FORM_PATTERN = r"\b(?:SUCCINATE|TOPROL)\b"

# S6.1 and S9: MC-MED medication classes for G0 D-AHT-home and D-ACAP.
G0_MCMED_HOME_AHT_CLASSES = (
    "ANGIOTENSIN RECEPTOR ANTAG.-THIAZIDE DIURETIC COMB",
    "CALCIUM CHANNEL BLOCKING AGENTS", "BETA-ADRENERGIC BLOCKING AGENTS",
    "ALPHA/BETA-ADRENERGIC BLOCKING AGENTS", "THIAZIDE AND RELATED DIURETICS",
    "POTASSIUM SPARING DIURETICS IN COMBINATION", "RENIN INHIBITORS",
)
G0_MCMED_HOME_AHT_COMBINATION_PATTERN = (
    r"^ANTIHYPERTENSIVES,|\bRENIN INHIBITORS?\b|"
    r"^(?:ACE INHIBITORS?|ANGIOTENSIN RECEPTOR ANTAG\.?|ANGIOTENSIN II RECEPTOR BLOCKERS?)\b.*\bCOMB(?:INATION)?S?\b"
)
G0_MCMED_HOME_AHT_SENSITIVITY_PATTERN = (
    r"^(?:POTASSIUM SPARING DIURETICS|ALPHA-ADRENERGIC BLOCKING AGENTS)$|"
    r"^BPH\b.*\bALPHA[- ]?1\b|\bPROSTATIC HYPERTROPHY AGENT\b.*\bALPHA[- ]1\b"
)
G0_MCMED_HOME_AHT_EXCLUSIONS = {
    "Inhaled adrenergic or ocular class excluded": r"\b(?:INHALED|MIOTICS|OPHTH(?:ALMIC)?|INTRAOCULAR PRESSURE REDUCERS)\b",
    "Loop diuretic excluded": r"\b(?:LOOP DIURETICS|DIURETIC - LOOP)\b",
    "Coronary vasodilator excluded": r"\b(?:VASODILATORS,\s*CORONARY|CORONARY VASODILATORS)\b",
    "Antidepressant, replacement, beta-3 or vasopressor class excluded": r"\b(?:ALPHA-2 RECEPTOR ANTAGONIST ANTIDEPRESSANTS|CALCIUM REPLACEMENT|BETA\s*-\s*3|VASOPRESSORS?)\b",
}
G0_MCMED_ACAP_ANTIPLATELET_CLASSES = ("PLATELET AGGREGATION INHIBITORS",)
G0_MCMED_ACAP_SALICYLATE_CLASS = "ANALGESIC/ANTIPYRETICS, SALICYLATES"
G0_MCMED_ACAP_SALICYLATE_SUBCLASSES = ("SALICYLATE ANALGESICS", "SALICYLATE ANALGESICS, BUFFERED")
G0_MCMED_ACAP_ANTICOAGULANT_CLASSES = (
    "ANTICOAGULANTS,COUMARIN TYPE", "ANTICOAGULANTS, COUMARIN TYPE",
    "DIRECT FACTOR XA INHIBITORS", "HEPARIN AND RELATED PREPARATIONS",
)
G0_MCMED_ACAP_THROMBIN_PATTERN = r"\b(?:DIRECT )?THROMBIN INHIBITORS?\b"
G0_MCMED_ACAP_EXCLUSIONS = {
    "Platelet-reducing agent excluded": r"\bPLATELET REDUCING AGENTS\b",
    "NSAID class excluded": r"\b(?:NSAIDS?|NON-STEROIDAL ANTI-INFLAMMATORY)\b",
}

# S6.1 and S9: MIMIC medication reconciliation classes for G0 D-AHT-home and D-ACAP.
G0_MIMIC_HOME_AHT_CLASSES = (
    "ACE INHIBITORS", "ACE INHIBITOR AND DIURETIC COMBINATIONS",
    "ACE INHIBITOR AND CALCIUM CHANNEL BLOCKER COMBINATIONS",
    "ANGIOTENSIN II RECEPTOR BLOCKERS (ARBS)",
    "ANGIOTENSIN II RECEPTOR BLOCKER (ARB)-DIURETIC COMBINATIONS",
    "ANGIOTENSIN II RECEPTOR BLOCKER (ARB)-CALCIUM CHANNEL BLOCKER COMB.",
    "ANGIOTENSIN II RECEPTOR BLOCKER (ARB)-CALCIUM CHANNEL BLOCKER-DIURETIC",
    "CALCIUM CHANNEL BLOCKERS - DIHYDROPYRIDINES",
    "CALCIUM CHANNEL BLOCKERS - BENZOTHIAZEPINES",
    "CALCIUM CHANNEL BLOCKERS - PHENYLAKYLAMINES",
    "ANTIHYPERLIPIDEMIC HMG COA REDUCT INHIB AND CALCIUM CHANNEL BLOCKER",
    "BETA BLOCKERS CARDIAC SELECTIVE", "BETA BLOCKERS NON-CARDIAC SELECTIVE",
    "BETA BLOCKERS CARDIAC SELECTIVE, INTRINSIC SYMPATHOMIMETIC ACTIVITY",
    "BETA BLOCKERS NON-CARDIAC SELECT., INTRINSIC SYMPATHOMIMETIC ACTIVITY",
    "CARDIAC SELECTIVE BETA BLOCKER-THIAZIDE DIURETIC AND RELATED COMB.",
    "NON-CARDIAC SELECTIVE BETA BLOCKER-THIAZIDE DIURETIC AND RELATED COMB.",
    "ALPHA-BETA BLOCKERS", "DIURETIC - THIAZIDES AND RELATED",
    "DIURETIC - POTASSIUM SPARING-THIAZIDE AND RELATED COMBINATIONS",
    "CENTRAL ALPHA-2 RECEPTOR AGONISTS",
    "CENTRAL ALPHA-2 AGONISTS-THIAZIDE DIURETIC AND RELATED COMB.",
    "DIRECT ACTING VASODILATORS", "RENIN INHIBITOR, DIRECT",
    "RENIN INHIBITOR, DIRECT AND DIURETIC COMBINATIONS",
)
G0_MIMIC_HOME_AHT_SENSITIVITY_PATTERN = (
    r"\bALDOSTERONE(?: RECEPTOR)? ANTAGONISTS?\b|"
    r"\b(?:PROSTATIC HYPERTROPHY AGENT|PERIPHERAL)\b.*\bALPHA[- ]1\b"
)
G0_MIMIC_HOME_AHT_EXCLUSIONS = {
    "Ophthalmic category excluded": r"\bOPHTH(?:ALMIC)?\b",
    "Loop or other excluded diuretic": r"\b(?:DIURETIC - LOOP|DIURETIC - POTASSIUM SPARING|CARBONIC ANHYDRASE INHIBITORS?)$",
    "ARNi category excluded": r"\b(?:ARNI|NEPRILYSIN)\b",
    "Other indication excluded": r"\b(?:PULMONARY ANTIHYPERTENSIVE|NITRATES?|NASSA|ANTIARRHYTHMIC\s*-\s*CLASS\s+(?:II|IV)|CEREBROVASCULAR SPECIFIC|VASOPRESSIN|VASOPRESSOR HORMONES)\b",
}
G0_MIMIC_ACAP_ANTIPLATELET_CLASSES = (
    "PLATELET AGGREGATION INHIBITORS - SALICYLATES", "SALICYLATE ANALGESICS",
    "SALICYLATE ANALGESICS, BUFFERED",
    "PLATELET AGGREGATION INHIBITORS - THIENOPYRIDINE AGENTS",
    "PLATELET AGGREGATION INHIB - CYCLOPENTYL-TRIAZOLO-PYRIMIDINES (CPTPS)",
    "PLATELET AGGREGATION INHIBITORS - PHOSPHODIESTERASE III INHIBITORS",
    "PLATELET AGGREGATION INHIBITOR COMBINATIONS",
    "PLATELET AGGREGATION INHIB-PDESTERASE AND ADENOSINE DEAMINASE INHIBITR",
)
G0_MIMIC_ACAP_ANTICOAGULANT_CLASSES = (
    "ANTICOAGULANTS - COUMARIN", "DIRECT FACTOR XA INHIBITORS",
    "THROMBIN INHIBITOR - SELECTIVE DIRECT AND REVERSIBLE",
    "LOW MOLECULAR WEIGHT HEPARINS", "INDIRECT FACTOR XA INHIBITORS", "HEPARINS",
)
G0_MIMIC_ACAP_EXCLUSIONS = {
    "Heparin flush excluded": r"\bHEPARIN FLUSH FORMULATIONS\b",
    "Other salicylate combination excluded": r"\bSALICYLATE\b.*\bCOMBINATIONS?\b",
    "Platelet-reducing quinazoline category excluded": r"\bQUINAZOLINE AGENTS\b",
}

# S7.2b: finalized G0 D-ESC and reversal candidates for G1 verification.
G0_ESC_LOCAL_ORDERS = ("EPISTAXIS MANAGEMENT", "BEDSIDE SUPPLIES FOR EPISTAXIS")
G0_ESC_ENT_ORDERS = ("CONSULT TO ENT",)
G0_ESC_ANESTHESIA_ORDERS = ("CONSULT TO ANESTHESIOLOGY", "CONSULT TO PERIOP ANESTHESIA")
G0_ESC_ENT_REFERRALS = ("REFERRAL TO ENT/OTOLARYNGOLOGY",)
G0_ESC_TRANSFUSION_ORDER_TYPES = ("BLOOD BANK", "NURSING TRANSFUSION", "TRANSFUSION COMMUNICATION")
G0_ESC_ADMISSION_DISPOSITIONS = ("Inpatient", "ICU")
G0_ESC_SILVER_PATTERN = r"\bSILVER NITRATE\b"
G0_ESC_TXA_PATTERN = r"\bTRANEXAMIC ACID\b"
G0_ESC_TXA_LOCAL_PATTERN = r"\b(?:TOPICAL|TP|NEB|NEBULI[ZS]ED|NEBULI[ZS]ATION|INHALATION)\b"
G0_ESC_TXA_SYSTEMIC_PATTERN = r"\b(?:IV|IVPB|IVP|INTRAVENOUS|PO|ORAL|TABS?|TABLETS?|INFUSION|MINIBAG)\b"
G0_ESC_TRANSFUSION_PATTERN = r"\b(?:RBCS?|RED BLOOD CELLS?|PLATELETS?|PLT|PLASMA|CRYO(?:PRECIPITATE)?|MTP|MASSIVE TRANSFUSION)\b"
G0_ESC_EXCLUDED_ORDER_TYPES = ("LAB", "IMAGING")
G0_ESC_EXCLUSIONS = {
    "Diagnostic angiography excluded": r"\b(?:CT|MR|CTA|MRA)\b.*\bANGIO(?:GRAM|GRAPHY)?\b|\b(?:CTA|MRA)\b",
    "Supportive or risk-assessment order excluded": r"\b(?:ICE PACK|VTE RISK)\b",
    "Regional anesthesia consult excluded": r"\bCONSULT TO REGIONAL ANESTHESIA\b",
}
G0_ESC_REVERSAL_PATTERNS = {
    "phytonadione": r"\b(?:PHYTONADIONE|VITAMIN K\s*1?|AQUAMEPHYTON)\b",
    "prothrombin_complex": r"\b(?:PROTHROMBIN COMPLEX|KCENTRA)\b",
    "idarucizumab": r"\b(?:IDARUCIZUMAB|PRAXBIND)\b",
    "andexanet": r"\b(?:ANDEXANET|ANDEXXA)\b",
    "desmopressin": r"\b(?:DESMOPRESSIN|DDAVP|STIMATE)\b",
}

# SAP v0.3 S7.2b: G1 approval for medication reversal labels and routes.
G1_ESC_REVERSAL_FAMILIES = ("phytonadione", "idarucizumab", "prothrombin_complex", "desmopressin")
G1_ESC_REVERSAL_MEDICATION_ORDER = "MEDICATIONS"
G1_ESC_VITAMIN_K_LAB_LABEL = "VITAMIN K1, SERUM"
G1_ESC_DESMOPRESSIN_ROUTE_PATTERN = r"\b(?:IVPB|INJ)\b"
G1_ESC_DESMOPRESSIN_EXCLUDED_ROUTE_PATTERN = r"\b(?:PO|ORAL|NASAL|INTRANASAL|SPRAY|TAB(?:S|LETS?)?)\b"

# S9: finalized G0 D-LAB, including point-of-care fallback labels.
G0_LAB_LABELS = {
    "HEMOGLOBIN (HGB)": ("hemoglobin", 1),
    "HGB": ("hemoglobin", 1),
    "POC:HGB (CALC), ISTAT": ("hemoglobin", 2),
    "PLATELET COUNT (PLT)": ("platelet_count", 1),
    "INR": ("inr", 1),
    "POC:INR, ISTAT": ("inr", 2),
}
G0_LAB_EXCLUSION_PATTERN = r"\b(?:MCHC?|CO-HEMOGLOBIN|HBA1C|A1C|GLYCATED HEMOGLOBIN|PLATELET MORPHOLOGY|PT|PROTHROMBIN TIME)\b"

# SAP v0.2 S3.2, S3.4: trauma and disposition definitions.
TRAUMA_DX9_PATTERN = r"^(?:802\.?[01]|873\.?[23])"
INCOMPLETE_DISPOSITIONS = (
    "LWBS", "ELOPED", "AMA", "LEFT WITHOUT BEING SEEN", "LEFT AGAINST MEDICAL ADVICE",
)
DIRECT_HOME_MCMED = "Discharge"
DIRECT_HOME_MIMIC = "HOME"

# SAP v0.2 S3.3 and S7.1–S7.2: historical and return windows.
PRIOR_EPISODES_LOOKBACK_DAYS = 365
RETURN_7D_HOURS = 168
RETURN_30D_HOURS = 720

# SAP v0.2 S4.3: finalized main deduplication window.
G0_BP_DUP_MINUTES = 3

# SAP v0.3 S5.1 and S5.4: BP thresholds and minimum count for ARV.
BP_ELEVATED_SBP, BP_ELEVATED_DBP = 140, 90
BP_SEVERE_SBP, BP_SEVERE_DBP = 180, 110
BP_MIN_ARV_READINGS = 4

# SAP v0.3 S4.3 and S11-14: widened window for MC-MED dense metrics only.
MC_DENSE_ALT_DUP_MINUTES = 5

# SAP v0.3 S5.5: Monte Carlo sparse sampling and descriptive percentiles.
SPARSE_BP_REPETITIONS = 500
SPARSE_MIN_POST_TRIAGE_READINGS = 1
SPARSE_BP_INTERVAL_QUANTILES = (0.025, 0.975)

# SAP v0.3 S12: Figure 2 pooled reading-level LOESS span.
FIG2_LOESS_FRACTION = 0.35

# SAP S12 and G5 6.1b/6.5: final Figure 2 display range and figure export resolution.
P6_FIG2_MAIN_MAX_HOURS = 12
P6_FIGURE_DPI = 600

# SAP v0.2 S7.2b: inpatient disposition and landmark timing.
ESC_ADMISSION_DISPOSITIONS = ("Inpatient", "ICU")
ESC_LANDMARK_MINUTES = 60

# SAP v0.2 S9 and G0 dictionary supplement: baseline diagnostic codes.
CKD_DX10_PATTERN = r"(?i)^N18"
CKD_DX9_PATTERN = r"^585"
LIVER_DX10_PATTERN = r"(?i)^K7[0-7]"
LIVER_DX9_PATTERN = r"^57[12]"
HHT_DX10_PATTERN = r"(?i)^I78\.?0$"
HHT_DX9_PATTERN = r"^448\.?0$"

# SAP v0.2 S9: triage acuity and conservative triage pain ascertainment.
TRIAGE_ACUITY_GROUPS = (("1–2", (1, 2)), ("3", (3,)), ("4–5", (4, 5)))
MCMED_TRIAGE_ACUITY_CODES = {
    "1-RESUSCITATION": 1,
    "2-EMERGENT": 2,
    "3-URGENT": 3,
    "4-SEMI-URGENT": 4,
    "5-NON-URGENT": 5,
}
TRIAGE_PAIN_MIN, TRIAGE_PAIN_MAX = 0, 10
TRIAGE_PAIN_MEASURES = ("PAIN", "PAIN SCORE")
TRIAGE_PAIN_LOOKUP_MINUTES = 10

# SAP v0.2 S10: two-sided 95% Wilson interval for descriptive proportions.
WILSON_95_Z = 1.959963984540054

# SAP v0.3 S10 and S12: recovery summary times, risk table, and figure sampling.
AIM1_KM_LANDMARKS_MINUTES = (60, 120)
AIM1_KM_TIMEPOINTS_MINUTES = (0, 60, 120, 180, 240, 300, 360)
AIM1_KM_PLOT_STEP_MINUTES = 1

# Frozen SAP v1.0 S6.2: patient-level Aim 2 follow-up windows and BP evidence.
P3_MCMED_FOLLOWUP_FIRST_DAY = 1
P3_MCMED_FOLLOWUP_STRICT_FIRST_DAY = 7
P3_MIMIC_FOLLOWUP_FIRST_DAY = 7
P3_FOLLOWUP_LAST_DAY = 365
P3_OMR_RESULT_NAMES = ("Blood Pressure", "Blood Pressure Sitting")
P3_OMR_BP_PATTERN = r"\s*([+-]?[0-9]+(?:[.][0-9]+)?)\s*/\s*([+-]?[0-9]+(?:[.][0-9]+)?)\s*"
P3_OMR_MIN_READINGS = 1
P3_OMR_STRICT_MIN_READINGS = 2
P3_OMR_STRICT_SBP, P3_OMR_STRICT_DBP = 130, 80

# Frozen SAP v1.0 S10: assess adjustment covariate missingness before fitting.
P3_COVARIATE_MI_TRIGGER_PERCENT = 5.0
P3_COVARIATE_IMPUTATIONS = 20

# Frozen SAP v1.0 S10 and G3 4.2: late-SBP exposure and display.
P4_SBP_UNIT = 10
P4_SBP_KNOT_QUANTILES = (0.10, 0.50, 0.90)
P4_SBP_CURVE_QUANTILES = (0.025, 0.975)
P4_SBP_CURVE_POINTS = 101
P4_COVARIATE_MI_TRIGGER_PERCENT = 5.0

# G3 4.4 and 4.5: burden bins and predeclared Aim 3 sensitivity rules.
P4_BURDEN_TERTILE_QUANTILES = (1 / 3, 2 / 3)
P4_S11_HIGH_SBP, P4_S11_HIGH_DBP = 160, 100
P4_PREDEPARTURE_MINUTES = 30
P4_S11_EXTRA_DISPOSITIONS = (
    "LWBS", "LEFT WITHOUT BEING SEEN", "AMA", "LEFT AGAINST MEDICAL ADVICE", "OBSERVATION",
)

# SAP S8.2-2 and G4 5.2: risk-set matching, response windows, and uncertainty.
P5_MATCH_ANCHOR_LOOKBACK_MINUTES = 60
P5_MATCH_FOLLOWUP_MINUTES = 150
P5_MATCH_TIME_CALIPER_MINUTES = 30
P5_MATCH_SBP_CALIPER_MMHG = 10
P5_MATCH_MAX_CONTROLS = 3
P5_MATCH_RESPONSE_WINDOWS = ((60, 30, 90), (120, 90, 150))
P5_MATCH_TARGETS_MINUTES = tuple(window[0] for window in P5_MATCH_RESPONSE_WINDOWS)
P5_MATCH_BOOTSTRAPS = 2000
P5_MATCH_CI_QUANTILES = (0.025, 0.975)

# SAP v0.2 S7.2c: conservative inpatient procedure title candidates.
INPATIENT_ENT_HEMOSTASIS_TITLE_PATTERN = (
    r"(?i)^(?:control of epistaxis|control (?:of )?bleeding (?:from|in) "
    r"(?:the )?nasal mucosa and soft tissue)\b"
)

# Frozen SAP S8.2, S11-15 and G4 5.1/5.3: treatment summaries and rare outcomes.
P5_TIMING_QUANTILES = (0.25, 0.5, 0.75)
# S8.2.1 and AGENTS.md 2A: avoid reproducing a singleton measurement as a quantile.
P5_TIMING_MIN_SUMMARY_N = 2
P5_TREATMENT_SBP_UNIT = 10
P5_RARE_TREATED_EVENT_COUNT = 10

# G7 8.1, post hoc relative to SAP S7/S9/S10: bleeding-severity sensitivity.
P8_PYXIS_CANDIDATE_PATTERN = (
    r"(?i)oxymetazoline|afrin|tranexamic|cyklokapron|silver nitrate|cocaine|"
    r"thrombin|gelfoam|surgicel|floseal|surgiflo|phenylephrine"
)
P8_PYXIS_TXA_PATTERN = r"(?i)tranexamic|cyklokapron"
P8_PYXIS_PHENYLEPHRINE_PATTERN = r"(?i)phenylephrine"
P8_PYXIS_LOCAL_FORM_PATTERN = r"(?i)\b(?:NASAL|INTRANASAL|NASL|NOSE|TOPICAL|TOP)\b"
P8_PYXIS_NONLOCAL_FORM_PATTERN = (
    r"(?i)\b(?:I\.?V\.?|IVPB|IVP|INTRAVENOUS|INFUSION|INJ|INJECTION|ORAL|PO|TABLETS?|TABS?|"
    r"OPHTH|OPHTHALMIC|EYE|INTRAOCULAR)\b"
)
P8_PYXIS_SYSTEMIC_THROMBIN_PATTERN = r"(?i)\b(?:ANTI[ -]?THROMBIN|PROTHROMBIN)\b"
P8_HEMOGLOBIN_FLUID = "blood"
P8_HEMOGLOBIN_LABEL = "hemoglobin"
P8_HEMOGLOBIN_UNITS_G_DL = ("G/DL", "GM/DL")
P8_HEMOGLOBIN_UNITS_G_L = ("G/L", "GM/L")
P8_HEMOGLOBIN_GL_PER_GDL = 10
P8_HEMOGLOBIN_MODEL_INCREMENT_G_DL = 1
P8_MIMIC_HEMOGLOBIN_MIN_AVAILABLE_PERCENT = 50
P8_MIN_STRATUM_EVENTS = 10
P8_MIN_EVENTS_PER_PARAMETER = 10

# G7 8.2: double-column submission figure width, raster resolution, and text floor.
P8_FIGURE_WIDTH_MM = 178
P8_FIGURE_DPI = 600
P8_FIGURE_MIN_TEXT_PT = 7
