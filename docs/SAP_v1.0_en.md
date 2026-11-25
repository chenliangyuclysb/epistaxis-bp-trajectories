# 01 Study Protocol and Statistical Analysis Plan (SAP), **frozen version v1.0** (2026-11-25, G2)

> English translation of the frozen Chinese original. Freeze date: 2026-11-25 (G2). SHA-256 of the exact original file: `8400c3b4eddc7a7fc460d0f8c04614aab3c663bf6bb403186c056d2a6f8de572`. The Chinese source is provided alongside this translation. Historical project file references are reproduced as written in the original.

Status: **Frozen**. Any subsequent change must be recorded in `review/deviations.md` and reported in the article as a deviation from the prespecified plan. Outcome association analyses may begin only after the freeze (AGENTS.md D). Earlier versions are archived in `review/SAP_v0.x_存档.md`. The revision history appears at the end. The dictionaries are in `docs/04_词典定稿_G0.md`.

## S0 Title (suggested)

- Chinese title, translated here: Elevated blood pressure in emergency department epistaxis: transient stress response or sustained hypertension? A two-center retrospective cohort study of repeated monitor-captured blood pressure trajectories, emergency antihypertensive treatment, and early hemostatic failure.
- English: Transient stress response or sustained hypertension? Serial monitor-captured blood pressure trajectories, emergency antihypertensive treatment, and early hemostatic failure in emergency department epistaxis: a two-center retrospective cohort study
- Avoid “high-frequency” and “continuous.” MC-MED blood pressure (BP) measurements are non-invasive and occur approximately every 15 minutes (docs/03 §2).

## S1 Research Questions

- **Aim 1 (main descriptive and physiological question):** How often is BP elevated on arrival in patients with epistaxis attending the emergency department (ED)? Among patients without antihypertensive treatment, how often and how quickly does BP normalize during the ED stay? Does this course differ by history of hypertension?
- **Aim 2 (anchoring persistent elevation):** Among patients without known hypertension, is persistent BP elevation during the ED stay, compared with transient elevation, associated with subsequent evidence of sustained hypertension?
- **Aim 3 (outcome):** Is the ED BP phenotype associated with early hemostatic failure within 7 days of discharge? **Primary contrast: P versus T** (S5.3).
- **Aim 4 (secondary, descriptive):** Frequency and determinants of ED antihypertensive treatment, BP change relative to untreated comparators matched on time, and exploratory association with hemostatic failure.

## S2 Design and Data

- Two-center retrospective cohort study.
- Center A: Stanford, MC-MED (2020–2022). The dense BP observations are used to construct trajectory phenotypes. (The MC-MED version number was awaiting confirmation from the user.)
- Center B: BIDMC, MIMIC-IV-ED 2.2 (2011–2019) plus MIMIC-IV 3.1 `hosp`. BP measurements are sparse, so replication uses the two-point phenotype. Outpatient BP from `omr` is used for Aim 2.
- Reporting guidelines: STROBE and RECORD.

## S3 Cohort

- **S3.1 Inclusion** (harmonized to mean that epistaxis was the main reason for the visit)
  - MC-MED: `Dx_ICD10` starts with R04.0, or `CC` matches `EPISTAXIS|NOSE ?BLEED` (case-insensitive).
  - MIMIC: `diagnosis` has `seq_num=1` with ICD-10 R040* or ICD-9 7847, or `triage.chiefcomplaint` matches the same regular expression.
  - Age ≥18 years. MIMIC age is `anchor_age + (year of intime − anchor_year)`, from `hosp/patients`. In the G0 audit, 100% of epistaxis stays matched the patients table.
- **S3.2 Primary-analysis exclusions**
  - Trauma: an index diagnosis in S02.2*, S01.2*, or S00.3*, or ICD-9 802.0*, 802.1*, 873.2*, or 873.3*. Use any diagnosis position in MIMIC and the primary diagnosis in MC-MED.
  - Invalid timing: departure is missing or is no later than arrival.
  - Postoperative or postprocedural presentation: `CC` contains POST-OP/POSTOPERATIVE, or the diagnosis name contains postsurgical/postoperative/postprocedural.
  - Incomplete encounters (for outcome analyses): LWBS, ELOPED, or AMA.
- **S3.3 Index visits (G0 revision):** The primary analysis includes all new episodes. A visit is a new episode if its arrival falls outside (0, 168 h] after departure from every previous epistaxis ED visit. A visit within this window is itself an outcome event and cannot be used as a new index. Multiple episodes per patient are handled with patient-clustered robust standard errors (SEs), with adjustment for whether a patient had an epistaxis ED visit in the previous 365 days within the available data. Sensitivity analysis: select only the first visit per patient.
  - Reason for the revision: The G0 reanalysis showed that selecting only the first visit per patient yielded 17 MC-MED and 88 MIMIC 7-day C2 events. Including all new episodes yielded 20 and 104 events, respectively (124 total). The primary model had approximately 10 degrees of freedom (df), giving about 12 events per variable. This revision was based solely on event counts before examining associations.
- **S3.4 Analysis subsets**
  - C1 (BP phenotype cohort): At least 2 BP readings during the ED visit, with an interval of at least 30 min between first and last.
  - C2 (outcome cohort): C1 visits discharged directly home from the ED (MC-MED `ED_dispo=="Discharge"`, MIMIC `disposition=="HOME"`). Observation is excluded from the primary analysis.
  - C3 (persistence anchoring cohort): C1 visits without known prior hypertension (S6.1).

## S4 BP Data Preparation

- **S4.1 Sources:** MC-MED uses triage (`Triage_SBP`/`DBP` in `visits`, timed at `Arrival_time`) plus SBP/DBP in `numerics` (Chart and Monitor). MIMIC uses `triage` (timed at `intime`) plus `vitalsign`.
- **S4.2 SBP–DBP pairing:** Pair SBP and DBP from the same CSN, Source, and Time as one reading.
- **S4.3 Deduplication (finalized at G0):**
  - Rule A: Remove Chart readings recorded within 10 min after arrival if both SBP and DBP exactly equal the triage values (P0: 307/314 readings).
  - Rule B: If Chart and Monitor readings are within ±3 min and |ΔSBP|≤3, retain Monitor.
  - Rule C: If triage and any other source are within ±3 min and |ΔSBP|≤3, retain triage as the arrival anchor and drop the other reading. In MIMIC, triage time is `intime`, and the same rule applies to `vitalsign`.
  - Sensitivity analysis: Widen the window to ±5 min for dense MC-MED measures only (average real variability, burden, and recovery time). C1 counts are insensitive to this window, as verified in P0.
  - Remove readings outside the ED stay (254 in MIMIC).
- **S4.4 Plausibility:** Require SBP 50–300, DBP 20–200, and SBP minus DBP ≥10. Remove readings that fail these criteria and report the proportion removed.
- **S4.5 Time window:** t=0 is arrival. The window is [arrival, departure] (MIMIC: [`intime`, `outtime`]), capped at 24 h.

## S5 BP Exposure and Phenotypes

- **S5.1 Thresholds:** Elevated BP is SBP ≥140 or DBP ≥90. Severe elevation is SBP ≥180 or DBP ≥110. Normal BP is SBP <140 and DBP <90.
- **S5.2:** Arrival BP is the first reading. Final BP is the last reading during the ED stay.
- **S5.3 Primary phenotypes (both centers):**
  - N: Normal at arrival.
  - T (transient): Elevated at arrival, normal at the last reading.
  - P (persistent): Elevated at arrival and still elevated at the last reading.
- **S5.4 Dense measures in MC-MED:**
  - Subdivide P into P-never (never normal) and P-rebound (at least 1 intermediate normal reading followed by elevation at the last reading).
  - Time to normalization: Kaplan–Meier analysis. The event is the first normal reading, with censoring at the last reading.
  - Time-weighted mean SBP by trapezoidal integration.
  - Hypertension burden: Area under the SBP curve above 140 divided by the observation time (mmHg).
  - Average real variability (ARV) of SBP, calculated when at least 4 readings are available.
  - Peak SBP and whether severe elevation occurred.
- **S5.5 Sparse-sampling agreement (redefined at G1):** The primary phenotype uses only the first and last points, so comparing “triage plus last” with the dense classification would necessarily yield perfect agreement. Instead, use 500 Monte Carlo thinning replicates in MC-MED C1 (seed 20261125). For each visit in each replicate, draw k from the empirical MIMIC C1 distribution of counts of post-triage readings, with k≥1 and no greater than the visit's observed count. Randomly sample k post-triage readings without replacement and reclassify using triage plus the last thinned reading. Compare with the full-data primary phenotype. Report mean κ and its 2.5th and 97.5th percentiles, the misclassification matrix, and sensitivity for any normalization during the ED stay. Also report the distribution of the time from the last reading to departure in both centers.
- **S5.6 Data-driven trajectories (MC-MED, optional sensitivity analysis):** Use R `lcmm::hlme` with natural splines over 0–6 h and k=1–4. Select according to BIC, smallest class ≥5%, and mean posterior probability ≥0.7. If R is unavailable, skip the analysis and record this.

## S6 Chronic Hypertension Anchor

- **S6.1 Prior hypertension** (before index)
  - MC-MED: `pmh` contains I10–I16 or 401–405, and the date of `Noted_date` is **before** the index arrival date. Same-day codes are excluded from prior hypertension because they may have been recorded during the index visit. A sensitivity analysis counts them as prior hypertension. Alternatively, an active outpatient antihypertensive (D-AHT-home) qualifies if the date of `Entry_date` precedes index, `Start_date` is missing or no later than index, and `End_date` is missing or no earlier than index. `Active` is mostly blank and is unused.
  - MC-MED outpatient medication records are incomplete: 121 of 378 epistaxis patients have no `meds` record. Medication status therefore means documented use. Retain a separate indicator for availability of an outpatient medication list and restrict to patients with a list in a sensitivity analysis.
  - MIMIC: An antihypertensive appears on index `medrecon`, or a hypertension diagnosis is recorded in a prior inpatient stay (`hosp diagnoses_icd`, `admittime < intime`) or a prior ED stay.
- **S6.2 Evidence of sustained hypertension after index** (C3 only)
  - MIMIC: Use `omr` readings with `result_name` equal to `Blood Pressure` or `Blood Pressure Sitting`, excluding postural measurements such as Standing or Lying. Parse `result_value` as SBP/DBP and apply S4.4. `chartdate` must fall on days 7–365 after the index departure date. Exclude days during any hospital admission (dates of `admittime` through `dischtime`, inclusive) or any ED visit. Average all readings in the window. Mean SBP ≥140 or mean DBP ≥90 indicates evidence of sustained hypertension. The primary analysis requires at least 1 reading. A sensitivity analysis requires at least 2 readings and applies a threshold of 130/80.
  - MC-MED: Use a newly recorded hypertension code (`Noted_date`) or newly started antihypertensive (`Entry_date`/`Start_date`) on days 1–365 after index. Exclude same-day records. A sensitivity analysis starts on day 7.
  - Analyze individuals with follow-up data. Follow-up availability for MIMIC requires at least 1 qualifying `omr` reading in the window. In MC-MED, any `pmh` `Noted_date`, `meds` `Entry_date`, or later MC-MED visit on days 1–365 qualifies. P3 reports follow-up proportions and compares patients with and without follow-up in each center, including phenotype distributions.
  - Unit of analysis (v1.0): Aim 2 uses patients, selecting the first C3-eligible visit per patient to avoid overlapping follow-up windows.

## S7 Outcomes

- **S7.1 Primary outcome, early hemostatic failure:** The next epistaxis ED visit, defined as in S3.1 without the trauma exclusion, arrives in (0, 168 h] after index departure. An overlapping stay whose arrival is at or before index departure is regarded as an error and is excluded as a return. The analysis population is C2.
- **S7.2 Secondary outcomes:**
  - (a) Return for epistaxis within 30 days.
  - (b) Escalation of care during the MC-MED index visit is D-ESC level 3: ENT consultation, transfusion, anesthesia consultation, anticoagulant reversal, or Inpatient/ICU admission with primary diagnosis R04.0. Use a landmark approach: exposure is early BP in the first 0–60 min and outcome is a level 3 event after 60 min. Exclude ED stays ≤60 min and visits with a level 3 event during 0–60 min from the denominator. The sensitivity analysis allows admission for any primary diagnosis.
  - (c) MIMIC nasal hemostatic ENT procedures after hospital admission (`procedures_icd`), exploratory.

## S8 ED Antihypertensive Treatment

- **S8.1 Exposure (G0 finalization, see docs/04 D-AHT-ED):** Antihypertensives given during the ED stay. In MC-MED, `orders` Medications require `First_admin_time` ∈ [arrival, departure]. In MIMIC, require `pyxis` `charttime` ∈ [`intime`, `outtime`]. The main exposure is class A acute antihypertensives, classified by drug regardless of route because MIMIC Pyxis names often omit route. The secondary exposure is class B, any antihypertensive. Anxiolytics and analgesics are separately recorded.
- **S8.2 Analyses:**
  1. Treatment proportion and BP at the first dose.
  2. BP response: Risk-set matching. Match each treated visit to up to 3 untreated visits at the same center with similar SBP (±10 mmHg) at similar time after arrival (±30 min). Compare SBP changes at 60 and 120 min using mixed models.
  3. Association with the primary outcome: If the treated group has <10 events, report descriptively only (counts and crude RR with exact CI). If the treated group has ≥10 events, use overlap weighting to estimate RR and label it exploratory.

## S9 Covariates

- Age, sex, race/ethnicity (for description), arrival mode, triage acuity (grouped as 1–2/3/4–5), anticoagulants and antiplatelets (recorded separately, primary model uses any antithrombotic, D-ACAP), prior hypertension, CKD, liver disease, HHT (coding in docs/04), epistaxis ED visit in the past 365 days, index pain and heart rate, and index level 2 hemostatic treatment (D-ESC).
- MC-MED laboratory values: First ED hemoglobin, platelet count, and INR.
- Comorbidity timing (G1): CKD and liver disease may come from prior or index documentation because circularity does not arise. In MIMIC, index ED secondary diagnoses can be used. In MC-MED, `Noted_date` may be on or before the index date. HHT has no time restriction. Only prior hypertension strictly requires documentation before index (S6.1).
- Differences in record completeness (G1): Prior hypertension was documented in 38.2% of MC-MED and 69.5% of MIMIC encounters. Any antithrombotic use was documented in 15.0% and 52.2%, respectively. Reviewer reanalysis confirmed the MIMIC values: Index medication reconciliation listed antihypertensives in 58%, prior hospital diagnosis in 45%, and prior ED diagnosis alone identified only 1.2%. The lower MC-MED values mainly reflect incomplete records (121/378 lacked medication lists, and over-the-counter aspirin may be missing). Thus, all analyses involving prior hypertension or antithrombotic use must be reported by center. Retain a center term in pooled models and perform S11-16 and S11-17.
- Pooled primary model (approximately 10 df): Phenotype (2), age (linear plus quadratic or a 3-knot spline, 2 df), sex, any antithrombotic, prior hypertension, triage acuity (2), and center. The MC-MED single-center model adjusts only for age and sex (approximately 20 events).

## S10 Statistical Analysis

- Describe results by center. Use Wilson 95% CIs for proportions.
- **Aim 1 (expanded at G1):**
  - For C1 in both centers, report proportions elevated on arrival, severely elevated on arrival, and normal at last reading, along with N/T/P distributions. MC-MED also reports P-never and P-rebound.
  - Cross-tabulate **phenotype by prior hypertension** by center. For visits with elevation on arrival, estimate the proportion with prior hypertension and the prevalence ratio for P relative to T using modified Poisson regression with clustered SEs. This compares exposures and involves no outcome, so it is permitted at P2.
  - “Spontaneous” normalization: Kaplan–Meier curves among arrivals with elevated BP stratified by prior hypertension and center. The event is the first normal reading. Censor at the last reading and **at the first class A acute antihypertensive dose**. MIMIC readings are sparse, so its Kaplan–Meier analysis is for approximate reference.
  - Linear mixed model (MC-MED, 0–6 h): `SBP ~ ns(time, df=3) × prior hypertension + age + sex`, with a random intercept and random linear slope, nested by visit. Patient-level clustering is an optional sensitivity analysis. Only readings before the first class A dose are used.
  - MIMIC: Paired change in SBP from arrival to last reading, stratified by prior hypertension.
  - Sensitivity analyses: Restrict to visits with medication lists in both centers and exclude visits with ED class A treatment.
- **Aim 2 (v1.0):** In C3 patients with no known hypertension and with follow-up, use one visit per patient. Estimate P vs T RRs, with T as the reference, using modified Poisson regression and robust SEs, using sustained hypertension evidence (S6.2) as the outcome and adjusting for age, sex, and arrival SBP. Describe the incidence in N and report PPV/NPV for P as a marker of subsequent hypertension. Estimate by center first. MIMIC provides the primary evidence. MC-MED is supportive because incomplete documentation may overestimate newly identified disease (S9). Pool with a center term if pooling is feasible. Sensitivities: ≥2 readings with the threshold of 130/80, a day 7 start for MC-MED, and IPW for selection into observed follow-up (model includes age, sex, acuity, arrival SBP, phenotype, and center).
- **Aim 3 (v1.0):** In C2 new episodes, estimate RRs with modified Poisson regression and a log link, using patient-clustered robust SEs (GEE with independence working correlation). The phenotypes are N/T/P and **P vs T is the primary contrast**. Pooled covariates are in S9 (approximately 10 df). The MC-MED single-center model adjusts only for age and sex. The secondary exposure is last SBP per 10 mmHg, using restricted cubic splines with 3 knots. Exploratory MC-MED analyses compare P-rebound with P-never and examine SBP burden, mainly descriptively. Test the phenotype-by-center interaction and supplement it with fixed-effect meta-analysis.
  - Interpretation: The G2 thinning simulation showed that sparse sampling misclassifies about 26% of T as P (κ=0.86). Discuss possible attenuation toward the null for the MIMIC P vs T estimate.
- **Multiplicity:** The sole primary contrast is P vs T for 7-day hemostatic failure. Label other analyses secondary and apply no adjustment for multiplicity.
- **Missing data:** If covariate missingness exceeds 5%, use multiple imputation (m=20). Compare characteristics of those included and excluded because they fail C1 criteria and conduct an IPW sensitivity analysis.
- **Event counts (G0 reanalysis of new episodes):** C2 includes 222 MC-MED visits with 20 events and 901 MIMIC visits with 104 events, a total of 124 events. Statistical power for Aim 3 comes mainly from MIMIC. MC-MED contributes detailed trajectories, and its single-center estimates will necessarily be imprecise. Wording should reflect this.

## S11 Sensitivity Analyses

1. Change the threshold to 160/100.
2. Retain only the first visit per patient (S3.3).
3. Use primary diagnosis only (Dx only for MC-MED and sequence-1 diagnosis only for MIMIC).
4. Exclude patients using anticoagulants or antiplatelets.
5. Exclude HHT.
6. Exclude visits with ED antihypertensive treatment.
7. Include LWBS, AMA, and Observation.
8. Use the 30-day outcome.
9. Replace last BP with the last reading recorded at least 30 min before departure.
10. LCMM trajectories (S5.6).
11. E-value.
12. Restrict MC-MED to patients with outpatient medication lists (S6.1).
13. Include same-day PMH codes as prior hypertension (S6.1).
14. Use a ±5 min deduplication window for dense MC-MED measures (S4.3).
15. Redefine exposure as any antihypertensive (class B), and then include alpha-1 blockers (S8).
16. Restrict both centers to visits with medication lists (S9 record completeness).
17. Restrict MIMIC prior hypertension to medication reconciliation plus previous hospital diagnoses, removing prior ED diagnoses as a source.

## S12 Figures and Tables

- Figure 1: Cohort flow for both centers.
- Figure 2: MC-MED SBP trajectories, with individual curves and LOESS by phenotype. The main x-axis is cropped to 0–12 h because the tail is sparse.
- Figure 3: Kaplan–Meier curves for normalization by prior hypertension.
- Figure 4: Forest plot of center-specific and pooled estimates.
- Figure 5: SBP changes after antihypertensive treatment and in matched untreated comparators.
- Table 1: Baseline characteristics. Table 2: Phenotypes and dense measures. Table 3: Aim 2. Table 4: Aim 3. Supplementary tables: Dictionaries, sensitivity analyses, and sparse-sampling agreement.

## S13 Versioning and Registration

- From v0.1 (this draft) to v1.0, frozen after G2.
- Suggested registering v1.0 in OSF before P4, while blinded to outcome associations.

## Revision History

- **v0.2 (2026-11-25, G0):** S3.1 added the MIMIC age calculation. S3.2 added ICD-9 trauma codes and exclusions for abnormal timing. **S3.3 changed primary inclusion to all new episodes** (rationale in that section). S4.3 finalized deduplication rules. S6.1 clarified same-day coding, active medication status, and incomplete medication list coverage. S6.2 moved the start to day 1. S7.1 set the return window to (0, 168 h]. S8.1 stratified treatment definitions by drug. S9 added covariates and revised model degrees of freedom. S10 updated event counts. S11 added analyses 12–15. All revisions occurred before computing any exposure–outcome association.
- **v0.3 (2026-11-25, G1):** S5.5 changed to Monte Carlo thinning. S7.2b clarified the level 3 escalation risk set. S9 clarified comorbidity timing and differences in record completeness. S10 Aim 1 clarified phenotype × prior hypertension, Kaplan–Meier censoring at acute antihypertensive treatment, and the predose reading restriction for the mixed model. S11 added analyses 16–17. No exposure–outcome association had yet been computed.
- **v1.0 (2026-11-25, G2, frozen):** S6.2 specified `omr` extraction, the definition of follow-up availability, and patient-level analysis for Aim 2. S10 specified models for Aims 2 and 3. The Figure 2 x-axis range was specified. When the plan was frozen, no exposure–outcome associations had yet been analyzed (blinding audits in P0–P2 had passed).
