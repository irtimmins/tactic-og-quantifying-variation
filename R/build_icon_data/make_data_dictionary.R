# =============================================================================
# Build the OG cohort data dictionary -> Excel
# -----------------------------------------------------------------------------
# Self-contained. The reviewed dictionary content is the dict tibble below, one
# row per variable in the cohort's column order. The script fills in the live
# type and completeness, assigns a tier, and writes a four-tab workbook:
# It produces ONE sheet. A short guidance preamble sits ABOVE the table (cohort,
# shading key, curative subsetting, core waiting times, method notes, pathways)
# in light blue; below it, one row per variable, shaded by tier:
#   core      (soft green)  essential for the epidemiological analysis
#   secondary (soft amber)  useful, but mostly feeds a core variable or is
#                           modality-specific
#   ignore    (soft grey)   raw inputs / helpers, not for analysis
#
# To change wording or tiers, edit the dict tibble / core_vars below.
# Needs: dplyr, openxlsx.
# =============================================================================

library(dplyr)
library(openxlsx)

source("R/build_icon_data/01_define_parameters.R")
og <- readRDS(f_cohort_cwt)

# -----------------------------------------------------------------------------
# Reviewed content: variable, definition, values, source, recommendation.
# Rows are in the cohort's column order. type and pct_complete are filled from
# the live data below so they never go stale.
# -----------------------------------------------------------------------------
dict <- tibble::tribble(
  ~variable, ~definition, ~values, ~source, ~recommendation,
  "pseudo_patientid", "Pseudonymised patient identifier. One row per patient in this cohort; the primary key for joining to any patient-level extract. A person can in rare cases hold more than one patient ID in the source.", "40695 distinct values (id / free-text)", "NCRAS Cancer registry (PATIENTID, pseudonymised)", "essential - primary key",
  "pseudo_tumourid", "Pseudonymised tumour identifier. A patient can have more than one tumour, so this is unique per tumour, not per patient; here it is effectively one-to-one with the patient because the cohort is one tumour per person.", "40695 distinct values (id / free-text)", "NCRAS Cancer registry (TUMOURID, pseudonymised)", "may be useful only if linking tumour-level data; unlikely to be needed",
  "diagmdy", "Date of diagnosis. The registry 'best' diagnosis date (DIAGNOSISDATEBEST), derived under UKIACR/ENCR rules: the date is taken from the highest-priority available source in a defined hierarchy (e.g. clinical, then pathological, etc.); where two candidate dates remain ambiguous it is set to the mid-point of diagnosisdate1 and diagnosisdate2. See DIAGNOSISDATEFLAG for how precisely the date is specified / whether it was imputed. It is the index date from which all waiting times and windows are measured.", "range 2015-01-02 to 2022-12-31", "NCRAS Cancer registry (DIAGNOSISDATEBEST; UKIACR/ENCR rules)", "essential - the index date",
  "ydiag", "Year of diagnosis. Calendar year extracted from diagmdy, for cohort selection and year-on-year breakdowns.", "range 2015 to 2022 (median 2018)", "Derived from diagmdy", "essential - year of diagnosis for trends and cohort selection",
  "cancer", "Cancer site in words: oesophagus or stomach. A plain-language label derived from the ICD-10 site code.", "oesophagus | stomach", "Derived from sitestr (NCRAS)", "useful - tumour_site_grp is the analysis split",
  "sitestr", "ICD-10 site code. Four-character code of the tumour's origin (C15x oesophagus, C16x stomach); the basis of the oesophageal/gastric split.", "C150 | C151 | C152 | C153 | C154 | C155 | C158 | C159 | C160 | C161 | C162 | C163 | C164 | C165 | C166 | C168 | C169", "NCRAS Cancer registry (SITE_ICD10_O2)", "useful - tumour_site_grp / cancer_subtype are the analysis splits",
  "typestr", "Morphology description in words. Free-text histology label (e.g. adenocarcinoma, squamous cell carcinoma).", "33 distinct values (id / free-text)", "NCRAS Cancer registry (HISTOLOGY_CODED)", "may be useful to look at, but cancer_subtype is the analysis variable",
  "basisofdiagnosis", "How the diagnosis was established. Registry code for the strongest evidence behind the diagnosis (e.g. 7 = histology of primary, 0 = death certificate only).", "1 | 2 | 4 | 5 | 6 | 7", "NCRAS Cancer registry (BASISOFDIAGNOSIS)", "may be useful to look at, but unlikely to be needed - mostly a quality check",
  "grade", "Tumour differentiation grade. How abnormal the cells look (G1 well differentiated to G4 undifferentiated; GX not assessable).", "G1 | G2 | G3 | G4 | GH | GI | GL | GX", "NCRAS Cancer registry (GRADE)", "useful - a tumour characteristic",
  "behav", "ICD-O behaviour code. Whether the tumour is malignant, in situ, etc.; here effectively constant (malignant) given the cohort definition.", "range 3 to 5 (median 3)", "NCRAS Cancer registry (BEHAVIOUR_ICD10_O2)", "ignore - near-constant in this cohort",
  "stage_best", "Registry best overall stage, raw. The registry's chosen stage including sub-stages (1A, 2B, etc.) and the staging system behind it; cleaned into stage_clean for analysis.", "1 | 1A | 1A1 | 1B | 2 | 2A | 2B | 2E | 3 | 3A | 3B | 3C | 3S", "NCRAS Cancer registry (STAGE_BEST)", "ignore for analysis - use stage_clean",
  "stage_best_system", "Staging system behind stage_best. Which edition/scheme was used (UICC 7/8, AJCC 7, etc.); needed to interpret the raw stage codes.", "| AJCC 7 | UICC 5 | UICC 6 | UICC 7 | UICC 8 | Unknown", "NCRAS Cancer registry (STAGE_BEST_SYSTEM)", "may be useful to look at, but unlikely to be needed - reference only",
  "t_best", "Best T category. Tumour local size/extent (UICC). One of the three components behind overall stage.", "| 0 | 1 | 1a | 1b | 1c | 2 | 2a | 2A | 2b | 2s | 3 | 3a | 3b | 3c | 4 | 4a | 4b", "NCRAS Cancer registry (T_BEST)", "may be useful to look at, but unlikely to be needed - a component of stage",
  "n_best", "Best N category. Regional lymph-node involvement (UICC).", "| 0 | 1 | 1a | 1c | 1mi | 2 | 2a | 2b | 2c | 3 | 3a | 3b | X", "NCRAS Cancer registry (N_BEST)", "may be useful to look at, but unlikely to be needed - a component of stage",
  "m_best", "Best M category. Distant metastasis (UICC).", "| 0 | 1 | X", "NCRAS Cancer registry (M_BEST)", "may be useful to look at, but unlikely to be needed - a component of stage",
  "t_path", "Pathological T category. T from pathology rather than imaging; often missing for non-surgical patients.", "| 0 | 1 | 1a | 1a2 | 1b | 1c | 2 | 2a | 2A | 2b | 2s | 3 | 3a | 3b | 4 | 4a | 4b | is | X", "NCRAS Cancer registry (T_PATH)", "may be useful to look at, but unlikely to be needed - surgical patients mainly",
  "n_path", "Pathological N category. N from pathology.", "| 0 | 1 | 1a | 1b | 1mi | 2 | 2a | 2b | 2c | 3 | 3a | 3b | X", "NCRAS Cancer registry (N_PATH)", "may be useful to look at, but unlikely to be needed - surgical patients mainly",
  "m_path", "Pathological M category. M from pathology.", "| 0 | 1 | 1a | 1b | 1c | X", "NCRAS Cancer registry (M_PATH)", "may be useful to look at, but unlikely to be needed - surgical patients mainly",
  "sex", "Sex. Self-stated gender at diagnosis (1 = male, 2 = female; 0/9 = not known/specified).", "range 1 to 2 (median 1)", "NCRAS Cancer registry (SEX)", "essential - core demographic",
  "agediag", "Age at diagnosis, in years. Patient age on the diagnosis date.", "range 19.2909 to 104.668 (median 73.12526)", "NCRAS Cancer registry (AGE)", "essential - core demographic",
  "birthmdy", "Date of birth. Registry best date of birth (day set to the 1st where only month/year known).", "range 1911-01-15 to 1999-11-15", "NCRAS Cancer registry (BIRTHDATEBEST)", "may be useful to look at, but agediag is normally what you want",
  "ethnicity_group_broad", "Broad ethnicity group. Grouped from the 16+1 census ethnicity categories (White, Asian, Black, Mixed, Chinese, Other, Not known); the grouping was applied by NDRS to this project's specification.", "Asian | Black | Chinese | Mixed | Not known | Other | White", "NCRAS Cancer registry (ETHNICTY_GROUP_BROAD, derived by NDRS)", "essential - core demographic",
  "lsoa11_code", "Lower-layer Super Output Area (2011) of residence at diagnosis. Small-area geography (~1,500 residents) used to attach deprivation and area measures; not used directly in analysis but could support deriving patient travel times to treatment.", "22207 distinct values (id / free-text)", "NCRAS Cancer registry (LSOA11_CODE, pseudonymised)", "may be useful later - e.g. to derive patient travel times to treatment; otherwise the IMD/area fields are derived from it",
  "NHSE_reversed_imd_quintile_lsoas", "Area-level deprivation quintile (REVERSED direction). The English Index of Multiple Deprivation 2019 - a composite area-level measure combining income, employment, education, health, crime, housing and environment, assigned by the patient's LSOA of residence, NOT an individual income measure. NDRS supplies IMD as 1 = least deprived to 5 = most deprived; THIS variable is reversed, so here 1 = MOST deprived and 5 = LEAST deprived. Check the direction before interpreting.", "1 - most deprived | 2 | 3 | 4 | 5 - least deprived", "NCRAS / English IMD 2019 (imd19_quintile_lsoas, reversed)", "essential - deprivation; note 1 = most deprived here",
  "canalliance_2024_code", "Cancer Alliance code (2024 geography). Code for the regional Cancer Alliance of the patient's area; the organisational geography for regional comparisons.", "E56000005 | ... | E56000035", "Derived (NCRAS geography)", "useful - for regional analysis",
  "canalliance_2024_name", "Cancer Alliance name. The readable name matching canalliance_2024_code.", "Cheshire and Merseyside | ... | West Yorkshire and Harrogate", "Derived (NCRAS geography)", "useful - for regional analysis",
  "diag_trust", "Trust of diagnosis (code). NHS trust that diagnosed the patient; the provider unit for diagnosis-level variation.", "132 distinct values (id / free-text)", "NCRAS Cancer registry (DIAG_TRUST, pseudonymised)", "essential - trust of diagnosis",
  "diag_trust_name", "Trust of diagnosis (name). Readable name for diag_trust.", "132 distinct values (id / free-text)", "NCRAS Cancer registry (DIAG_TRUST_NAME)", "useful - readable trust of diagnosis",
  "first_trust", "Trust of first recorded event (code). Trust for the patient's first event in the pathway, which may differ from the diagnosing trust.", "136 distinct values (id / free-text)", "NCRAS Cancer registry (FIRST_TRUST, pseudonymised)", "may be useful to look at, but diag_trust is normally the reference",
  "first_trust_name", "Trust of first event (name). Readable name for first_trust.", "136 distinct values (id / free-text)", "NCRAS Cancer registry (FIRST_TRUST_NAME)", "may be useful to look at, but diag_trust is normally the reference",
  "first_hosp_date", "Date of first recorded event. Date of the patient's first event in the registry pathway.", "2946 distinct values (id / free-text)", "NCRAS Cancer registry (FIRST_HOSP_DATE)", "may be useful to look at, but unlikely to be needed",
  "diag_hosp", "Diagnosing hospital site. Hospital (site-level) of diagnosis - finer than diag_trust where site-level detail is wanted.", "392 distinct values (id / free-text)", "NCRAS Cancer registry (DIAG_HOSP)", "essential - diagnosing hospital site",
  "route_bjc", "Route to diagnosis, BJC grouping. A regrouping of the route algorithm into the categories used in the Routes to Diagnosis (BJC) work; an INPUT used to build route_combined, not the field to report.", "Emergency presentation | GP referral | Inpatient elective | Other outpatient | Screening | TWW", "NCRAS Routes to Diagnosis", "ignore - input to route_combined; report final_route",
  "final_route", "Published registry route to diagnosis. The finalised NCRAS route with all dataset types accounted for; feeds route_combined, which is the variable to use for analysis here. route_code and route_bjc are earlier inputs.", "| Emergency presentation | GP referral | Inpatient elective | Other outpatient | TWW | Unknown", "NCRAS Routes to Diagnosis (FINAL_ROUTE)", "useful - feeds route_combined (the variable to use)",
  "route_code", "Route to diagnosis, internal algorithm code. The raw code the routes algorithm assigns; an input behind final_route, not for reporting.", "63 distinct values (id / free-text)", "NCRAS Routes to Diagnosis (ROUTE_CODE)", "ignore - input to final_route",
  "sg_flag", "Registry surgery flag (raw). NCRAS's own indicator that the patient had surgery to remove the primary tumour (1 = yes, 0 = none recorded). This is the REGISTRY flag, not this build's surgery ascertainment; the build uses HES (had_surgery / curative_surgery), which is the gold standard.", "| 0 | 1", "NCRAS Cancer registry (SG_FLAG)", "ignore for analysis - use had_surgery / curative_surgery",
  "rt_flag", "Registry radiotherapy flag (raw). NCRAS indicator of any radiotherapy (1/0). Superseded for analysis by the RTDS-based had_rt / had_curative_rt.", "| 0 | 1", "NCRAS Cancer registry (RT_FLAG)", "ignore for analysis - use had_rt / had_curative_rt",
  "ct_flag", "Registry chemotherapy flag (raw). NCRAS indicator of any chemotherapy (1/0). Superseded for analysis by the SACT-based had_sact.", "| 0 | 1", "NCRAS Cancer registry (CT_FLAG)", "ignore for analysis - use had_sact",
  "dead", "Death indicator (raw registry). Registry vital-status flag; the cleaned analysis version is 'died'.", "range 0 to 1 (median 1)", "NCRAS Cancer registry (vital status)", "ignore for analysis - use died",
  "finmdy", "Final follow-up / death date. Date of death or last known follow-up.", "range 2015-01-10 to 2023-12-31", "NCRAS / ONS", "needed only for survival analysis - follow-up/death date",
  "dco", "Death-certificate-only flag. Whether the cancer was known only from the death certificate (Y/N); DCO cases are typically excluded.", "N", "NCRAS Cancer registry (DCO)", "may be useful to look at, but unlikely to be needed - quality/exclusion check",
  "morphology_num", "Numeric ICD-O morphology code. The numeric histology code (8000-9990) behind cancer_subtype.", "range 8010 to 8576 (median 8140)", "NCRAS Cancer registry (MORPH_ICD10_O2)", "ignore for analysis - use cancer_subtype",
  "tumour_site_grp", "Site group for analysis: oesophageal or gastric. The two-way split used throughout, derived from the ICD-10 site.", "gastric | oesophageal", "Derived from sitestr", "essential - the oesophageal/gastric split",
  "cancer_subtype", "Tumour subtype: Oes ACA, Oes SCC, or Gast. Oesophageal adenocarcinoma vs squamous cell carcinoma (which behave differently), with gastric as one group; derived from morphology.", "Gast | Oes ACA | Oes SCC", "Derived from morphology (NCRAS)", "essential - the subtype split",
  "stage_clean", "Cleaned analysis stage, restricted to 1-3. Sub-stages collapsed to whole numbers and the cohort limited to stages 1-3; this IS the analysis-stage variable and part of the cohort definition.", "1 | 2 | 3", "Derived from stage_best", "essential - the analysis stage and part of the cohort definition",
  "final_route_chr", "Cleaned label of final_route. A tidied character version of final_route for display.", "Emergency presentation | GP referral | Inpatient elective | Other outpatient | TWW | Unknown", "Derived from final_route", "useful - same content as final_route",
  "route_bjc_chr", "Cleaned label of route_bjc. Tidied character version of the BJC route input.", "Emergency presentation | GP referral | Inpatient elective | Other outpatient | Screening | TWW", "Derived from route_bjc", "ignore - input only",
  "route_combined", "Route to diagnosis - the variable to use. The combined route assembled for analysis, including an emergency-presentation category; use this for route-to-diagnosis work and to exclude emergency cases. final_route and the raw route fields are inputs behind it.", "Emergency presentation | GP referral | Inpatient elective | Other outpatient | Screening | TWW | Unknown", "Derived from final_route / route inputs", "essential - the route-to-diagnosis variable (use to exclude emergency cases)",
  "emergency_admission", "Diagnosed via an emergency admission (flag). 1 if the route to diagnosis was an emergency presentation (NOGCA Performance Indicator 1). route_combined already carries an emergency-presentation category, so use that to exclude emergency cases rather than this flag.", "range 0 to 1 (median 0)", "Derived from final_route", "useful - but route_combined already flags emergency cases",
  "surv_from_dx_days", "Days survived from diagnosis. Interval from diagnosis to death or censoring; the basis for overall survival.", "range 1 to 3285 (median 537)", "Derived (diagnosis to ONS death)", "needed only for survival analysis - the survival time",
  "died", "Death indicator for analysis. Cleaned 0/1 death flag paired with surv_from_dx_days.", "range 0 to 1 (median 1)", "Derived", "needed only for survival analysis - the event flag (pair with surv_from_dx_days)",
  "ps_num", "WHO performance status at diagnosis (0-4). How well the patient functions day-to-day (0 fully active to 4 completely disabled); in many contexts a key risk-adjustment variable but often incomplete.", "range 0 to 4 (median 1)", "NCRAS COSD (tumour_performancestatus)", "useful - risk-adjustment variable, but note completeness",
  "cnsinvolved", "Clinical nurse specialist involvement (flag). Whether a CNS was recorded as involved (NOGCA PI4); incomplete.", "range 0 to 1 (median 1)", "NCRAS COSD (clinicalnursespecialist)", "may be useful to look at, but unlikely to be needed - low completeness",
  "endoscopy_date", "Date of the diagnostic endoscopy. First diagnostic upper-GI endoscopy near diagnosis, from HES (admitted or outpatient); the clock start behind the endoscopy-based waits.", "range 2014-12-08 to 2022-12-31", "HES APC / OP (OPCS-4 endoscopy codes, NOGCA Appendix 6)", "may be useful to look at, but the derived endoscopy waits are normally what you need",
  "days_endo_to_dx", "Days from endoscopy to diagnosis. Interval between the diagnostic endoscopy and the diagnosis date; useful context for the diagnostic interval but not a core analysis wait.", "range 0 to 30 (median 0)", "Derived", "may be useful to look at, but not expected to be a core analysis variable",
  "emresd_date", "Date of EMR/ESD endotherapy. Endoscopic mucosal/submucosal resection - a curative treatment for early disease; presence of a date means the patient had EMR/ESD in window.", "range 2015-01-05 to 2023-06-24", "HES APC (OPCS-4 EMR/ESD codes, NOGCA Appendix 7)", "essential - a curative treatment date",
  "emresd_provider", "Provider (trust) of the EMR/ESD. Trust that performed the endoscopic resection; used to assign tx_trust for EMR-only patients.", "131 distinct values (id / free-text)", "HES APC (PROCODE3 on the EMR/ESD episode)", "may be useful to look at, but unlikely to be needed - trust attribution only",
  "days_dx_to_emresd", "Days from diagnosis to EMR/ESD. Interval to the endoscopic resection.", "range -29 to 275 (median 46)", "Derived", "useful - interval to EMR/ESD",
  "surgery_date", "Date of major OG resection. First major oesophageal/gastric resection in window, from HES; the surgical clock-stop. HES is the gold-standard source for surgery.", "range 2015-01-01 to 2023-09-07", "HES APC (OPCS-4 resection codes, NOGCA Appendix 8)", "essential - the surgery date",
  "surgery_type", "Resection type. Oesophagectomy, total gastrectomy, or partial gastrectomy; derived from the OPCS-4 code.", "oesophagectomy | partial_gastrectomy | total_gastrectomy", "Derived from OPCS-4 (HES APC)", "useful - resection type",
  "surgery_class", "Resection class. The two-way oesophagectomy vs gastrectomy split.", "gastrectomy | oesophagectomy", "Derived from OPCS-4 (HES APC)", "useful - resection class",
  "opcs_primary", "Primary OPCS-4 resection code. The specific procedure code matched as the resection.", "32 distinct values (id / free-text)", "HES APC (OPERTN)", "may be useful to look at, but unlikely to be needed - reference",
  "PROCODE3", "Provider (trust) of the surgery episode. Three-character trust code of the resecting provider; used for tx_trust on surgical pathways.", "78 distinct values (id / free-text)", "HES APC (PROCODE3)", "may be useful to look at, but unlikely to be needed - trust attribution only",
  "SITETRET", "Site of treatment for the surgery episode. Site-level provider code for the resection.", "105 distinct values (id / free-text)", "HES APC (SITETRET)", "ignore for analysis - PROCODE3 is the trust",
  "days_dx_to_surg", "Days from diagnosis to surgery. Interval to the resection.", "range -16 to 275 (median 149)", "Derived", "useful - interval to surgery",
  "curative_surgery", "Curative-intent resection (flag). TRUE where the resection counts as curative; follows NOGCA by excluding stage-4 partial gastrectomies. Only TRUE values appear because the cohort is stage 1-3.", "TRUE", "Derived (NOGCA rule)", "useful - feeds received_curative_tx_audit (the headline flag)",
  "sact_date", "Date of first systemic anti-cancer therapy. First chemo/immunotherapy in window (from SACT, supplemented by HES delivery codes); the chemo clock-stop.", "range 2015-01-14 to 2023-07-27", "SACT (+ HES APC delivery codes)", "essential - the chemo date",
  "days_dx_to_sact", "Days from diagnosis to SACT. Interval to first systemic therapy.", "range -29 to 275 (median 58)", "Derived", "useful - interval to chemo",
  "chemo_source", "Where the chemo record came from. sact, hes, or both - shows whether the chemo was found in SACT, in HES delivery codes, or both (a data-provenance flag).", "both | hes | sact", "Derived", "may be useful to look at, but unlikely to be needed - provenance, not clinical",
  "hes_chemo_date", "Chemo date from HES delivery codes. Chemo identified from HES OPCS-4 delivery codes; supplements SACT where SACT is missing it.", "range 2015-01-15 to 2023-08-08", "HES APC (OPCS-4 SACT delivery codes, NOGCA Appendix 9)", "may be useful to look at, but unlikely to be needed - supplements sact_date",
  "BENCHMARK_GROUP", "SACT benchmark regimen group (raw input). SACT's mapping of the drug regimen into a high-level benchmark group. A RAW SACT field used to decide whether chemo is palliative (per NOGCA Appendix 10), not an analysis output.", "112 distinct values (id / free-text)", "SACT (Benchmark_Group)", "ignore - raw SACT input used in derivation",
  "benchmark_group_lwr", "Lower-cased BENCHMARK_GROUP. Helper version for matching; not for analysis.", "112 distinct values (id / free-text)", "Derived from BENCHMARK_GROUP", "ignore - helper only",
  "INTENT_OF_TREATMENT_V3", "SACT drug treatment intent (raw input). SACT's coded intent of the regimen (adjuvant / neoadjuvant / curative / palliative / disease-modification), sometimes multiple codes per patient. A RAW SACT field; the build's tx_intent_audit is the analysis intent variable.", "01 | 01|02 | ... | 99", "SACT (Intent_of_Treatment_v3)", "ignore - raw SACT input; use tx_intent_audit",
  "CHEMO_RADIATION", "SACT chemo-radiation indicator (raw input). SACT v2 flag that a regimen was given with radiation; a raw input to the chemoRT logic, not an output.", "n | N | y | Y", "SACT (Chemo_Radiation)", "ignore - raw SACT input",
  "ORGANISATION_CODE_OF_PROVIDER", "Provider of SACT administration. Trust/organisation that gave the systemic therapy.", "170 distinct values (id / free-text)", "SACT (Organisation_Code_of_Provider)", "may be useful to look at, but unlikely to be needed - SACT-provider only",
  "rt_date", "Date of first radiotherapy. First radiotherapy prescription in window (from RTDS); the radiotherapy clock-stop.", "range 2015-01-26 to 2023-09-06", "RTDS (TREATMENTSTARTDATE / APPTDATE)", "essential - the radiotherapy date",
  "rt_curative", "Curative-intent radiotherapy (flag). TRUE where the dose/fractionation matches a curative schedule (per the NOGCA dose list); distinguishes radical from palliative RT.", "FALSE | TRUE", "Derived from RTDS dose/fractions", "useful - feeds received_curative_tx_audit (the headline flag)",
  "rt_dose", "Radiotherapy total dose (Gray). Total prescribed dose; an input to rt_curative.", "range 2.75 to 70 (median 45)", "RTDS (RTPRESCRIBEDDOSE)", "may be useful to look at, but unlikely to be needed - input to rt_curative",
  "rt_fractions", "Radiotherapy number of fractions. Number of treatment fractions; an input to rt_curative.", "range 1 to 50 (median 22)", "RTDS (PRESCRIBEDFRACTIONS)", "may be useful to look at, but unlikely to be needed - input to rt_curative",
  "days_dx_to_rt", "Days from diagnosis to radiotherapy. Interval to first RT.", "range 0 to 275 (median 84)", "Derived", "useful - interval to radiotherapy",
  "ORGCODEPROVIDER", "Provider of radiotherapy. Organisation that delivered the RT.", "55 distinct values (id / free-text)", "RTDS (ORGCODEPROVIDER)", "may be useful to look at, but unlikely to be needed - RT-provider only",
  "rcs_ch_score", "RCS Charlson comorbidity score. Weighted score of pre-existing conditions from HES in the year before diagnosis, using the Royal College of Surgeons Charlson method; a continuous risk-adjustment measure.", "range 0 to 3 (median 0)", "Derived from HES APC (RCS Charlson)", "useful - continuous comorbidity; cci_group is the banded version, normally preferred",
  "cci_any", "Any Charlson condition present (0/1). Whether at least one of the 14 Charlson conditions was found.", "range 0 to 1 (median 0)", "Derived from HES APC", "may be useful to look at, but cci_group is richer",
  "cci_group", "Charlson group for analysis: 0, 1, 2, 3+. Banded comorbidity count used as the standard risk-adjustment categorical (matches NOGCA).", "0 | 1 | 2 | 3+", "Derived from HES APC", "essential - the standard comorbidity variable",
  "cci_n_conditions", "Count of Charlson conditions. How many of the 14 conditions were present.", "range 0 to 8 (median 0)", "Derived from HES APC", "may be useful to look at, but cci_group is the banded version",
  "cci_conditions", "Which Charlson conditions were found (text). Names of the specific conditions detected; for checking, not modelling.", "643 distinct values (id / free-text)", "Derived from HES APC", "ignore for modelling - reference only",
  "had_emresd", "Had EMR/ESD in window (flag). TRUE if an endoscopic resection was recorded; a building block of tx_pathway.", "FALSE | TRUE", "Derived", "useful - building block of tx_pathway",
  "had_surgery", "Had a major resection in window (flag). TRUE if a HES resection was recorded; building block of tx_pathway.", "FALSE | TRUE", "Derived", "useful - building block of tx_pathway",
  "had_curative_surgery", "Had a curative-intent resection (flag). The curative subset of had_surgery; a building block behind received_curative_tx_audit.", "FALSE | TRUE", "Derived", "useful - feeds received_curative_tx_audit (the headline flag)",
  "had_sact", "Had systemic therapy in window (flag). TRUE if chemo/immunotherapy was recorded.", "FALSE | TRUE", "Derived", "useful - building block of tx_pathway",
  "had_rt", "Had radiotherapy in window (flag). TRUE if any RT was recorded.", "FALSE | TRUE", "Derived", "useful - building block of tx_pathway",
  "had_curative_rt", "Had curative-intent radiotherapy (flag). Curative subset of had_rt; a building block behind received_curative_tx_audit.", "FALSE | TRUE", "Derived", "useful - feeds received_curative_tx_audit (the headline flag)",
  "had_palliative_rt", "Had palliative radiotherapy (flag). Non-curative RT.", "FALSE | TRUE", "Derived", "may be useful to look at, but unlikely to be needed",
  "had_chemo_for_chemort", "Chemo timed to count toward chemoRT (flag). Chemo present and close enough to RT to form a chemoRT pairing; internal to pathway logic.", "FALSE | TRUE", "Derived", "ignore - internal pathway logic",
  "sact_before_surgery", "Chemo before surgery (flag). Marks the neoadjuvant timing pattern; internal to pathway logic.", "FALSE | TRUE", "Derived", "ignore - internal pathway logic",
  "sact_after_surgery", "Chemo after surgery (flag). Marks the adjuvant timing pattern; internal to pathway logic.", "FALSE | TRUE", "Derived", "ignore - internal pathway logic",
  "rt_before_surgery", "RT before surgery (flag). Neoadjuvant RT timing; internal to pathway logic.", "FALSE | TRUE", "Derived", "ignore - internal pathway logic",
  "rt_after_surgery", "RT after surgery (flag). Adjuvant RT timing; internal to pathway logic.", "FALSE | TRUE", "Derived", "ignore - internal pathway logic",
  "concurrent_chemo_rt", "Concurrent chemo+RT (flag). Chemo and RT within the concurrency window, defining definitive chemoRT; a building block behind received_curative_tx_audit.", "FALSE | TRUE", "Derived", "useful - feeds received_curative_tx_audit (the headline flag)",
  "received_curative_tx", "Received any curative treatment (build-level flag). The internal curative flag from the pathway build; the audit version (received_curative_tx_audit) is the one to report.", "FALSE | TRUE", "Derived", "ignore - use received_curative_tx_audit",
  "tx_pathway", "Treatment pathway - the main treatment variable. One of 14 mutually exclusive categories summarising what curative/non-curative treatment the patient received and in what combination (see the Pathways tab). The headline treatment classification for most analyses.", "Curative RT only | ... | Surgery only", "Derived (treatment anchors + timing)", "essential - the main treatment variable",
  "first_tx_date", "Date of first curative treatment. The clock-stop date for curative pathways (earliest of EMR/surgery/chemo/RT as applicable); NA for non-curative pathways.", "range 2015-01-01 to 2023-08-02", "Derived", "useful - clock-stop date for curative waiting times",
  "tx_trust", "Treating trust for the curative clock-stop. Trust credited with the curative treatment (surgery, EMR or RT provider as appropriate); NA for non-curative patients. This is TRUST-level, not hospital site: SACT and RTDS report only the 5-digit site code inconsistently, so treatment is resolved to trust to keep it comparable across modalities. The unit for trust-level variation analysis.", "133 distinct values (id / free-text)", "Derived (PROCODE3 / emresd_provider / ORGCODEPROVIDER)", "essential - the treating trust for variation analysis (trust-level, not hospital site)",
  "wt_dx_to_tx", "Waiting time: diagnosis to first treatment (days). Days from diagnosis to first curative treatment.", "range -29 to 274 (median 57)", "Derived", "useful - diagnosis-to-treatment wait, but unlikely to be the primary interval",
  "wt_endo_to_tx", "Waiting time: endoscopy to first treatment (days). NOGCA-style wait from diagnostic endoscopy to first treatment.", "range -29 to 292 (median 60)", "Derived", "useful - endoscopy-to-treatment wait, but unlikely to be the primary interval",
  "wt_dx_to_surg", "Waiting time: diagnosis to surgery (days).", "range -16 to 275 (median 149)", "Derived", "may be useful to look at, but unlikely to be needed - modality-specific wait",
  "wt_endo_to_surg", "Waiting time: endoscopy to surgery (days).", "range 0 to 292 (median 151)", "Derived", "may be useful to look at, but unlikely to be needed - modality-specific wait",
  "wt_dx_to_sact", "Waiting time: diagnosis to chemo (days).", "range -29 to 275 (median 58)", "Derived", "may be useful to look at, but unlikely to be needed - modality-specific wait",
  "wt_endo_to_sact", "Waiting time: endoscopy to chemo (days).", "range -29 to 287 (median 60)", "Derived", "may be useful to look at, but unlikely to be needed - modality-specific wait",
  "wt_sact_to_surg", "Waiting time: chemo to surgery (days). Interval between neoadjuvant chemo and resection.", "range -243 to 239 (median 101)", "Derived", "may be useful to look at, but unlikely to be needed - neoadjuvant patients only",
  "wt_surg_to_sact", "Waiting time: surgery to chemo (days). Interval for adjuvant chemo.", "range -239 to 243 (median -101)", "Derived", "may be useful to look at, but unlikely to be needed - adjuvant patients only",
  "wt_dx_to_rt", "Waiting time: diagnosis to radiotherapy (days).", "range 0 to 275 (median 84)", "Derived", "may be useful to look at, but unlikely to be needed - modality-specific wait",
  "wt_endo_to_rt", "Waiting time: endoscopy to radiotherapy (days).", "range 2 to 298 (median 86)", "Derived", "may be useful to look at, but unlikely to be needed - modality-specific wait",
  "wt_rt_to_surg", "Waiting time: radiotherapy to surgery (days). Sparse - only the few neoadjuvant-RT-then-surgery patients.", "range -250 to 228 (median 88)", "Derived", "may be useful to look at, but unlikely to be needed - very low completeness",
  "surv_from_surg_days", "Days survived from surgery. Survival measured from the resection date; surgical patients only.", "range 0 to 3265 (median 927)", "Derived (ONS)", "useful - post-surgery survival, surgical patients only",
  "alive_90d_post_surg", "Alive 90 days after surgery (flag). NOGCA Performance Indicator 7 (post-operative mortality).", "FALSE | TRUE", "Derived", "useful - post-operative mortality (surgical patients only)",
  "alive_1yr_post_surg", "Alive 1 year after surgery (flag). NOGCA Performance Indicator 8.", "FALSE | TRUE", "Derived", "useful - 1-year post-surgery survival (surgical patients only)",
  "cwt_dtt_date", "CWT decision-to-treat date. Date the decision to treat was recorded in Cancer Waiting Times; the start of the decision-to-treatment wait.", "range 2014-12-19 to 2023-08-11", "CWT (TREAT_PERIOD_START)", "useful - decision-to-treat date",
  "cwt_mdt_date", "CWT MDT meeting date. Date the case was discussed at the multidisciplinary team meeting.", "range 2000-01-01 to 2021-03-17", "CWT (MDT_DATE)", "may be useful to look at, but unlikely to be needed - incomplete",
  "cwt_treat_date", "CWT treatment-start date. Treatment start date as recorded in CWT. Note: CWT is used for DATES only - NOT to decide whether treatment happened (HES/SACT/RTDS do that).", "range 2014-12-29 to 2023-08-29", "CWT (TREAT_START)", "useful - treatment date (for waits, not for ascertaining treatment)",
  "cwt_modality", "CWT treatment modality code. The kind of treatment CWT recorded (01/23/24 surgery, 02/04/05 chemo/RT, 07/08/09 palliative care, etc.); used to align the CWT date with the right treatment, not to define treatment.", "01 | 02 | ... | 97", "CWT (MODALITY)", "may be useful to look at, but unlikely to be needed - alignment only, see Sources tab",
  "wt_endo_to_dtt", "Waiting time: endoscopy to decision-to-treat (days). The core decision-to-treat wait, measured from the diagnostic endoscopy - the audit-aligned clock start, preferred over the diagnosis-date version (wt_dx_to_dtt).", "range -30 to 287 (median 43)", "Derived", "essential - the core endoscopy-to-decision-to-treat wait (audit-aligned clock start)",
  "wt_dtt_to_tx", "Waiting time: decision-to-treat to first treatment (days). The core decision-to-treatment interval.", "range -268 to 274 (median 13)", "Derived", "essential - a core waiting-time interval (decision-to-treat to treatment)",
  "wt_dx_to_dtt", "Waiting time: diagnosis to decision-to-treat (days). The interval from the registry diagnosis date to the decision to treat; wt_endo_to_dtt (from the diagnostic endoscopy) is the better clock start and aligns with the audit standard, so prefer that for the core analysis.", "range -30 to 270 (median 42)", "Derived", "useful - but prefer wt_endo_to_dtt (endoscopy is the audit-aligned clock start)",
  "dtt_valid", "Decision-to-treat record is usable (flag). TRUE/FALSE where the CWT decision-to-treat can be trusted for waiting times; NA where there is no curative first treatment to validate against. Filter to dtt_valid == TRUE before reporting decision-to-treatment waits.", "FALSE | TRUE", "Derived", "essential - the filter for clean decision-to-treat waits",
  "tx_modality_audit", "Audit treatment-modality category. A higher-level grouping that matches NOGCA's reporting categories (Surgery only, Surgery plus SACT/RT, Definitive chemoRT, EMR/ESD, Curative RT only, Chemo/RT only (non-curative), No treatment recorded). For figures comparable to the published audit; tx_pathway gives the finer 14-way detail.", "Chemo/RT only (non-curative) | ... | Surgery plus SACT/RT", "Derived", "useful - for audit-comparable reporting; tx_pathway gives the detail",
  "tx_intent_audit", "Audit treatment intent: Curative / Non-curative / No treatment. A three-way summary of treatment intent.", "Curative | No treatment | Non-curative", "Derived", "useful - simple curative/non-curative split; received_curative_tx_audit is the headline flag",
  "received_any_tx", "Received any treatment (flag). TRUE if the patient had any treatment at all (curative or not); FALSE only for 'No treatment recorded'.", "FALSE | TRUE", "Derived", "useful - any-treatment flag",
  "received_curative_tx_audit", "Received curative treatment - the headline curative flag. TRUE if the patient received any curative treatment (EMR/ESD, curative surgery, definitive chemoRT, or curative RT). The single flag for defining the curative cohort; matches the audit's curative definition.", "FALSE | TRUE", "Derived", "essential - the headline curative-cohort flag")

# core = essential for the epidemiological analysis (green). Everything not core
# is "secondary" (amber) unless its recommendation starts "ignore" (grey).
core_vars <- c("pseudo_patientid", "diagmdy", "ydiag", "sex", "agediag", "ethnicity_group_broad", "NHSE_reversed_imd_quintile_lsoas", "route_combined", "tumour_site_grp", "cancer_subtype", "stage_clean", "diag_trust", "diag_hosp", "cci_group", "tx_pathway", "tx_trust", "wt_endo_to_dtt", "wt_dtt_to_tx", "dtt_valid", "received_curative_tx_audit")

# -----------------------------------------------------------------------------
# guard: the reviewed content must still match the live cohort exactly
# -----------------------------------------------------------------------------
missing_from_dict <- setdiff(names(og), dict$variable)
extra_in_dict     <- setdiff(dict$variable, names(og))
if (length(missing_from_dict))
  stop("cohort variables missing from dict: ",
       paste(missing_from_dict, collapse = ", "), call. = FALSE)
if (length(extra_in_dict))
  stop("dict has variables not in the cohort: ",
       paste(extra_in_dict, collapse = ", "), call. = FALSE)

dict <- dict[match(names(og), dict$variable), ]

# live type and completeness
dict$type         <- vapply(og[dict$variable], function(x) class(x)[1], character(1))
dict$pct_complete <- round(100 * vapply(og[dict$variable],
                                        function(x) mean(!is.na(x)), numeric(1)), 1)

# tier
dictionary <- dict %>%
  mutate(tier = case_when(
    variable %in% core_vars                              ~ "core",
    grepl("^ignore", recommendation, ignore.case = TRUE) ~ "ignore",
    TRUE                                                 ~ "secondary")) %>%
  select(variable, tier, definition, type, pct_complete, values,
         source, recommendation)

# -----------------------------------------------------------------------------
# preamble: short guidance blocks shown ABOVE the table, in sentence case.
# Written as heading + text pairs into the top rows, light blue, so they orient
# the reader without taking a column slot or pushing the variables off-screen.
preamble <- tibble::tribble(
  ~heading, ~text,
  "About this cohort",
  "Oesophago-gastric cancer (ICD-10 C15 oesophagus, C16 stomach), stage 1-3, diagnosed 2015-2022, England. One row per patient. Sources: NCRAS cancer registry (demographics, tumour, stage, route, survival); HES APC/OP (surgery, EMR/ESD, endoscopy, comorbidity); SACT (chemo); RTDS (radiotherapy); CWT (waiting-time dates only).",
  
  "How to read the shading",
  "Green = core: essential for the epidemiological analysis, reach for these first. Amber = secondary: useful but mostly feeds a core variable or is modality-specific. Grey = ignore: raw inputs or helpers, not for analysis. The recommendation column gives the per-variable steer.",
  
  "Subsetting to the curative cohort",
  "Simplest: filter received_curative_tx_audit == TRUE (matches the audit's curative definition). Equivalently tx_intent_audit == 'Curative', or keep the curative tx_pathway categories. Curative = EMR/ESD, curative-intent surgery (excludes stage-4 partial gastrectomy per NOGCA), definitive chemoRT, or curative-dose RT. Do not define curative treatment from the raw registry flags (sg_flag/rt_flag/ct_flag) or raw SACT inputs.",
  
  "Core waiting times",
  "The core intervals are wt_endo_to_dtt (endoscopy to the decision to treat - the audit-aligned clock start) and wt_dtt_to_tx (decision to treat to treatment). Filter dtt_valid == TRUE before reporting decision-to-treatment waits. The diagnosis-based version wt_dx_to_dtt, the modality-specific waits, the other diagnosis-to-X waits and days_endo_to_dx are secondary - useful context but not core.",
  
  "Key method notes",
  "Surgery is ascertained from HES APC (gold standard; OPCS-4 matches NOGCA Appendix 8) - CWT surgery codes without a HES resection are not counted. Treatments counted up to 275 days (~9 months) after diagnosis. IMD here is reversed (1 = most deprived). Use route_combined for route to diagnosis, and exclude emergency cases by dropping its 'Emergency presentation' category. References: NOGCA State of the Nation 2025 Methodology Supplement; NDRS DARS data dictionary; NHS England HES Data Dictionary.",
  
  "Treatment pathways (tx_pathway)",
  "Curative: surgery only; surgery + neoadjuvant chemo / RT / chemoRT; surgery + adjuvant chemo; surgery + other; EMR/ESD only; EMR/ESD then surgery; definitive chemoRT; curative RT only. Non-curative: palliative chemo + RT; SACT only; palliative RT only. And: no treatment recorded."
)

# -----------------------------------------------------------------------------
# write the single-sheet workbook: preamble block, then the variable table
# -----------------------------------------------------------------------------
out <- file.path(dir_icon, "og_cohort_data_dictionary.xlsx")
wb  <- createWorkbook()
addWorksheet(wb, "Dictionary")

ncol_tbl    <- ncol(dictionary)               # table is this many columns wide
table_start <- nrow(preamble) + 2             # leave a blank row under the preamble

# --- preamble styles (light blue) ---
pre_head <- createStyle(fgFill = "#EAF1F8", fontColour = "#1F4E79",
                        textDecoration = "bold", valign = "top", wrapText = TRUE,
                        border = "TopBottomLeftRight", borderColour = "#D2DCEA")
pre_text <- createStyle(fgFill = "#F4F8FC", fontColour = "#33475B",
                        valign = "top", wrapText = TRUE,
                        border = "TopBottomLeftRight", borderColour = "#D2DCEA")

for (i in seq_len(nrow(preamble))) {
  writeData(wb, "Dictionary", preamble$heading[i], startCol = 1, startRow = i, colNames = FALSE)
  writeData(wb, "Dictionary", preamble$text[i],    startCol = 2, startRow = i, colNames = FALSE)
  addStyle(wb, "Dictionary", pre_head, rows = i, cols = 1, stack = TRUE)
  # the text spans the remaining columns for readability
  mergeCells(wb, "Dictionary", cols = 2:ncol_tbl, rows = i)
  addStyle(wb, "Dictionary", pre_text, rows = i, cols = 2:ncol_tbl,
           gridExpand = TRUE, stack = TRUE)
}

# --- the variable table, starting below the preamble ---
hdr  <- createStyle(textDecoration = "bold", valign = "top", halign = "left",
                    fgFill = "#E8EEF4", border = "bottom")
writeData(wb, "Dictionary", dictionary, startRow = table_start, headerStyle = hdr)
setColWidths(wb, "Dictionary", cols = seq_len(ncol_tbl),
             widths = c(variable = 26, tier = 11, definition = 82, type = 10,
                        pct_complete = 12, values = 40, source = 38,
                        recommendation = 50))

# wrap + light dashed bottom border on every table body row, for legibility
body_rows <- (table_start + 1):(table_start + nrow(dictionary))
row_line  <- createStyle(wrapText = TRUE, valign = "top",
                         border = "bottom", borderColour = "#D9D9D9",
                         borderStyle = "dashed")
addStyle(wb, "Dictionary", row_line, rows = body_rows,
         cols = seq_len(ncol_tbl), gridExpand = TRUE, stack = TRUE)

# freeze just below the table header so the preamble scrolls away but the column
# headers stay visible while reading the variables
freezePane(wb, "Dictionary", firstActiveRow = table_start + 1)

# soft traffic-light shading (muted so amber does not read as a warning)
core_fill   <- createStyle(fgFill = "#DDEBD8", wrapText = TRUE, valign = "top")  # soft green
sec_fill    <- createStyle(fgFill = "#FBF0DD", wrapText = TRUE, valign = "top")  # soft amber
ignore_fill <- createStyle(fgFill = "#EFEFEF", wrapText = TRUE, valign = "top")  # soft grey
shade <- function(t, style) {
  r <- which(dictionary$tier == t) + table_start   # offset to the table body
  if (length(r)) addStyle(wb, "Dictionary", style, rows = r,
                          cols = seq_len(ncol_tbl), gridExpand = TRUE, stack = TRUE)
}
shade("core", core_fill)
shade("secondary", sec_fill)
shade("ignore", ignore_fill)

saveWorkbook(wb, out, overwrite = TRUE)
cat("data dictionary written ->", out, "\n")
cat("variables:", nrow(dictionary),
    "| core (green):", sum(dictionary$tier == "core"),
    " secondary (amber):", sum(dictionary$tier == "secondary"),
    " ignore (grey):", sum(dictionary$tier == "ignore"), "\n")