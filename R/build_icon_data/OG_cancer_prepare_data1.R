# =============================================================================
# OG Cancer Treatment Waiting Times Pipeline
# =============================================================================

library(arrow)
library(haven)
library(dplyr)
library(tidyr)
library(stringr)
library(forcats)
library(purrr)
library(lubridate)


# =============================================================================
# 1. NCRAS: OG cancer cohort
#    ICD-10: C15.x (oesophagus) or C16.x (stomach/gastric)
#    Field names match the NCRAS clean parquet used in the colon pipeline
#    Inclusion/exclusion criteria per NOGCA Methodology Supplement ?3
# =============================================================================

# --- 1a. Column selection (mirrors colon script ncras_cols) -----------------
ncras_cols <- c(
  # Identity and linkage
  "pseudo_patientid", "pseudo_tumourid",
  
  # Diagnosis date and year
  "diagmdy", "ydiag",
  
  # Tumour characteristics
  "cancer", "sitestr", "typestr",
  "basisofdiagnosis",
  "grade", "behav",
  
  # Morphology - field name varies by parquet version:
  # may be 'morphology', 'tumour_morphology', or encoded within 'typestr'
  # Adjust the coalesce logic in ?1f if needed
  "morphology",
  
  # Staging
  "stage_best", "stage_best_system",
  "t_best", "n_best", "m_best",
  "t_path", "n_path", "m_path",
  
  # Patient characteristics
  "sex", "agediag", "birthmdy",
  "ethnicity_group_broad",
  
  # Geography and deprivation
  "lsoa11_code",
  "NHSE_reversed_imd_quintile_lsoas",
  "canalliance_2024_code", "canalliance_2024_name",
  
  # Organisation
  "diag_trust", "diag_trust_name",
  "first_trust", "first_trust_name", "first_hosp_date",
  "diag_hosp",
  
  # Pathway and route
  "route_bjc", "final_route", "route_code",
  
  # Performance status (risk-adjustment for PI7/PI8; subgroup for PI3/PI6)
  "tumour_performancestatus",
  
  # CNS involvement (PI4) - derived by NDRS/NCRAS from COSD
  "clinicalnursespecialist",
  "firstmdtmeetingdate",
  
  # Treatment flags (broad)
  "sg_flag", "rt_flag", "ct_flag",
  
  # Survival
  "dead", "finmdy", "dco"#,
  
  # Comorbidity fallback (Charlson, pre-calculated; else derive from HES)
#  "chrl_tot_27_03"
)

# --- 1b. ICD-10 site codes --------------------------------------------------
og_icd10 <- c(
  "C15","C150","C151","C152","C153","C154","C155","C158","C159",
  "C16","C160","C161","C162","C163","C164","C165","C166","C168","C169"
)

# --- 1c. Epithelial morphology codes (Appendix 4, NOGCA Methodology 2025) ---
#   Morphology must be 8001-9989 AND in this list (generic 8000 excluded)
morph_epithelial <- c(
  8005, 8010, 8020, 8021, 8032, 8033, 8050, 8051, 8052,
  8070, 8071, 8072, 8073, 8074, 8075, 8076, 8077, 8078,
  8083, 8084,
  8140, 8141, 8142, 8143, 8144, 8145,
  8190, 8210, 8211, 8213, 8214, 8231,
  8255, 8260, 8261, 8262, 8263,
  8310, 8323,
  8430, 8440,
  8480, 8481, 8490,
  8510, 8512,
  8560, 8562,
  8570, 8571, 8572, 8573, 8574, 8576,
  8982
)

# --- 1d. Neuroendocrine morphology codes (Appendix 5) - EXCLUSION -----------
morph_neuroendocrine <- c(
  8013, 8041, 8042, 8043, 8044, 8045,
  8150, 8151, 8152, 8153, 8154, 8155, 8156, 8157, 8158,
  8240, 8241, 8242, 8243, 8244, 8245, 8246, 8247, 8249,
  9091
)


# --- 1e. Read NCRAS partioned parquet ----------------------------------------
# 
#  ncras_raw <- read_parquet(
#    "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/NCRAS/NCRAS_clean_1995_2022_route.parquet"
#  )
# 
# ncras_raw  <- ncras_raw %>%
#   mutate(across(where(is.labelled), ~as.character(as_factor(.x))))
# 
#  # Write partitioned dataset
#  ncras_raw %>%
#    write_dataset(
#      path         = "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/NCRAS/NCRAS_clean_1995_2022_route_sitestr/",
#      format       = "parquet",
#      partitioning = "sitestr"
#    )
#  
 # 
# #  
#  test <- open_dataset(
#     "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/NCRAS/NCRAS_clean_1995_2022_route_sitestr/"
#   ) %>%
#    filter(sitestr %in% og_icd10) %>%
#     collect()
# 
# test$dco[1:100]
# summary(as.factor(test$dco))
# names_ncras <- names(test)
# # #names_ncras[substr(names_ncras, start = 1, stop = 3) == "cli"]

ncras_og_raw <- open_dataset(
  "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/NCRAS/NCRAS_clean_1995_2022_route_sitestr/"
) %>%
  filter(sitestr %in% og_icd10) %>%
  collect()  %>%
  select(any_of(ncras_cols)) %>%
  mutate(
    pseudo_patientid = as.character(pseudo_patientid),
    diagmdy          = as.Date(diagmdy),
    finmdy           = as.Date(finmdy),
    morphology_num   = as.integer(as.character(typestr))
  ) %>% 
  # ---- Derived fields ------------------------------------------------------
mutate(
  
  # Tumour site group
  tumour_site_grp = case_when(
    str_starts(as.character(sitestr), "C15") ~ "oesophageal",
    str_starts(as.character(sitestr), "C16") ~ "gastric",
    TRUE ~ NA_character_
  ),
  
  # Cancer subtype (SCC / ACA; Appendix 3 morphology codes)
  cancer_subtype = case_when(
    tumour_site_grp == "oesophageal" &
      morphology_num %in% c(8033, 8051, 8052,
                            8070, 8071, 8072, 8073, 8074, 8075,
                            8076, 8077, 8078, 8083, 8084) ~ "Oes SCC",
    tumour_site_grp == "oesophageal" &
      morphology_num %in% c(8005, 8140, 8141, 8142, 8143, 8144, 8145,
                            8190, 8210, 8211, 8213, 8214,
                            8255, 8260, 8261, 8262, 8263,
                            8310, 8323, 8440,
                            8480, 8481,
                            8570, 8571, 8572, 8573, 8574, 8576) ~ "Oes ACA",
    tumour_site_grp == "gastric"                              ~ "Gast",
    TRUE                                                      ~ NA_character_
  ),
  
  # Stage derivation (NOGCA ?4 / AJCC v8):
  #   Stage 0 -> 1
  #   Missing -> impute as 4 if metastatic behaviour code or metastasis biopsy
  stage_clean = case_when(
    as.character(stage_best) == "0"           ~ "1",
    str_starts(as.character(stage_best), "1") ~ "1",
    str_starts(as.character(stage_best), "2") ~ "2",
    str_starts(as.character(stage_best), "3") ~ "3",
    str_starts(as.character(stage_best), "4") ~ "4",
    TRUE                                      ~ NA_character_
    # Note: "", "X", "U", "?" and "6" all map to NA (unknown stage)
    # behav /6 not present in this extract; basisofdiagnosis imputation catches 0 cases
  ),
  
  # Route: combine final_route and route_bjc (precedence: final_route)
  # Mirrors colon script route_combined derivation
  final_route_chr = na_if(as.character(final_route), ""),
  route_bjc_chr   = na_if(as.character(route_bjc),   ""),
  route_combined  = factor(
    coalesce(final_route_chr, route_bjc_chr, "Unknown")
  ),
  
  # Emergency admission flag (PI1)
  emergency_admission = as.integer(
    as.character(route_combined) == "Emergency presentation"
  ),
  
  # Performance status (numeric; for risk-adjustment)
  # Need to bring in from COSD.
  #ps_num = as.integer(as.character(tumour_performancestatus)),
  
  # CNS involvement flag (PI4): any "Yes" response in clinicalnursespecialist
  # Full NOGCA derivation also requires firstmdtmeetingdate within 90d;
  # apply that filter below if the field is populated
  #  cnsinvolved = as.integer(
  #    str_detect(tolower(as.character(clinicalnursespecialist)), "yes")
  #  ),
  
  # Survival from diagnosis
  surv_from_dx_days = as.integer(finmdy - diagmdy),
  died = as.integer(dead)

)


# Checks:

# --- Tumour site group -------------------------------------------------------
cat("tumour_site_grp:\n")
print(count(ncras_og_raw, tumour_site_grp))

# --- Cancer subtype ----------------------------------------------------------
cat("\ncancer_subtype:\n")
print(count(ncras_og_raw, cancer_subtype, sort = TRUE))

# --- Stage -------------------------------------------------------------------
cat("\nstage_best (raw):\n")
print(count(ncras_og_raw, stage_best, sort = TRUE) %>% print(n = 20))

cat("\nstage_clean (derived):\n")
print(count(ncras_og_raw, stage_clean))

# Cross-check: NAs in stage_clean - what are the raw stage_best values?
cat("\nstage_best values where stage_clean is NA:\n")
ncras_og_raw %>%
  filter(is.na(stage_clean)) %>%
  count(stage_best, sort = TRUE) %>%
  print(n = 20)

# --- Route -------------------------------------------------------------------
cat("\nfinal_route (raw):\n")
print(count(ncras_og_raw, final_route, sort = TRUE))

cat("\nroute_combined (derived):\n")
print(count(ncras_og_raw, route_combined, sort = TRUE))

cat("\nemergency_admission:\n")
print(count(ncras_og_raw, emergency_admission))

# --- Survival ----------------------------------------------------------------
cat("\ndied:\n")
print(count(ncras_og_raw, died))

cat("\ndead (raw):\n")
print(count(ncras_og_raw, dead))

# Sanity check: anyone with finmdy but died == 0?
cat("\nfinmdy present but died == 0:\n")
ncras_og_raw %>%
  filter(!is.na(finmdy), died == 0) %>%
  nrow() %>%
  cat("\n")

# Survival distribution
cat("\nsurv_from_dx_days summary:\n")
summary(ncras_og_raw$surv_from_dx_days)

# --- Morphology --------------------------------------------------------------
cat("\nmorphology_num - NAs (typestr couldn't be parsed):\n")
cat(sum(is.na(ncras_og_raw$morphology_num)), "\n")

cat("\ntop morphology codes:\n")
print(count(ncras_og_raw, morphology_num, sort = TRUE) %>% print(n = 20))

# --- behav (for stage imputation cross-check) --------------------------------
cat("\nbehav:\n")
print(count(ncras_og_raw, behav, sort = TRUE))

# How many would get stage 4 imputed via basisofdiagnosis?
cat("\nbasisofdiagnosis containing 'metastasis' with missing stage:\n")
ncras_og_raw %>%
  filter(is.na(stage_clean),
         str_detect(tolower(as.character(basisofdiagnosis)), "metastasis")) %>%
  nrow() %>%
  cat("\n")

ncras_og_raw %>%
  filter(!is.na(finmdy), died == 0) %>%
  select(diagmdy, finmdy, dead, surv_from_dx_days) %>%
  head(20)

# --- 1g. Apply inclusion / exclusion criteria (NOGCA ?3) --------------------
ncras_og <- ncras_og_raw %>%
      
# ---- Inclusion: diagnosis year 2015+ ------------------------------------
  filter(ydiag >= 2015) %>%
    
# ---- Inclusion: adults --------------------------------------------------
  filter(agediag >= 18 | is.na(agediag)) %>%
      
# ---- Inclusion: histological diagnosis (morphology 8001-9989) -----------
  filter(!is.na(morphology_num),
           morphology_num >= 8001,
           morphology_num <= 9989) %>%
      
# ---- Inclusion: epithelial tumour (Appendix 4) --------------------------
  filter(morphology_num %in% morph_epithelial) %>%
      
# ---- Exclusion: neuroendocrine morphology (Appendix 5) ------------------
    filter(!morphology_num %in% morph_neuroendocrine) %>%
      
# ---- Exclusion: death certificate only ----------------------------------
# basis=9: unknown/death certificate (264 cases)
# diagmdy==finmdy & died: registered from death, not basis=9 coded (658 cases)
    filter(
      as.integer(basisofdiagnosis) != 9L,
      !(died == 1L & diagmdy == finmdy)
    ) %>%
      
# ---- Retain earliest diagnosis per patient (first primary OG) -----------
    arrange(pseudo_patientid, diagmdy) %>%
    distinct(pseudo_patientid, .keep_all = TRUE) %>%
      
# ---- Inclusion: stage 1-3 only (curative treatment analysis) ------------
    filter(stage_clean %in% c("1", "2", "3"))

cat("After ydiag >= 2015:        ", nrow(filter(ncras_og_raw, ydiag >= 2015)), "\n")
cat("After epithelial morphology:", nrow(filter(ncras_og_raw, ydiag >= 2015, morphology_num %in% morph_epithelial)), "\n")
cat("After DCO exclusions:       ", nrow(ncras_og), "\n")
cat("After deduplication:        ", n_distinct(ncras_og$pseudo_patientid), "\n")

cat("OG cancer patients (NCRAS, 2015+):", n_distinct(ncras_og$pseudo_patientid), "\n")
cat("  Oesophageal:", sum(ncras_og$tumour_site_grp == "oesophageal", na.rm = TRUE), "\n")
cat("  Gastric:    ", sum(ncras_og$tumour_site_grp == "gastric",     na.rm = TRUE), "\n")
cat("  By year:\n")
print(table(ncras_og$ydiag))

cat("\ncancer_subtype:\n")
print(count(ncras_og, cancer_subtype, sort = TRUE))

cat("After stage 1-3 restriction:", nrow(ncras_og), "\n")
cat("  Stage 1:", sum(ncras_og$stage_clean == "1"), "\n")
cat("  Stage 2:", sum(ncras_og$stage_clean == "2"), "\n")
cat("  Stage 3:", sum(ncras_og$stage_clean == "3"), "\n")
cat("  By site:\n")
print(count(ncras_og, tumour_site_grp, stage_clean))

ncras_og_ids <- ncras_og %>% distinct(pseudo_patientid) %>% pull()

# saveRDS(ncras_og,
#         "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/ncras_og_2015_2022.rds")
# 


# =============================================================================
#  2. COSD linkage.
# =============================================================================

cosd <- read_dta("E:/Data_PHE/Raw data files received from PHE READ ONLY/NCRAS/Stata files/18_COSD_data.dta")
#cosd  

# How complete is PS for your OG cohort?

cosd %>%
  filter(pseudo_patientid %in% ncras_og_ids) %>%
  count(performancestatus, sort = TRUE)

# CNS involvement
cosd %>%
  filter(pseudo_patientid %in% ncras_og_ids) %>%
  count(clinicalnursespecialist, sort = TRUE)

# Multiple rows per patient?
cosd %>%
  filter(pseudo_patientid %in% ncras_og_ids) %>%
  count(pseudo_patientid, sort = TRUE) %>%
  filter(n > 1) %>%
  nrow()  

cosd %>%
  filter(pseudo_patientid %in% ncras_og_ids) %>%
  group_by(pseudo_patientid) %>%
  filter(n() > 1) %>%
  select(pseudo_patientid, pseudo_tumourid, 
         performancestatus, clinicalnursespecialist) %>%
  head(20)

cosd_og <- cosd %>%
  filter(pseudo_patientid %in% ncras_og_ids) %>%
  mutate(
    # Performance status: recode 9/"" as missing
    ps_num = case_when(
      performancestatus %in% c("0","1","2","3","4") ~ as.integer(performancestatus),
      TRUE                                          ~ NA_integer_
    ),
    # CNS involvement
    cnsinvolved = case_when(
      str_starts(clinicalnursespecialist, "Y") ~ 1L,
      clinicalnursespecialist == "NN"          ~ 0L,
      TRUE                                     ~ NA_integer_
    )
  ) %>%
  select(pseudo_patientid, pseudo_tumourid, ps_num, cnsinvolved)

# Join on both IDs - exact tumour match
ncras_og <- ncras_og %>%
  left_join(cosd_og, by = c("pseudo_patientid", "pseudo_tumourid"))

# Check completeness
cat("PS completeness:\n")
ncras_og %>%
  summarise(
    n       = n(),
    n_ps    = sum(!is.na(ps_num)),
    pct_ps  = round(100 * n_ps / n, 1),
    n_cns   = sum(!is.na(cnsinvolved)),
    pct_cns = round(100 * n_cns / n, 1)
  ) %>%
  print()

cat("\nPS distribution:\n")
count(ncras_og, ps_num)

cat("\nCNS distribution:\n")
count(ncras_og, cnsinvolved)

# Overwrite saved object with COSD variables included
saveRDS(ncras_og,
        "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/ncras_og_2015_2022.rds")

#-------------------------------------------------------
#  OPCS-4 CODE DEFINITIONS
# ----------------------------------------------------

# --- 2a. Diagnostic endoscopy codes (Appendix 6) ----------------------------
#   First qualifying endoscopy within 30 days BEFORE diagnosis date
#   = start of PI3 time-to-treatment interval

opcs_diagnostic_endoscopy <- c(
  # Fibreoptic oesophageal - therapeutic codes present in diagnostic pathway
  "G142","G143","G145","G147",
  "G152","G153","G154","G156","G157","G158","G159",
  # Fibreoptic oesophageal - diagnostic
  "G161","G162","G168","G169",
  # Rigid oesophagoscope - therapeutic
  "G172","G173","G188","G189",
  # Rigid oesophagoscope - diagnostic
  "G191","G198","G199",
  # Oesophageal haemostasis / other
  "G201","G202","G208","G209",
  "G214","G215","G218","G219",
  # Upper GI fibreoptic - therapeutic
  "G422","G432","G433","G435",
  "G441","G443","G445","G446","G448","G449",
  # Upper GI fibreoptic - diagnostic
  "G451","G452","G454","G458","G459",
  # Upper GI haemostasis / other
  "G462","G463","G468","G469"
)

# --- 2b. EMR / ESD codes (Appendix 7) ---------------------------------------
#   Curative endotherapy; window: 30d before to 9 months after diagnosis
opcs_emresd <- c(
  # Oesophageal EMR/ESD
  "G121","G128","G129",
  "G141","G146","G148","G149",
  "G171","G178","G179",
  # Upper GI EMR/ESD
  "G421","G423","G428","G429",
  "G431","G438","G439",
  # Ablation codes (HGD; included to avoid missed curative procedures)
  "G143","G145","G433","G435"
)

# --- 2c. Major OG resection codes (Appendix 8) ------------------------------
#   Curative intent = any code EXCEPT stage 4 + partial gastrectomy (G28x)
opcs_og_surgery <- list(
  
  oesophagectomy = c(
    # Oesophagogastrectomy - oesophageal cancer flank
    "G011","G018","G019",
    # Total oesophagectomy
    "G021","G022","G023","G024","G025","G028","G029",
    # Partial oesophagectomy
    "G031","G032","G033","G034","G035","G036","G038","G039"
  ),
  
  oesophagogastrectomy_jejunum = c(
    # Oesophagogastrectomy with jejunal reconstruction - gastric cancer flank
    "G012","G013"
  ),
  
  total_gastrectomy = c(
    "G271","G272","G273","G274","G275","G278","G279"
  ),
  
  partial_gastrectomy = c(
    "G281","G282","G283","G288","G289"
  )
)

opcs_og_surgery_all    <- unique(unlist(opcs_og_surgery))
opcs_oesophagectomy    <- unique(c(opcs_og_surgery$oesophagectomy,
                                   opcs_og_surgery$oesophagogastrectomy_jejunum))
opcs_gastrectomy_total <- opcs_og_surgery$total_gastrectomy
opcs_gastrectomy_part  <- opcs_og_surgery$partial_gastrectomy

opcs_og_surgery_lookup <- c(
  setNames(rep("oesophagectomy",     length(opcs_oesophagectomy)),    opcs_oesophagectomy),
  setNames(rep("total_gastrectomy",  length(opcs_gastrectomy_total)), opcs_gastrectomy_total),
  setNames(rep("partial_gastrectomy",length(opcs_gastrectomy_part)),  opcs_gastrectomy_part)
)

# Emergency admission method codes (HES)
admimeth_emerg <- c("21","22","23","24","25","28","2A","2B","2C","2D")

# Treatment window: 9 months post-diagnosis (NOGCA standard; ~275 days)
tx_window_days <- 275


# =============================================================================
# 3. READ HES APC
#    Restricted to OG cohort patients
#    2014-2024 to capture pre-diagnosis endoscopies (30d window)
#    and full post-diagnosis treatment follow-up through to 2023 diagnoses
# =============================================================================

ncras_og     <- readRDS("E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/ncras_og_2015_2022.rds")
ncras_og_ids <- ncras_og %>% distinct(pseudo_patientid) %>% pull()

hes_apc_file_list <- list.files(
  "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/HES/APC/",
  pattern    = "FILE*",
  full.names = TRUE
) %>%
  keep(~{
    yr <- str_extract(.x, "(?<=HES_APC_)\\d{4}") %>% as.integer()
    !is.na(yr) && yr %in% 2014:2024
  })

stopifnot(length(hes_apc_file_list) > 0)

op_cols     <- paste0("OPERTN_", str_pad(1:24, 2, pad = "0"))
opdate_cols <- paste0("OPDATE_", str_pad(1:24, 2, pad = "0"))

hes_cols_select <- c(
  "STUDY_ID", "ADMIDATE", "ADMIMETH", "PROCODE3", "SITETRET",
  "EPISTART", "EPIORDER", "EPITYPE",
  op_cols, opdate_cols
)

hes_apc_raw <- map_dfr(
  hes_apc_file_list,
  ~{
    read_parquet(.x, col_select = all_of(hes_cols_select)) %>%
      filter(STUDY_ID %in% ncras_og_ids) %>%
      mutate(
        STUDY_ID = as.character(STUDY_ID),
        ADMIMETH = as.character(ADMIMETH),
        EPISTART = as.Date(EPISTART),
        ADMIDATE = as.Date(ADMIDATE),
        across(all_of(op_cols),     as.character),
        across(all_of(opdate_cols), as.Date)
      )
  },
  .progress = TRUE
)


saveRDS(
  hes_apc_raw,
  "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/hes_apc_og_2014_2022.rds"
)

# Then subsequent runs start from:
hes_apc_raw <- readRDS(
  "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/hes_apc_og_2014_2022.rds"
)

names(hes_apc_raw)
View(hes_apc_raw[1:100,])
# Helper: normalise OPCS codes for consistent matching
normalise_opcs <- function(x) str_replace_all(str_to_upper(as.character(x)), "\\.", "")

cat("HES APC rows (OG cohort):", nrow(hes_apc_raw), "\n")
cat("Patients:                ", n_distinct(hes_apc_raw$STUDY_ID), "\n")


# =============================================================================
# 4. SHARED HELPER: PIVOT OPERATIONS + DATES LONG, MATCH OPCS LIST
#    Avoids repeating the pivot pattern for endoscopy, EMR/ESD, and surgery
# =============================================================================

match_opcs_episodes <- function(hes_data, opcs_list, op_cols, opdate_cols) {
  
  ops_long <- hes_data %>%
    filter(!is.na(OPERTN_01), OPERTN_01 != "-") %>%
    pivot_longer(cols      = all_of(op_cols),
                 names_to  = "op_position",
                 values_to = "opcs_code") %>%
    filter(!is.na(opcs_code), opcs_code != "-") %>%
    mutate(
      opcs4         = normalise_opcs(opcs_code),
      op_position_n = as.integer(str_extract(op_position, "[0-9]+"))
    ) %>%
    filter(opcs4 %in% opcs_list)
  
  if (nrow(ops_long) == 0) return(tibble())
  
  dates_long <- hes_data %>%
    pivot_longer(cols      = all_of(opdate_cols),
                 names_to  = "opdate_position",
                 values_to = "op_date") %>%
    mutate(
      op_position_n = as.integer(str_extract(opdate_position, "[0-9]+")),
      op_date       = as.Date(op_date)
    ) %>%
    select(STUDY_ID, EPISTART, EPIORDER, op_position_n, op_date)
  
  ops_long %>%
    left_join(dates_long,
              by           = c("STUDY_ID","EPISTART","EPIORDER","op_position_n"),
              relationship = "many-to-many") %>%
    rename(pseudo_patientid = STUDY_ID)
}

# =============================================================================
# 5a. IDENTIFY DIAGNOSTIC ENDOSCOPY HES-APC (Appendix 6)
#    First qualifying endoscopy within 30 days BEFORE diagnosis date
#    = PI3 clock start
# =============================================================================

hes_endoscopy <- match_opcs_episodes(
  hes_apc_raw, opcs_diagnostic_endoscopy, op_cols, opdate_cols
)

endoscopy_anchor <- ncras_og %>%
  select(pseudo_patientid, diagmdy) %>%
  left_join(
    hes_endoscopy %>% select(pseudo_patientid, EPISTART, EPIORDER, opcs4, op_date),
    by = "pseudo_patientid"
  ) %>%
  mutate(days_endo_to_dx = as.integer(diagmdy - op_date)) %>%
  filter(!is.na(days_endo_to_dx),
         days_endo_to_dx >= 0,    # endoscopy on or before diagnosis
         days_endo_to_dx <= 30) %>%
  arrange(pseudo_patientid, op_date) %>%
  distinct(pseudo_patientid, .keep_all = TRUE) %>%
  rename(endoscopy_date = op_date) %>%
  select(pseudo_patientid, endoscopy_date, days_endo_to_dx)

cat("Patients with diagnostic endoscopy <=30d before diagnosis:",
    n_distinct(endoscopy_anchor$pseudo_patientid), "\n")

# =============================================================================
# 5b. DIAGNOSTIC ENDOSCOPY ANCHOR - HES-OP
#     Supplementary endoscopy capture from outpatient records
#     Attended appointments only; OPCS codes per Appendix 6
#     Date: APPTDATE (yyyy-mm-dd format in HES-OP)
# =============================================================================

hes_op_file_list <- list.files(
  "E:/Data_PHE/Raw data files received from PHE READ ONLY/HES/OP/",
  pattern    = "*.txt",
  full.names = TRUE
) %>%
  keep(~{
    yr <- str_extract(.x, "(?<=HES_OP_)\\d{4}") %>% as.integer()
    !is.na(yr) && yr >= 2014
  })

cat("HES-OP files (2014+):", length(hes_op_file_list), "\n")

op_cols_select <- c(
  "STUDY_ID", "APPTDATE", "ATTENDED",
  paste0("OPERTN_0", 1:9),
  paste0("OPERTN_", 10:24),
  "PROCODET", "TRETSPEF", "MAINSPEF"
)

hes_op_raw <- map_dfr(
  hes_op_file_list,
  ~{
    read_delim(
      .x,
      delim          = "|",
      col_select     = any_of(op_cols_select),
      col_types      = cols(.default = col_character()),
      show_col_types = FALSE
    ) %>%
      filter(
        STUDY_ID %in% ncras_og_ids,
        ATTENDED %in% c("5", "6")
      ) %>%
      mutate(
        STUDY_ID  = as.character(STUDY_ID),
        appt_date = as.Date(APPTDATE)
      )
  },
  .progress = TRUE
)

saveRDS(
  hes_op_raw,
  "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/hes_op_og_2014_2022.rds"
)

hes_op_raw <- readRDS( "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/hes_op_og_2014_2022.rds")
hes_op_raw %>% View()


cat("HES-OP rows (OG cohort, attended):", nrow(hes_op_raw), "\n")
cat("Patients:                         ", n_distinct(hes_op_raw$STUDY_ID), "\n")

# --- Match endoscopy OPCS codes --------------------------------------------
op_op_cols <- names(hes_op_raw)[str_starts(names(hes_op_raw), "OPERTN_")]

hes_op_endoscopy <- hes_op_raw %>%
  pivot_longer(
    cols      = all_of(op_op_cols),
    names_to  = "op_position",
    values_to = "opcs_code"
  ) %>%
  filter(!is.na(opcs_code), opcs_code != "-") %>%
  mutate(opcs4 = normalise_opcs(opcs_code)) %>%
  filter(opcs4 %in% opcs_diagnostic_endoscopy) %>%
  select(STUDY_ID, appt_date, opcs4)

cat("HES-OP endoscopy records:         ", nrow(hes_op_endoscopy), "\n")
cat("Patients with OP endoscopy:       ", 
    n_distinct(hes_op_endoscopy$STUDY_ID), "\n")


# --- Build OP endoscopy anchor ---------------------------------------------
# First qualifying endoscopy within 30 days BEFORE diagnosis date
endoscopy_anchor_op <- ncras_og %>%
  select(pseudo_patientid, diagmdy) %>%
  left_join(
    hes_op_endoscopy %>%
      rename(pseudo_patientid = STUDY_ID,
             endoscopy_date   = appt_date),
    by = "pseudo_patientid"
  ) %>%
  mutate(days_endo_to_dx = as.integer(diagmdy - endoscopy_date)) %>%
  filter(
    !is.na(days_endo_to_dx),
    days_endo_to_dx >= 0,
    days_endo_to_dx <= 30
  ) %>%
  arrange(pseudo_patientid, endoscopy_date) %>%
  distinct(pseudo_patientid, .keep_all = TRUE) %>%
  select(pseudo_patientid, endoscopy_date, days_endo_to_dx) %>%
  mutate(endo_source = "OP")

cat("OP endoscopy anchor patients:     ",
    n_distinct(endoscopy_anchor_op$pseudo_patientid), "\n")


cat("Already in APC anchor:", 
    sum(endoscopy_anchor_op$pseudo_patientid %in% 
          endoscopy_anchor$pseudo_patientid), "\n")

cat("New from OP only:     ",
    sum(!endoscopy_anchor_op$pseudo_patientid %in% 
          endoscopy_anchor$pseudo_patientid), "\n")


endoscopy_anchor_combined <- endoscopy_anchor %>%
  mutate(endo_source = "APC") %>%
  bind_rows(
    endoscopy_anchor_op %>%
      filter(!pseudo_patientid %in% endoscopy_anchor$pseudo_patientid)
  ) %>%
  arrange(pseudo_patientid, endoscopy_date) %>%
  distinct(pseudo_patientid, .keep_all = TRUE)

cat("Combined endoscopy anchor:", 
    n_distinct(endoscopy_anchor_combined$pseudo_patientid), "\n")
cat("  APC: ", sum(endoscopy_anchor_combined$endo_source == "APC"), "\n")
cat("  OP:  ", sum(endoscopy_anchor_combined$endo_source == "OP"),  "\n")

saveRDS(endoscopy_anchor_combined, "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/OG_endoscopy_anchor_combined.rds")

# =============================================================================
# 6. IDENTIFY EMR / ESD (Appendix 7)
#    Window: 30 days before to 9 months after diagnosis
# =============================================================================

hes_emresd <- match_opcs_episodes(
  hes_apc_raw, opcs_emresd, op_cols, opdate_cols
)

emresd_anchor <- ncras_og %>%
  select(pseudo_patientid, diagmdy) %>%
  left_join(
    hes_emresd %>% select(pseudo_patientid, EPISTART, EPIORDER, opcs4, op_date),
    by = "pseudo_patientid"
  ) %>%
  mutate(days_dx_to_emresd = as.integer(op_date - diagmdy)) %>%
  filter(!is.na(days_dx_to_emresd),
         days_dx_to_emresd >= -30,
         days_dx_to_emresd <= tx_window_days) %>%
  arrange(pseudo_patientid, op_date) %>%
  distinct(pseudo_patientid, .keep_all = TRUE) %>%
  rename(emresd_date = op_date) %>%
  select(pseudo_patientid, emresd_date, days_dx_to_emresd)


saveRDS(emresd_anchor, "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/OG_emresd_anchor.rds")


cat("Patients with EMR/ESD:", n_distinct(emresd_anchor$pseudo_patientid), "\n")

# What stage and subtype are EMR/ESD patients?
ncras_og %>%
  filter(pseudo_patientid %in% emresd_anchor$pseudo_patientid) %>%
  count(cancer_subtype, stage_clean)

# How many have both endoscopy and EMR/ESD?
sum(emresd_anchor$pseudo_patientid %in% endoscopy_anchor$pseudo_patientid)

# Time from diagnosis to EMR/ESD distribution
summary(emresd_anchor$days_dx_to_emresd)

# =============================================================================
# 7. IDENTIFY MAJOR OG SURGICAL RESECTION (Appendix 8)
#    Elective admissions for performance analyses
#    Window: 30 days before to 9 months after diagnosis
#    Curative intent: excludes stage 4 + partial gastrectomy (G28x)
# =============================================================================

hes_og_surgery <- match_opcs_episodes(
  hes_apc_raw, opcs_og_surgery_all, op_cols, opdate_cols
)

proc_priority_og <- c("oesophagectomy","total_gastrectomy","partial_gastrectomy")

first_or_na      <- function(x) if (length(x) == 0) NA_character_ else dplyr::first(x)

# Collapse to one row per episode; hierarchy: oesophagectomy > total > partial
hes_og_surgery_episodes <- hes_og_surgery %>%
  mutate(surgery_type = unname(opcs_og_surgery_lookup[opcs4])) %>%
  arrange(pseudo_patientid, EPISTART, op_position_n) %>%
  group_by(pseudo_patientid, ADMIDATE, EPISTART, EPIORDER,
           EPITYPE, PROCODE3, SITETRET, ADMIMETH) %>%
  summarise(
    surgery_type = first_or_na(intersect(proc_priority_og, unique(surgery_type))),
    opcs_primary = first_or_na(opcs4),
    all_og_opcs  = paste(unique(opcs4), collapse = "; "),
    surgery_date = as.Date(first_or_na(as.character(op_date[!is.na(op_date)]))),
    .groups      = "drop"
  ) %>%
  mutate(emergency = ADMIMETH %in% admimeth_emerg)

cat("OG resection episodes:", nrow(hes_og_surgery_episodes),
    "| Patients:", n_distinct(hes_og_surgery_episodes$pseudo_patientid), "\n")

# Who has multiple resection episodes?
hes_og_surgery_episodes %>%
  count(pseudo_patientid, sort = TRUE) %>%
  filter(n > 1) %>%
  count(n, name = "n_patients")

# What surgery types are involved in duplicates?
hes_og_surgery_episodes %>%
  add_count(pseudo_patientid) %>%
  filter(n > 1) %>%
  count(surgery_type, sort = TRUE)

# How far apart are the duplicate episodes?
hes_og_surgery_episodes %>%
  group_by(pseudo_patientid) %>%
  filter(n() > 1) %>%
  arrange(pseudo_patientid, surgery_date) %>%
  mutate(days_between = as.integer(surgery_date - lag(surgery_date))) %>%
  filter(!is.na(days_between)) %>%
  summary(days_between)

hes_og_surgery_episodes %>%
  group_by(pseudo_patientid) %>%
  filter(n() > 1) %>%
  arrange(pseudo_patientid, surgery_date) %>%
  mutate(days_between = as.integer(surgery_date - lag(surgery_date))) %>%
  filter(!is.na(days_between), days_between == 0) %>%
  select(pseudo_patientid, EPISTART, EPIORDER, surgery_type, 
         opcs_primary, all_og_opcs, surgery_date)

# Anchor to NCRAS; earliest elective episode within treatment window
surgery_anchor <- ncras_og %>%
  select(pseudo_patientid, diagmdy, stage_clean) %>%
  left_join(
    hes_og_surgery_episodes %>% filter(!emergency),
    by = "pseudo_patientid"
  ) %>%
  mutate(days_dx_to_surg = as.integer(surgery_date - diagmdy)) %>%
  filter(!is.na(days_dx_to_surg),
         days_dx_to_surg >= -30,
         days_dx_to_surg <= tx_window_days) %>%
  # Curative intent flag: stage 4 + partial gastrectomy = NOT curative
  mutate(
    curative_surgery = !(stage_clean == "4" &
                           surgery_type == "partial_gastrectomy"),
    surgery_class = case_when(
      surgery_type == "oesophagectomy"                        ~ "oesophagectomy",
      surgery_type %in% c("total_gastrectomy",
                          "partial_gastrectomy")              ~ "gastrectomy",
      TRUE ~ NA_character_
    )
  ) %>%
  arrange(pseudo_patientid, surgery_date, EPIORDER) %>%
  distinct(pseudo_patientid, .keep_all = TRUE) %>%
#names(surgery_anchor)  
  select(pseudo_patientid, surgery_date, surgery_type, surgery_class,
         opcs_primary, all_og_opcs, PROCODE3, SITETRET,
         days_dx_to_surg, curative_surgery, "ADMIMETH", "emergency")
names(surgery_anchor)

saveRDS(surgery_anchor, "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/og_surgery_anchor_2015_2022.rds")

cat("Patients with elective OG surgery (any intent):",
    n_distinct(surgery_anchor$pseudo_patientid), "\n")
cat("  of which curative intent:",
    sum(surgery_anchor$curative_surgery, na.rm = TRUE), "\n")

cat("  Oesophagectomy:  ", sum(surgery_anchor$surgery_class == "oesophagectomy",  na.rm = TRUE), "\n")
cat("  Gastrectomy:     ", sum(surgery_anchor$surgery_class == "gastrectomy",      na.rm = TRUE), "\n")
cat("  Curative intent: ", sum(surgery_anchor$curative_surgery,                    na.rm = TRUE), "\n")

# Surgery type breakdown
count(surgery_anchor, surgery_type, surgery_class, curative_surgery)

# =============================================================================
# Check: Examine diagnosis to treatment waiting times
# Note that only <10% have surgery alone, so these are expected to be long.
# =============================================================================

og_cohort <- ncras_og %>%
  
  # --- Diagnostic endoscopy --------------------------------------------------
left_join(endoscopy_anchor %>%
            select(pseudo_patientid, endoscopy_date, days_endo_to_dx),
          by = "pseudo_patientid") %>%
  
  # --- EMR/ESD ---------------------------------------------------------------
left_join(emresd_anchor, by = "pseudo_patientid") %>%
  
  # --- Surgery ---------------------------------------------------------------
left_join(surgery_anchor %>%
            select(pseudo_patientid, surgery_date, surgery_type,
                   surgery_class, opcs_primary, PROCODE3, SITETRET, 
                   ADMIMETH , emergency,
                   days_dx_to_surg, curative_surgery),
          by = "pseudo_patientid") %>%
  
  # --- Waiting time variables ------------------------------------------------
mutate(
  # Primary: diagnosis to surgery
  wt_dx_to_surg = as.integer(surgery_date - diagmdy),
  
  # PI3 equivalent: endoscopy to surgery
  wt_endo_to_surg = as.integer(surgery_date - endoscopy_date)
) %>%
  
  # --- Restrict to surgical patients for now --------------------------------
filter(!is.na(surgery_date))

cat("Surgical cohort:", nrow(og_cohort), "\n")

saveRDS(
  og_cohort,
  "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/og_cohort_ncras_hes_2015_2022.rds"
)

test <- readRDS("E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/og_cohort_ncras_hes_2015_2022.rds")
names(test)
# --- Waiting time summaries -------------------------------------------------
og_cohort %>%
  summarise(
    n              = n(),
    median_dx_surg = median(wt_dx_to_surg, na.rm = TRUE),
    p25_dx_surg    = quantile(wt_dx_to_surg, 0.25, na.rm = TRUE),
    p75_dx_surg    = quantile(wt_dx_to_surg, 0.75, na.rm = TRUE),
    pct_over_62    = round(100 * mean(wt_dx_to_surg > 62,  na.rm = TRUE), 1),
    pct_over_104   = round(100 * mean(wt_dx_to_surg > 104, na.rm = TRUE), 1),
    median_endo_surg = median(wt_endo_to_surg, na.rm = TRUE),
    p25_endo_surg    = quantile(wt_endo_to_surg, 0.25, na.rm = TRUE),
    p75_endo_surg    = quantile(wt_endo_to_surg, 0.75, na.rm = TRUE)
  ) %>%
  print()

# --- By surgery type --------------------------------------------------------
og_cohort %>%
  group_by(surgery_class) %>%
  summarise(
    n              = n(),
    median_dx_surg = median(wt_dx_to_surg, na.rm = TRUE),
    p25            = quantile(wt_dx_to_surg, 0.25, na.rm = TRUE),
    p75            = quantile(wt_dx_to_surg, 0.75, na.rm = TRUE),
    pct_over_62    = round(100 * mean(wt_dx_to_surg > 62, na.rm = TRUE), 1),
    .groups        = "drop"
  ) %>%
  print()

# --- By cancer subtype ------------------------------------------------------
og_cohort %>%
  group_by(cancer_subtype) %>%
  summarise(
    n              = n(),
    median_dx_surg = median(wt_dx_to_surg, na.rm = TRUE),
    p25            = quantile(wt_dx_to_surg, 0.25, na.rm = TRUE),
    p75            = quantile(wt_dx_to_surg, 0.75, na.rm = TRUE),
    .groups        = "drop"
  ) %>%
  print()

# --- By year ----------------------------------------------------------------
og_cohort %>%
  group_by(ydiag) %>%
  summarise(
    n              = n(),
    median_dx_surg = median(wt_dx_to_surg, na.rm = TRUE),
    p25            = quantile(wt_dx_to_surg, 0.25, na.rm = TRUE),
    p75            = quantile(wt_dx_to_surg, 0.75, na.rm = TRUE),
    pct_over_62    = round(100 * mean(wt_dx_to_surg > 62, na.rm = TRUE), 1),
    .groups        = "drop"
  ) %>%
  print()

# --- Negative or implausible waiting times ----------------------------------
cat("\nNegative wt_dx_to_surg:\n")
cat(sum(og_cohort$wt_dx_to_surg < 0, na.rm = TRUE), "\n")

cat("\nDistribution of negative values:\n")
og_cohort %>%
  filter(wt_dx_to_surg < 0) %>%
  pull(wt_dx_to_surg) %>%
  summary()

