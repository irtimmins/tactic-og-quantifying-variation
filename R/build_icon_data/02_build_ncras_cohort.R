# =============================================================================
# 02  Build the NCRAS cohort
# -----------------------------------------------------------------------------
# Reads the NCRAS registry extract and the COSD linkage, derives the tumour and
# patient fields, and applies the NOGCA inclusion / exclusion criteria to produce
# the analysis cohort (one row per patient, stage 1-3, 2015+). Writes the cohort
# and the patient-id list that every patient-restricted extract downstream uses.
#
# The NCRAS read is the only ultra-raw read here and is gated by refresh_raw: it
# re-reads the parquet only when the extract is missing or a refresh is asked for.
#
# Reads : NCRAS parquet (gated), COSD dta
# Writes: Data/ICON/ncras_og_2015_2022.rds
# =============================================================================

source("R/build_icon_data/01_define_parameters.R")

# -----------------------------------------------------------------------------
# Read NCRAS (gated): the registry rows for OG sites, with the derived fields
# -----------------------------------------------------------------------------
f_ncras_raw <- file.path(dir_icon, "ncras_og_raw.rds")

if (refresh_raw || !file.exists(f_ncras_raw)) {
  reason <- if (!file.exists(f_ncras_raw)) "extract not found - first run"
  else "refresh_raw = TRUE - rebuilding"
  message(sprintf("Reading NCRAS parquet from raw (%s); this is slow...", reason))
  ncras_og_raw <- open_dataset(path_ncras_parquet) %>%
    filter(sitestr %in% og_icd10) %>%
    collect() %>%
    select(any_of(ncras_cols)) %>%
    mutate(
      pseudo_patientid = as.character(pseudo_patientid),
      diagmdy          = as.Date(diagmdy),
      finmdy           = as.Date(finmdy),
      morphology_num   = as.integer(as.character(typestr))
    ) %>%
    mutate(
      tumour_site_grp = case_when(
        str_starts(as.character(sitestr), "C15") ~ "oesophageal",
        str_starts(as.character(sitestr), "C16") ~ "gastric",
        TRUE ~ NA_character_),
      cancer_subtype = case_when(
        tumour_site_grp == "oesophageal" & morphology_num %in% morph_oes_scc ~ "Oes SCC",
        tumour_site_grp == "oesophageal" & morphology_num %in% morph_oes_aca ~ "Oes ACA",
        tumour_site_grp == "gastric"                                         ~ "Gast",
        TRUE                                                                 ~ NA_character_),
      # stage (NOGCA / AJCC v8): 0 -> 1; X/U/blank -> NA
      stage_clean = case_when(
        as.character(stage_best) == "0"           ~ "1",
        str_starts(as.character(stage_best), "1") ~ "1",
        str_starts(as.character(stage_best), "2") ~ "2",
        str_starts(as.character(stage_best), "3") ~ "3",
        str_starts(as.character(stage_best), "4") ~ "4",
        TRUE                                      ~ NA_character_),
      # route: final_route takes precedence over route_bjc
      final_route_chr = na_if(as.character(final_route), ""),
      route_bjc_chr   = na_if(as.character(route_bjc),   ""),
      route_combined  = factor(coalesce(final_route_chr, route_bjc_chr, "Unknown")),
      emergency_admission = as.integer(
        as.character(route_combined) == "Emergency presentation"),
      surv_from_dx_days = as.integer(finmdy - diagmdy),
      died = as.integer(dead)
    )
  saveRDS(ncras_og_raw, f_ncras_raw)
} else {
  ncras_og_raw <- readRDS(f_ncras_raw)
}

# -----------------------------------------------------------------------------
# Apply inclusion / exclusion (NOGCA section 3)
# -----------------------------------------------------------------------------
ncras_og <- ncras_og_raw %>%
  filter(ydiag >= 2015) %>%                                  # diagnosis 2015+
  filter(agediag >= 18 | is.na(agediag)) %>%                 # adults
  filter(!is.na(morphology_num),                             # histological dx
         morphology_num >= 8001, morphology_num <= 9989) %>%
  filter(morphology_num %in% morph_epithelial) %>%           # epithelial (App 4)
  filter(!morphology_num %in% morph_neuroendocrine) %>%      # exclude NE (App 5)
  filter(as.integer(basisofdiagnosis) != 9L,                 # exclude DCO
         !(died == 1L & diagmdy == finmdy)) %>%
  arrange(pseudo_patientid, diagmdy) %>%                     # earliest primary
  distinct(pseudo_patientid, .keep_all = TRUE) %>%
  filter(stage_clean %in% c("1", "2", "3"))                 # stage 1-3 only

cat("Cohort after inclusion/exclusion:", n_distinct(ncras_og$pseudo_patientid),
    "patients\n")

# -----------------------------------------------------------------------------
# COSD linkage: performance status and CNS involvement
# -----------------------------------------------------------------------------
ncras_og_ids <- ncras_og %>% distinct(pseudo_patientid) %>% pull()

cosd_og <- read_dta(path_cosd_dta) %>%
  filter(pseudo_patientid %in% ncras_og_ids) %>%
  mutate(
    ps_num = case_when(
      performancestatus %in% c("0","1","2","3","4") ~ as.integer(performancestatus),
      TRUE                                          ~ NA_integer_),
    cnsinvolved = case_when(
      str_starts(clinicalnursespecialist, "Y") ~ 1L,
      clinicalnursespecialist == "NN"          ~ 0L,
      TRUE                                     ~ NA_integer_)
  ) %>%
  select(pseudo_patientid, pseudo_tumourid, ps_num, cnsinvolved)

ncras_og <- ncras_og %>%
  left_join(cosd_og, by = c("pseudo_patientid", "pseudo_tumourid"))

saveRDS(ncras_og, f_ncras_cohort)
cat("Saved", f_ncras_cohort, "(", nrow(ncras_og), "patients ).",
    "Next: 03_extract_raw_sources.R\n")

# ---- optional checks (uncomment to inspect) ---------------------------------
# print(count(ncras_og, tumour_site_grp, stage_clean))
# print(count(ncras_og, cancer_subtype, sort = TRUE))
# ncras_og %>% summarise(pct_ps = round(100*mean(!is.na(ps_num)),1),
#                        pct_cns = round(100*mean(!is.na(cnsinvolved)),1)) %>% print()