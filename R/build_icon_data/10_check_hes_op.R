# =============================================================================
# 10  Check HES OP  (optional, print-only)
# -----------------------------------------------------------------------------
# An exploratory check, not part of the build: it does not derive or save
# anything. It confirms the HES outpatient endoscopy supplement behaves as
# expected - in particular the ATTENDED filter (5 = seen, 6 = seen having
# arrived late) used when building the OP endoscopy anchor in script 04.
#
# Reads the saved HES OP extract (no slow raw read), so it is cheap to run. Left
# out of 00_master by default; source it directly when you want to inspect OP
# endoscopy coverage.
#
# Reads: Data/ICON/hes_op_og_2014_2024.rds, ncras cohort
# =============================================================================

source("R/build_icon_data/01_define_parameters.R")

hes_op   <- readRDS(f_hes_op_extract)
ncras_og <- readRDS(f_ncras_cohort)

cat("HES OP extract:", nrow(hes_op), "rows,",
    n_distinct(hes_op$STUDY_ID), "patients\n")

# endoscopy OPCS records in OP, with their ATTENDED status. The build keeps
# ATTENDED 5/6 (attended); this shows how many endoscopy-coded rows sit under
# each status, so the filter can be sanity-checked.
op_op_cols <- names(hes_op)[str_starts(names(hes_op), "OPERTN_")]
op_endoscopy <- hes_op %>%
  pivot_longer(all_of(op_op_cols), names_to = "op_position", values_to = "opcs_code") %>%
  filter(!is.na(opcs_code), opcs_code != "-") %>%
  mutate(opcs4 = normalise_opcs(opcs_code)) %>%
  filter(opcs4 %in% opcs_diagnostic_endoscopy)

cat("\nOP endoscopy-coded rows by ATTENDED status",
    "(5 = attended, 6 = attended late; others not counted in the build):\n")
op_endoscopy %>% count(ATTENDED, sort = TRUE) %>% print()

cat("\nOP endoscopy records (attended only):",
    nrow(filter(op_endoscopy, ATTENDED %in% c("5","6"))),
    "across",
    n_distinct(filter(op_endoscopy, ATTENDED %in% c("5","6"))$STUDY_ID),
    "patients\n")

cat("\n10 check complete.\n")
