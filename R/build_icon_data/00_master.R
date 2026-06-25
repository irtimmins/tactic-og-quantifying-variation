# =============================================================================
# 00  Master  -  run the OG data build end to end
# -----------------------------------------------------------------------------
# Runs the build scripts in order. Set refresh_raw below to control whether the
# slow ultra-raw reads (script 03, and the NCRAS read in 02) re-pull from the raw
# sources or reuse the saved extracts under Data/ICON.
#
#   refresh_raw = FALSE  reuse the saved extracts; skip the multi-minute reads.
#                        Use this for any rebuild where the source data is
#                        unchanged - it is much faster.
#   refresh_raw = TRUE   re-read NCRAS, HES APC, HES OP, SACT and RTDS from the
#                        raw sources and rewrite the extracts. Use after a new
#                        data drop.
#
# The synthetic-data pipeline (R/build_synthetic_data) is informed by the cohort
# this build produces: run this first, then profile the cohort there.
# =============================================================================

refresh_raw <- FALSE   # set TRUE to re-read the raw sources
dir_icon    <- "Data/ICON"   # authoritative: a build always targets the real folder
suppressWarnings(rm(read_cwt))   # clear any test seam left in the session

dir_build <- "R/build_icon_data"
step <- function(file) {
  message("\n========== ", file, " ==========")
  source(file.path(dir_build, file), local = new.env())
}

# 01 is sourced by each step below; run the build in order
step("02_build_ncras_cohort.R")    # NCRAS + COSD -> cohort
step("03_extract_raw_sources.R")   # HES APC/OP, SACT, RTDS extracts (gated)
step("04_derive_hes_treatments.R") # HES APC/OP -> endoscopy, EMR/ESD, surgery
step("05_derive_comorbidities.R")  # HES APC (+ surgery/EMR dates) -> Charlson / CCI
step("06_derive_sact_rtds.R")      # SACT + RTDS (+ HES chemo) -> chemo, RT anchors
step("07_build_pathways.R")        # assemble -> flags, tx_pathway, audit cats
step("08_merge_cwt.R")             # CWT -> DTT node, waiting times -> final cohort

# validation
step("09_validation_of_build_logic.R")  # logic tests on fixtures - proves the code
step("10_full_validation.R")            # hard assertions on the real cohort - proves the output

# diagnostics - optional, print-only; consult when a check in 10 fails or to
# investigate a specific question. None of these change the cohort. The three CWT
# scripts are a sequence: 13 describes the group, 14 inspects the raw record, 15
# drills into the specialist-centre residual.
 step("11_check_treatment_classification.R")  # leakage, neoadjuvant clock-stop, HES-chemo timing
 step("12_check_hes_op.R")                     # HES-OP endoscopy coverage
 step("13_check_cwt_surgery.R")                # CWT surgery vs HES resection - characterise
 step("14_check_cwt_records.R")               # what the raw CWT record actually is
 step("15_check_cwt_residual.R")              # specialist-centre residual, HES contents

message("\nBuild complete. Final cohort: ", file.path("Data/ICON",
                                                      "og_cohort_cwt_2015_2022.rds"))