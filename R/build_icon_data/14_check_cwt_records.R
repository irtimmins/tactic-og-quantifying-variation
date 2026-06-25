# =============================================================================
# 14  CWT surgery field check - raw fields + specialist-centre cross-check
# -----------------------------------------------------------------------------
# Second of three CWT-surgery diagnostics (13 -> 14 -> 15), read-only. Run 13
# first to see the group; this script explains what the CWT record actually is.
#
# Goes back to the raw CWT partition for the patients flagged in 13 (No treatment
# + CWT surgery code + no HES resection) and shows EVERY raw CWT field for those
# records. The merge (08) keeps only modality and a few dates; the raw partition
# holds the rest (referral source, priority, treatment-period detail, treating
# org, care setting), which is where the trigger shows up. It then cross-checks
# the treating organisation against the NOGCA Appendix 12 specialist OG centres -
# major resections happen only at those centres, so a high non-specialist share
# is strong evidence the coded "surgery" is a non-resection procedure.
#
# Reads : Data/ICON/og_cohort_cwt_2015_2022.rds, og_surgery_anchor.rds,
#         the partitioned CWT dataset
# Writes: nothing (prints a report). Optional raw dump if write_cwt_trigger_raw.
# =============================================================================

source("R/build_icon_data/01_define_parameters.R")

og   <- readRDS(f_cohort_cwt)
surg <- readRDS(f_surgery_anchor)

# the flagged group, by pseudo_patientid
flagged_ids <- og %>%
  filter(tx_pathway == "No treatment recorded", !is.na(cwt_treat_date),
         cwt_modality %in% c("01","23","24"),
         !pseudo_patientid %in% surg$pseudo_patientid) %>%
  pull(pseudo_patientid) %>% unique()

cat("flagged patients:", length(flagged_ids), "\n\n")

# pull the FULL raw CWT rows for these patients (all columns, surgery modalities)
cwt_raw <- open_dataset(path_cwt_partition) %>%
  filter(site_icd10 %in% og_icd10,
         pseudo_patientid %in% flagged_ids,
         modality %in% c("01","23","24")) %>%
  collect()

cat("raw CWT surgery rows for these patients:", nrow(cwt_raw), "\n")
cat("\n=== all columns available in the raw CWT record ===\n")
print(names(cwt_raw))

# show the distribution of every low-cardinality field - these reveal the trigger
cat("\n=== value distributions of categorical CWT fields ===\n")
for (col in names(cwt_raw)) {
  v <- cwt_raw[[col]]
  # only tabulate fields with a manageable number of distinct values
  if (is.character(v) || is.factor(v)) {
    u <- length(unique(v))
    if (u > 1 && u <= 40) {
      cat("\n--", col, "(", u, "distinct ) --\n")
      print(sort(table(v), decreasing = TRUE))
    } else if (u == 1) {
      cat("\n--", col, ": constant =", unique(v)[1], "\n")
    } else {
      cat("\n--", col, ":", u, "distinct values (too many to tabulate)\n")
    }
  }
}

# a handful of full example rows to eyeball the raw record end to end
cat("\n=== 10 example raw rows (transposed for readability) ===\n")
ex <- cwt_raw %>% slice_head(n = 10)
print(t(as.matrix(ex)))

# -----------------------------------------------------------------------------
# Specialist-centre cross-check (NOGCA State of the Nation 2025, Appendix 12)
# -----------------------------------------------------------------------------
# Major OG resections in England are performed only at the designated specialist
# surgical centres. If the CWT treating organisation for these records is mostly
# NOT one of those centres, the coded "surgery" cannot be a major resection -
# strong, near-definitional evidence the record is a non-resection procedure.
og_specialist_centres <- c(
  "R0D","RA2","RA7","RAE","RAJ","REM","RF4","RGT","RHM","RHQ","RHU","RJ1","RJE",
  "RK9","RKB","RM1","RM3","RPY","RR8","RRK","RRV","RTD","RTE","RTG","RTH","RTR",
  "RWA","RWE","RX1","RXN","RYJ","RYR")

# the treating-org field on the CWT record (first 3 chars = trust). Try the
# known CWT org column names; fall back gracefully if absent.
org_col <- intersect(c("org_treat_start","org_dec_to_treat","org_code_treat",
                       "treatment_org"), names(cwt_raw))[1]
if (!is.na(org_col)) {
  trust3 <- substr(as.character(cwt_raw[[org_col]]), 1, 3)
  is_specialist <- trust3 %in% og_specialist_centres
  cat("\n=== treating org vs NOGCA Appendix 12 specialist OG centres (using ",
      org_col, ") ===\n", sep = "")
  cat("records at a specialist OG surgical centre:    ", sum(is_specialist, na.rm = TRUE), "\n")
  cat("records NOT at a specialist centre:            ", sum(!is_specialist, na.rm = TRUE), "\n")
  cat(sprintf("specialist-centre share: %.1f%%\n",
              100 * mean(is_specialist, na.rm = TRUE)))
  cat("\n-- treating trusts NOT on the specialist list (top 15) --\n")
  print(head(sort(table(trust3[!is_specialist]), decreasing = TRUE), 15))
  cat("\nReading: a high non-specialist share means these CWT 'surgery' records\n")
  cat("are very unlikely to be major resections (which only occur at the centres\n")
  cat("above). The specialist-centre subgroup that remains is examined in\n")
  cat("15_check_cwt_residual.R, which inspects the actual HES episode contents.\n")
} else {
  cat("\n(no recognised treating-org column found in the raw CWT record;\n")
  cat(" available columns are listed above - adjust org_col to match)\n")
}

# optional line-list of the raw records for offline inspection
if (!exists("write_cwt_trigger_raw")) write_cwt_trigger_raw <- FALSE
if (isTRUE(write_cwt_trigger_raw)) {
  out <- file.path(dir_icon, "check_cwt_trigger_raw.rds")
  saveRDS(cwt_raw, out)
  cat("\nraw records saved ->", out, "\n")
}