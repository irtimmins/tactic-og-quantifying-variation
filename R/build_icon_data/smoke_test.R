# =============================================================================
# Smoke test  (optional, not part of the build)
# -----------------------------------------------------------------------------
# Exercises the derivation chain (04 -> 08) on tiny synthetic fixtures that match
# each source's schema, so the build logic - joins, flags, the pathway ladder,
# the CWT merge - can be checked in seconds without the multi-minute raw reads.
# It is a structural test, not a clinical one: it confirms the scripts run end to
# end, the column contracts line up, and a few hand-built patients land on the
# pathway they should.
#
# Run from the project root:  Rscript R/build_icon_data/smoke_test.R
# It points dir_icon at a temp folder and feeds the CWT step a fixture, so it
# never touches Data/ICON or the real sources.
# =============================================================================

suppressPackageStartupMessages(library(tidyverse))
dir_build <- "R/build_icon_data"

# redirect the build at a temp folder BEFORE any script defines its defaults
tmp <- file.path(tempdir(), "og_smoke")
dir.create(tmp, showWarnings = FALSE, recursive = TRUE)
dir_icon <- tmp          # 01 honours a pre-set dir_icon
refresh_raw <- FALSE     # skip the gated raw reads (fixtures already "exist")
source(file.path(dir_build, "01_define_parameters.R"))

# -----------------------------------------------------------------------------
# Fixtures: five patients with known intended pathways
#   p1 surgery only   p2 neoadjuvant chemo + surgery   p3 definitive chemoRT
#   p4 EMR/ESD only   p5 no treatment
# -----------------------------------------------------------------------------
dx <- as.Date("2018-06-01")
ncras_og <- tibble(
  pseudo_patientid = paste0("p", 1:5), pseudo_tumourid = paste0("t", 1:5),
  diagmdy = dx, ydiag = 2018L, sex = "1", agediag = 68L,
  ethnicity_group_broad = "White", NHSE_reversed_imd_quintile_lsoas = "3",
  tumour_site_grp = c("oesophageal","oesophageal","oesophageal","gastric","gastric"),
  cancer_subtype = c("Oes ACA","Oes ACA","Oes SCC","Gast","Gast"),
  stage_clean = c("2","3","2","1","3"), route_combined = "Unknown",
  diag_trust = "R01", diag_hosp = "R01A", finmdy = dx + 500L, died = 0L)
saveRDS(ncras_og, f_ncras_cohort)

mk_apc <- function(id, opcs, opdate, diag2 = "-") {
  row <- tibble(STUDY_ID = id, ADMIDATE = opdate, ADMIMETH = "11",
                PROCODE3 = "R01", SITETRET = "R01A",
                EPISTART = opdate, EPIORDER = 1L, EPITYPE = "1")
  for (i in 1:24) {
    row[[sprintf("OPERTN_%02d", i)]] <- if (i == 1) opcs else "-"
    row[[sprintf("OPDATE_%02d", i)]] <- if (i == 1) opdate else as.Date(NA)
  }
  for (i in 1:20) row[[sprintf("DIAG_4_%02d", i)]] <- if (i == 2) diag2 else "-"
  row
}
saveRDS(bind_rows(
  mk_apc("p1", "G021", dx + 60L, "E119"),   # oesophagectomy + diabetes
  mk_apc("p2", "G021", dx + 110L),          # oesophagectomy after neoadj chemo
  mk_apc("p4", "G121", dx + 20L)            # EMR
), f_hes_apc_extract)

hes_op <- tibble(STUDY_ID = character(), APPTDATE = as.Date(character()),
                 ATTENDED = character(), appt_date = as.Date(character()))
for (i in 1:24) hes_op[[sprintf("OPERTN_%02d", i)]] <- character()
saveRDS(hes_op, f_hes_op_extract)

saveRDS(tibble(
  pseudo_patientid = c("p2","p3"),
  sact_regimen_date = c(dx + 20L, dx + 30L), sact_cycle_date = c(dx + 20L, dx + 30L),
  BENCHMARK_GROUP = "FOLFOX", benchmark_group_lwr = "folfox",
  INTENT_OF_TREATMENT_V3 = "1", CYCLE_NUMBER = "1", cycle_number = 1L,
  ORGANISATION_CODE_OF_PROVIDER = "R01", CHEMO_RADIATION = "Y"), f_sact_extract)

saveRDS(tibble(
  pseudo_patientid = "p3", rt_start_date = dx + 32L, rt_decision_date = dx + 25L,
  rt_dose = 50, rt_fractions = 25L, ORGCODEPROVIDER = "R02"), f_rtds_extract)

# CWT fixture, fed to 08 via the read_cwt seam (modality matched to pathway)
fmt <- function(d) format(d, "%d/%m/%Y")
cwt_fixture <- tibble(
  pseudo_patientid = c("p1","p2","p3","p4"), site_icd10 = c("C159","C159","C159","C169"),
  modality = c("01","02","04","01"),
  treat_period_start = fmt(c(dx+50L, dx+15L, dx+25L, dx+12L)),
  treat_start        = fmt(c(dx+60L, dx+20L, dx+30L, dx+20L)),
  crtp_date          = fmt(rep(dx-10L, 4)),
  date_first_seen    = fmt(rep(dx-5L, 4)), mdt_date = fmt(rep(dx, 4)))
read_cwt <- function() cwt_fixture

# -----------------------------------------------------------------------------
# Run 04..08 in this environment (paths and read_cwt persist via global sourcing)
# -----------------------------------------------------------------------------
for (s in c("04_derive_hes_treatments.R", "05_derive_comorbidities.R",
            "06_derive_sact_rtds.R", "07_build_pathways.R", "08_merge_cwt.R")) {
  message("---- ", s)
  source(file.path(dir_build, s), local = FALSE)
}

# -----------------------------------------------------------------------------
# Checks
# -----------------------------------------------------------------------------
og <- readRDS(f_cohort_cwt)
cat("\n================  smoke test results  ================\n")
cat("rows:", nrow(og), "(expected 5)\n")

need_cols <- c("tx_pathway","first_tx_date","tx_trust","cwt_dtt_date","dtt_valid",
               "tx_modality_audit","received_any_tx","received_curative_tx_audit",
               "cci_group","wt_dx_to_dtt","wt_dtt_to_tx")
miss <- setdiff(need_cols, names(og))
cat("missing expected columns:",
    if (length(miss)) paste(miss, collapse = ", ") else "none", "\n")

expected <- c(p1 = "Surgery only", p2 = "Surgery + neoadjuvant chemo",
              p3 = "Definitive chemoRT", p4 = "EMR/ESD only",
              p5 = "No treatment recorded")
got <- setNames(og$tx_pathway, og$pseudo_patientid)[names(expected)]
cat("\npathway checks:\n")
for (i in names(expected))
  cat(sprintf("  %-3s %-30s %s\n", i, got[i],
              if (isTRUE(got[i] == expected[i])) "OK" else paste("EXPECTED", expected[i])))
cat("all pathways as expected:", all(got == expected, na.rm = TRUE), "\n")

trust <- setNames(og$tx_trust, og$pseudo_patientid)
cat("\np1 tx_trust == surgery provider R01:", identical(trust[["p1"]], "R01"), "\n")
cat("p3 tx_trust == RT provider R02:",      identical(trust[["p3"]], "R02"), "\n")
cat("p5 tx_trust missing (no treatment):",  is.na(trust[["p5"]]), "\n")
cat("p1 has diabetes comorbidity:",
    og$cci_any[og$pseudo_patientid == "p1"] == 1, "\n")
cat("=====================================================\n")
