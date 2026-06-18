# =============================================================================
# 09  Validation of build logic
# -----------------------------------------------------------------------------
# Tests that the build CODE is correct, independent of the real data. It drives
# controlled fixtures - patients whose correct classification is known by
# construction - through the actual build scripts, and asserts the output is
# exactly what the logic should produce. Every assertion is pass/fail; the
# script stops with a non-zero count if any fail. Nothing here touches Data/ICON.
#
# Two parts:
#   A. plumbing - a handful of patients through the full chain (04 -> 08), to
#      confirm the anchor derivations, the joins and the merge run end to end.
#   B. classification - anchors built directly for one patient per pathway plus
#      the awkward edge cases (the HES-chemo guard, the concurrency boundary, the
#      curative-intent flag), through 07 and 08, asserting tx_pathway,
#      first_tx_date, tx_trust and the audit category exactly.
#
# Run from the project root:  Rscript R/build_icon_data/09_validation_of_build_logic.R
# =============================================================================

suppressPackageStartupMessages(library(tidyverse))
dir_build <- "R/build_icon_data"

# This script repoints dir_icon at temp folders to run on fixtures. Save the
# caller's session state and restore it before returning, so running the test
# never leaves dir_icon (or the read_cwt / refresh_raw seams) pointing at a
# fixture - otherwise a later build in the same session would read the temp
# extracts. restore_session() is called on every exit path (success or stop).
.saved <- mget(c("dir_icon","refresh_raw","read_cwt"),
               ifnotfound = list(NULL), envir = globalenv())
restore_session <- function() {
  for (nm in names(.saved)) {
    if (is.null(.saved[[nm]])) suppressWarnings(rm(list = nm, envir = globalenv()))
    else assign(nm, .saved[[nm]], envir = globalenv())
  }
}

# --- tiny test harness -------------------------------------------------------
.tests <- new.env(); .tests$rows <- list()
expect <- function(label, cond) {
  ok <- isTRUE(cond)
  .tests$rows[[length(.tests$rows) + 1]] <- list(label = label, ok = ok)
  cat(if (ok) "  pass  " else "  FAIL  ", label, "\n")
}
fmt <- function(d) format(d, "%d/%m/%Y")

# build the six anchor / cohort files a classification run needs, in dir_icon,
# from a per-patient spec of treatment offsets (days from diagnosis; NA = absent)
build_fixture_inputs <- function(spec, dx = as.Date("2018-06-01")) {
  off <- function(x) if_else(is.na(x), as.Date(NA), dx + x)
  ncras <- spec %>% transmute(
    pseudo_patientid = id, pseudo_tumourid = paste0("t_", id),
    diagmdy = dx, ydiag = 2018L, sex = "1", agediag = 68L,
    ethnicity_group_broad = "White", NHSE_reversed_imd_quintile_lsoas = "3",
    tumour_site_grp = site, cancer_subtype = subtype, stage_clean = stage,
    route_combined = "Unknown", diag_trust = "R01", diag_hosp = "R01A",
    finmdy = dx + 500L, died = 0L)
  saveRDS(ncras, f_ncras_cohort)
  
  saveRDS(spec %>% transmute(pseudo_patientid = id,
                             endoscopy_date = dx - 5L, days_endo_to_dx = 5L),
          f_endoscopy_anchor)
  
  saveRDS(spec %>% filter(!is.na(emresd)) %>%
            transmute(pseudo_patientid = id, emresd_date = off(emresd),
                      days_dx_to_emresd = emresd), f_emresd_anchor)
  
  saveRDS(spec %>% filter(!is.na(surgery)) %>%
            transmute(pseudo_patientid = id, surgery_date = off(surgery),
                      surgery_type = "oesophagectomy", surgery_class = "oesophagectomy",
                      opcs_primary = "G021", PROCODE3 = "R01", SITETRET = "R01A",
                      days_dx_to_surg = surgery, curative_surgery = cur_surg),
          f_surgery_anchor)
  
  saveRDS(spec %>% filter(!is.na(sact) | !is.na(hes_chemo)) %>%
            transmute(pseudo_patientid = id,
                      chemo_date = off(coalesce(sact, hes_chemo)),
                      days_dx_to_chemo = coalesce(sact, hes_chemo),
                      chemo_source = chemo_src, hes_chemo_date = off(hes_chemo),
                      BENCHMARK_GROUP = "FOLFOX", benchmark_group_lwr = "folfox",
                      INTENT_OF_TREATMENT_V3 = "1", CHEMO_RADIATION = "Y",
                      ORGANISATION_CODE_OF_PROVIDER = "R09"),  # never the trust
          f_chemo_anchor)
  
  saveRDS(spec %>% filter(!is.na(rt)) %>%
            transmute(pseudo_patientid = id, rt_date = off(rt),
                      rt_curative = rt_cur, rt_dose = 50, rt_fractions = 25L,
                      days_dx_to_rt = rt, ORGCODEPROVIDER = "R02"), f_rt_anchor)
  
  saveRDS(spec %>% transmute(pseudo_patientid = id, rcs_ch_score = 0L,
                             cci_any = 0L, cci_group = factor("0", levels = c("0","1","2","3+")),
                             cci_n_conditions = 0L, cci_conditions = "none"), f_cci)
}

# a CWT fixture: one record per treated patient, modality matched to pathway,
# DTT shortly before and treatment on the earliest treatment date
make_cwt <- function(spec, dx = as.Date("2018-06-01")) {
  treated <- spec %>% rowwise() %>%
    mutate(first_off = suppressWarnings(min(c(emresd, surgery, sact, rt, hes_chemo),
                                            na.rm = TRUE))) %>%
    ungroup() %>% filter(is.finite(first_off))
  modality <- function(p) dplyr::case_when(
    p %in% c("EMR/ESD only","EMR/ESD then surgery","Surgery only","Surgery + other",
             "Surgery + adjuvant chemo") ~ "01",
    p %in% c("Surgery + neoadjuvant chemo","Palliative chemo + RT","SACT only") ~ "02",
    p %in% c("Surgery + neoadjuvant chemoRT","Definitive chemoRT") ~ "04",
    p %in% c("Surgery + neoadjuvant RT","Curative RT only","Palliative RT only") ~ "05",
    TRUE ~ "07")
  treated %>% transmute(
    pseudo_patientid = id, site_icd10 = "C159", modality = modality(expect_pathway),
    treat_period_start = fmt(dx + first_off - 10L),
    treat_start = fmt(dx + first_off),
    crtp_date = fmt(dx - 20L), date_first_seen = fmt(dx - 15L), mdt_date = fmt(dx - 8L))
}

# =============================================================================
# Part A: plumbing - full chain on raw fixtures (04 -> 08)
# =============================================================================
cat("\n==========  A. plumbing: full chain (04 -> 08)  ==========\n")
tmpA <- file.path(tempdir(), "og_logicA"); dir.create(tmpA, showWarnings = FALSE)
dir_icon <- tmpA; refresh_raw <- FALSE
source(file.path(dir_build, "01_define_parameters.R"))

# safety: never write fixtures to the real data folder. If the path did not move
# to the temp dir, stop before any saveRDS clobbers a real extract.
if (!startsWith(normalizePath(f_hes_apc_extract, mustWork = FALSE),
                normalizePath(tmpA, mustWork = FALSE)))
  stop("fixture paths did not redirect to the temp folder - aborting to protect ",
       "Data/ICON. (dir_icon = ", dir_icon, ")", call. = FALSE)

dxA <- as.Date("2018-06-01")
saveRDS(tibble(pseudo_patientid = paste0("p", 1:5), pseudo_tumourid = paste0("t", 1:5),
               diagmdy = dxA, ydiag = 2018L, sex = "1", agediag = 68L,
               ethnicity_group_broad = "White", NHSE_reversed_imd_quintile_lsoas = "3",
               tumour_site_grp = c("oesophageal","oesophageal","oesophageal","gastric","gastric"),
               cancer_subtype = c("Oes ACA","Oes ACA","Oes SCC","Gast","Gast"),
               stage_clean = c("2","3","2","1","3"), route_combined = "Unknown",
               diag_trust = "R01", diag_hosp = "R01A", finmdy = dxA + 500L, died = 0L), f_ncras_cohort)

mk_apc <- function(id, opcs, opdate, diag2 = "-") {
  row <- tibble(STUDY_ID = id, ADMIDATE = opdate, ADMIMETH = "11",
                PROCODE3 = "R01", SITETRET = "R01A", EPISTART = opdate,
                EPIORDER = 1L, EPITYPE = "1")
  for (i in 1:24) { row[[sprintf("OPERTN_%02d", i)]] <- if (i == 1) opcs else "-"
  row[[sprintf("OPDATE_%02d", i)]] <- if (i == 1) opdate else as.Date(NA) }
  for (i in 1:20) row[[sprintf("DIAG_4_%02d", i)]] <- if (i == 2) diag2 else "-"
  row
}
saveRDS(bind_rows(mk_apc("p1","G021",dxA+60L,"E119"), mk_apc("p2","G021",dxA+110L),
                  mk_apc("p4","G121",dxA+20L),
                  mk_apc("p1","G459",dxA-3L), mk_apc("p2","G459",dxA-2L),
                  mk_apc("p3","G459",dxA-4L)),   # diagnostic endoscopies pre-dx
        f_hes_apc_extract)
hes_op <- tibble(STUDY_ID = character(), APPTDATE = as.Date(character()),
                 ATTENDED = character(), appt_date = as.Date(character()))
for (i in 1:24) hes_op[[sprintf("OPERTN_%02d", i)]] <- character()
saveRDS(hes_op, f_hes_op_extract)
saveRDS(tibble(pseudo_patientid = c("p2","p3"), sact_regimen_date = c(dxA+20L, dxA+30L),
               sact_cycle_date = c(dxA+20L, dxA+30L), BENCHMARK_GROUP = "FOLFOX",
               benchmark_group_lwr = "folfox", INTENT_OF_TREATMENT_V3 = "1", CYCLE_NUMBER = "1",
               cycle_number = 1L, ORGANISATION_CODE_OF_PROVIDER = "R01", CHEMO_RADIATION = "Y"),
        f_sact_extract)
saveRDS(tibble(pseudo_patientid = "p3", rt_start_date = dxA+32L, rt_decision_date = dxA+25L,
               rt_dose = 50, rt_fractions = 25L, ORGCODEPROVIDER = "R02"), f_rtds_extract)
cwtA <- tibble(pseudo_patientid = c("p1","p2","p3","p4"), site_icd10 = "C159",
               modality = c("01","02","04","01"),
               treat_period_start = fmt(c(dxA+50L,dxA+15L,dxA+25L,dxA+12L)),
               treat_start = fmt(c(dxA+60L,dxA+20L,dxA+30L,dxA+20L)),
               crtp_date = fmt(rep(dxA-10L,4)), date_first_seen = fmt(rep(dxA-5L,4)),
               mdt_date = fmt(rep(dxA,4)))
read_cwt <- function() cwtA

for (s in c("04_derive_hes_treatments.R","05_derive_comorbidities.R",
            "06_derive_sact_rtds.R","07_build_pathways.R","08_merge_cwt.R"))
  source(file.path(dir_build, s), local = FALSE)

ogA <- readRDS(f_cohort_cwt)
pa <- setNames(ogA$tx_pathway, ogA$pseudo_patientid)
expect("chain completes with 5 patients", nrow(ogA) == 5)
expect("p1 surgery only",            pa[["p1"]] == "Surgery only")
expect("p2 neoadjuvant chemo",       pa[["p2"]] == "Surgery + neoadjuvant chemo")
expect("p3 definitive chemoRT",      pa[["p3"]] == "Definitive chemoRT")
expect("p4 EMR/ESD only",            pa[["p4"]] == "EMR/ESD only")
expect("p5 no treatment",            pa[["p5"]] == "No treatment recorded")
expect("all audit/merge columns present",
       length(setdiff(c("tx_modality_audit","received_any_tx","dtt_valid",
                        "received_curative_tx_audit","wt_dx_to_dtt"), names(ogA))) == 0)
expect("p1 comorbidity picked up from HES (diabetes)",
       ogA$cci_any[ogA$pseudo_patientid == "p1"] == 1)

# =============================================================================
# Part B: classification - one patient per pathway + edge cases (07 -> 08)
# =============================================================================
cat("\n==========  B. classification: every pathway + edge cases  ==========\n")
tmpB <- file.path(tempdir(), "og_logicB"); dir.create(tmpB, showWarnings = FALSE)
dir_icon <- tmpB
source(file.path(dir_build, "01_define_parameters.R"))   # recompute f_ paths at tmpB

# spec: offsets in days from diagnosis (NA = that treatment absent)
spec <- tribble(
  ~id,    ~expect_pathway,                 ~expect_first, ~expect_trust, ~stage, ~site, ~subtype,
  ~emresd, ~surgery, ~cur_surg, ~sact, ~rt, ~rt_cur, ~chemo_src, ~hes_chemo,
  "emr",  "EMR/ESD only",                  "emresd","R01","1","oesophageal","Oes ACA",  20, NA, NA,    NA, NA, NA,    NA,    NA,
  "emrs", "EMR/ESD then surgery",          "emresd","R01","1","oesophageal","Oes ACA",  20, 90, TRUE,  NA, NA, NA,    NA,    NA,
  "ncrt", "Surgery + neoadjuvant chemoRT", "min",   "R02","3","oesophageal","Oes ACA",  NA, 120,TRUE,  20, 22, TRUE,  "sact",NA,
  "nc",   "Surgery + neoadjuvant chemo",   "sact",  "R01","3","oesophageal","Oes ACA",  NA, 120,TRUE,  20, NA, NA,    "sact",NA,
  "nrt",  "Surgery + neoadjuvant RT",      "rt",    "R02","3","oesophageal","Oes SCC",  NA, 120,TRUE,  NA, 20, TRUE,  NA,    NA,
  "adj",  "Surgery + adjuvant chemo",      "surgery","R01","2","gastric",   "Gast",     NA, 60, TRUE,  100,NA, NA,    "sact",NA,
  "so",   "Surgery only",                  "surgery","R01","2","gastric",   "Gast",     NA, 60, TRUE,  NA, NA, NA,    NA,    NA,
  "soth", "Surgery + other",               "surgery","R01","2","oesophageal","Oes ACA", NA, 60, TRUE,  60, NA, NA,    "sact",NA,
  "dcrt", "Definitive chemoRT",            "min",   "R02","2","oesophageal","Oes SCC",  NA, NA, NA,    30, 32, TRUE,  "sact",NA,
  "crt",  "Curative RT only",              "rt",    "R02","2","oesophageal","Oes SCC",  NA, NA, NA,    NA, 30, TRUE,  NA,    NA,
  "pcrt", "Palliative chemo + RT",         NA,      NA,   "3","gastric",   "Gast",      NA, NA, NA,    30, 40, FALSE, "sact",NA,
  "sacto","SACT only",                     NA,      NA,   "3","gastric",   "Gast",      NA, NA, NA,    30, NA, NA,    "sact",NA,
  "prt",  "Palliative RT only",            NA,      NA,   "3","oesophageal","Oes ACA",  NA, NA, NA,    NA, 30, FALSE, NA,    NA,
  "none", "No treatment recorded",         NA,      NA,   "2","gastric",   "Gast",      NA, NA, NA,    NA, NA, NA,    NA,    NA,
  # edge cases for the HES-chemo guard and the concurrency boundary
  "hesN", "Definitive chemoRT",            "min",   "R02","2","oesophageal","Oes SCC",  NA, NA, NA,    NA, 30, TRUE,  "hes", 35,  # HES chemo within 28d of RT -> counts
  "hesF", "Curative RT only",              "rt",    "R02","2","oesophageal","Oes SCC",  NA, NA, NA,    NA, 30, TRUE,  "hes", 100, # HES chemo far from RT -> does not
  "conc", "Definitive chemoRT",            "min",   "R02","2","oesophageal","Oes SCC",  NA, NA, NA,    30, 44, TRUE,  "sact",NA   # sact-rt gap exactly 14d -> concurrent
)

build_fixture_inputs(spec)
read_cwt <- function() make_cwt(spec)
source(file.path(dir_build, "07_build_pathways.R"), local = FALSE)
source(file.path(dir_build, "08_merge_cwt.R"),      local = FALSE)

og <- readRDS(f_cohort_cwt)
dx <- as.Date("2018-06-01")
row_of <- function(pid) og[og$pseudo_patientid == pid, ]

# expected first_tx_date by code -> offset
first_off <- function(s) dplyr::case_when(
  s$expect_first == "emresd"  ~ s$emresd,
  s$expect_first == "surgery" ~ s$surgery,
  s$expect_first == "sact"    ~ s$sact,
  s$expect_first == "rt"      ~ s$rt,
  s$expect_first == "min"     ~ pmin(s$sact, s$rt, na.rm = TRUE),
  TRUE                        ~ NA_real_)

for (i in seq_len(nrow(spec))) {
  s <- spec[i, ]; r <- row_of(s$id)
  expect(paste0(s$id, ": pathway = ", s$expect_pathway),
         nrow(r) == 1 && r$tx_pathway == s$expect_pathway)
  # first_tx_date
  fo <- first_off(s)
  if (is.na(s$expect_first))
    expect(paste0(s$id, ": first_tx_date is NA"), is.na(r$first_tx_date))
  else
    expect(paste0(s$id, ": first_tx_date on the right event"),
           !is.na(r$first_tx_date) && r$first_tx_date == dx + fo)
  # tx_trust
  if (is.na(s$expect_trust))
    expect(paste0(s$id, ": tx_trust NA (non-curative)"), is.na(r$tx_trust))
  else
    expect(paste0(s$id, ": tx_trust = ", s$expect_trust), identical(r$tx_trust, s$expect_trust))
}

# audit category mapping follows tx_pathway
expect("audit: surgery-only -> 'Surgery only'",
       row_of("so")$tx_modality_audit == "Surgery only")
expect("audit: neoadjuvant -> 'Surgery plus SACT/RT'",
       row_of("nc")$tx_modality_audit == "Surgery plus SACT/RT")
expect("audit: definitive chemoRT -> 'Definitive chemoRT'",
       row_of("dcrt")$tx_modality_audit == "Definitive chemoRT")
expect("audit: EMR -> 'EMR/ESD'", row_of("emr")$tx_modality_audit == "EMR/ESD")
expect("audit: SACT only -> 'Chemo/RT only (non-curative)'",
       row_of("sacto")$tx_modality_audit == "Chemo/RT only (non-curative)")
expect("audit: none -> 'No treatment recorded'",
       row_of("none")$tx_modality_audit == "No treatment recorded")

# intent + curative flag
expect("intent: definitive chemoRT is Curative",
       row_of("dcrt")$tx_intent_audit == "Curative")
expect("intent: palliative chemo+RT is Non-curative",
       row_of("pcrt")$tx_intent_audit == "Non-curative")
expect("curative flag TRUE for surgery-only", row_of("so")$received_curative_tx_audit)
expect("curative flag FALSE for SACT only", !row_of("sacto")$received_curative_tx_audit)

# dtt_valid is NA for EMR pathways, computed elsewhere
expect("dtt_valid NA for EMR/ESD only", is.na(row_of("emr")$dtt_valid))
expect("dtt_valid not NA for surgery-only", !is.na(row_of("so")$dtt_valid))

# the HES-chemo guard, stated as its own assertions
expect("HES guard: chemo within 28d of RT -> Definitive chemoRT",
       row_of("hesN")$tx_pathway == "Definitive chemoRT")
expect("HES guard: chemo far from RT -> Curative RT only",
       row_of("hesF")$tx_pathway == "Curative RT only")
expect("concurrency boundary at 14d is inclusive (-> Definitive chemoRT)",
       row_of("conc")$tx_pathway == "Definitive chemoRT")

# =============================================================================
# Summary
# =============================================================================
oks <- vapply(.tests$rows, function(x) x$ok, logical(1))
cat(sprintf("\n==========  %d passed, %d failed, %d total  ==========\n",
            sum(oks), sum(!oks), length(oks)))
restore_session()
if (any(!oks)) {
  cat("failed:\n")
  for (x in .tests$rows[!oks]) cat("  -", x$label, "\n")
  stop(sum(!oks), " build-logic test(s) failed", call. = FALSE)
}
cat("build logic validated.\n")