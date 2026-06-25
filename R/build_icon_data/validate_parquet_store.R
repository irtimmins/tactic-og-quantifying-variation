# =============================================================================
# Validate the HES APC / OP parquet stores
# -----------------------------------------------------------------------------
# Standalone, read-only. Checks the year-partitioned parquet stores that
# build_hes_parquet_store.R produces and that 03_extract_raw_sources.R reads.
# It does NOT rebuild anything. Run it after a (re)build of the stores, or any
# time you want to confirm the parquet layer is intact before a cohort build.
#
# Three tiers, mirroring 10_full_validation.R:
#   structural  - stores exist, open, have rows and the expected columns
#   coverage    - full year span present, no missing/empty year partitions
#   linkage     - the cohort's patients actually find rows in each store
#
# Every check is collected; prints a tier-by-tier result and stops non-zero on
# any failure. Reads paths and the column lists from 01_define_parameters.R.
# =============================================================================

source("R/build_icon_data/01_define_parameters.R")

.checks <- list()
chk <- function(tier, label, cond, detail = "") {
  .checks[[length(.checks) + 1]] <<-
    list(tier = tier, label = label, ok = isTRUE(cond), detail = detail)
}

# columns each store must carry for 03 to build its extracts
apc_required <- c("STUDY_ID","ADMIDATE","ADMIMETH","PROCODE3","SITETRET",
                  "EPISTART","EPIORDER","EPITYPE", op_cols, opdate_cols, "year")
op_required  <- c("STUDY_ID","ATTENDED", op_cols, "year")   # plus a date, checked below

# 03 reads from 2014 onward; the store floor is earlier. Require an unbroken run
# from this year to the most recent full year present.
min_year_needed <- 2014L

# -----------------------------------------------------------------------------
# helper: open a store safely and return the dataset (or NULL)
# -----------------------------------------------------------------------------
open_store <- function(path) {
  if (!dir.exists(path)) return(NULL)
  tryCatch(open_dataset(path), error = function(e) NULL)
}

apc <- open_store(path_hes_apc_dir)
op  <- open_store(path_hes_op_dir)

# =============================================================================
# Structural
# =============================================================================
chk("structural", "APC store exists and opens", !is.null(apc),
    if (is.null(apc)) path_hes_apc_dir else "")
chk("structural", "OP store exists and opens", !is.null(op),
    if (is.null(op)) path_hes_op_dir else "")

if (!is.null(apc)) {
  apc_cols <- names(apc$schema)
  miss <- setdiff(apc_required, apc_cols)
  chk("structural", "APC has all columns 03 needs", length(miss) == 0,
      if (length(miss)) paste("missing:", paste(head(miss, 8), collapse = ", ")) else "")
  chk("structural", "APC has rows", nrow(apc) > 0,
      format(nrow(apc), big.mark = ","))
}

if (!is.null(op)) {
  op_cols_present <- names(op$schema)
  miss <- setdiff(op_required, op_cols_present)
  chk("structural", "OP has all columns 03 needs", length(miss) == 0,
      if (length(miss)) paste("missing:", paste(head(miss, 8), collapse = ", ")) else "")
  # 03 accepts either APPTDATE or appt_date for the attendance date
  has_date <- any(c("APPTDATE","appt_date") %in% op_cols_present)
  chk("structural", "OP has an attendance-date column (APPTDATE/appt_date)", has_date)
  chk("structural", "OP has rows", nrow(op) > 0,
      format(nrow(op), big.mark = ","))
}

# =============================================================================
# Coverage - year partitions present and non-empty
# =============================================================================
check_years <- function(ds, name) {
  yr <- ds %>% count(year) %>% collect() %>% arrange(year)
  present <- yr$year
  most_recent_full <- max(present[present < as.integer(format(Sys.Date(), "%Y"))],
                          na.rm = TRUE)
  needed <- min_year_needed:most_recent_full
  missing_years <- setdiff(needed, present)
  chk("coverage", paste(name, "- unbroken year run", min_year_needed, "to",
                        most_recent_full),
      length(missing_years) == 0,
      if (length(missing_years)) paste("missing:", paste(missing_years, collapse = ", ")) else "")
  # flag any year whose count is a tiny fraction of the median (truncated write)
  med <- median(yr$n)
  thin <- yr %>% filter(year >= min_year_needed, n < med * 0.05)
  chk("coverage", paste(name, "- no suspiciously thin year partition"),
      nrow(thin) == 0,
      if (nrow(thin)) paste(thin$year, collapse = ", ") else "")
}
if (!is.null(apc)) check_years(apc, "APC")
if (!is.null(op))  check_years(op,  "OP")

# =============================================================================
# Linkage - the cohort's patients find rows in each store
# =============================================================================
ncras <- tryCatch(readRDS(f_ncras_cohort), error = function(e) NULL)
if (!is.null(ncras) && !is.null(apc)) {
  ids <- unique(as.character(ncras$pseudo_patientid))
  # sample to keep the scan cheap; a real linkage failure shows up in any sample
  s <- head(ids, 2000)
  apc_hit <- apc %>% filter(STUDY_ID %in% s) %>%
    summarise(p = n_distinct(STUDY_ID)) %>% collect() %>% pull(p)
  chk("linkage", "APC: cohort sample finds matching rows", apc_hit > 0,
      sprintf("%d of %d sampled ids matched", apc_hit, length(s)))
}
if (!is.null(ncras) && !is.null(op)) {
  ids <- unique(as.character(ncras$pseudo_patientid))
  s <- head(ids, 2000)
  op_hit <- op %>% filter(STUDY_ID %in% s) %>%
    summarise(p = n_distinct(STUDY_ID)) %>% collect() %>% pull(p)
  chk("linkage", "OP: cohort sample finds matching rows", op_hit > 0,
      sprintf("%d of %d sampled ids matched", op_hit, length(s)))
}

# =============================================================================
# Report and stop on any failure
# =============================================================================
res <- tibble(tier = map_chr(.checks, "tier"), label = map_chr(.checks, "label"),
              ok = map_lgl(.checks, "ok"), detail = map_chr(.checks, "detail"))
for (t in c("structural","coverage","linkage")) {
  rows <- res %>% filter(tier == t)
  if (nrow(rows) == 0) next
  cat("\n==========  ", t, "  ==========\n")
  rows %>%
    transmute(result = if_else(ok, "pass", "FAIL"), label,
              detail = if_else(detail == "", "", paste0("(", detail, ")"))) %>%
    pwalk(function(result, label, detail)
      cat(sprintf("  %-5s %s %s\n", result, label, detail)))
}
cat(sprintf("\n==========  %d passed, %d failed, %d total  ==========\n",
            sum(res$ok), sum(!res$ok), nrow(res)))
if (any(!res$ok))
  stop(sum(!res$ok), " parquet check(s) failed - see the FAIL lines above. ",
       "Rebuild the affected store with build_hes_parquet_store.R.", call. = FALSE)
cat("parquet stores validated - APC and OP are sound.\n")