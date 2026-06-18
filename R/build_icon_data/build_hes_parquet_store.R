# =============================================================================
# Build the partitioned HES Parquet store  (one-off, programme-wide)
# -----------------------------------------------------------------------------
# Converts the raw HES APC and HES OP sources into year-partitioned Parquet
# stores under Extracts/#2045_ICON_TACTIC/HES/, sized for fast cohort reads.
# Run once for the whole research programme, not per project. After it, point
# 01_define_parameters.R at the new stores and the cohort extracts in 03 read
# partitioned Parquet (column-projected, predicate-pushed) instead of re-parsing
# multi-GB text - seconds rather than minutes.
#
# Design:
#   - partition by YEAR only. The cohort filter is STUDY_ID %in% ids, which no
#     semantic sub-partition speeds up; year gives date pushdown and that is
#     enough. A second semantic level (provider/month) would just make many tiny
#     files, which scan slower.
#   - target ~256 MB Parquet files within each year (max_rows_per_file), so reads
#     parallelise and conversion never holds a whole 5-12 GB year in memory.
#   - STREAM the OP text in row chunks: a 12 GB file is never loaded whole.
#   - keep ALL patients and a generous column superset, so the store serves every
#     project in the programme; per-project, per-tumour slicing stays in 03.
#
# This writes a lot of data and takes a while, but only once. It is safe to
# re-run: it writes to fresh *_parquet folders and does not touch existing data.
# =============================================================================

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(stringr)
  library(readr)
  library(purrr)
})

dir_raw <- "E:/Data_PHE"
hes_out <- file.path(dir_raw, "Extracts/#2045_ICON_TACTIC/HES")
dir.create(hes_out, recursive = TRUE, showWarnings = FALSE)

# ~256 MB target files: tune down if a year still spikes memory on write
ROWS_PER_FILE  <- 2e6      # rows per Parquet file within a year partition
OP_CHUNK_ROWS  <- 2e6      # rows read per streaming chunk from the OP text
COMPRESSION    <- "zstd"   # better ratio than snappy; fast enough to read
MIN_YEAR       <- 2012L    # skip files before this (cohort is 2015-2022 dx; HES
# lookback starts 2014). Widen if a project needs more.

# source locations
apc_dir <- file.path(dir_raw, "Extracts/#2045_ICON_TACTIC/HES/APC")          # parquet
op_dir  <- file.path(dir_raw, "Raw data files received from PHE READ ONLY/HES/OP")  # txt

# -----------------------------------------------------------------------------
# Column supersets - generous, so the store serves the whole programme. Keep
# everything a project might plausibly need; column projection at read time means
# carrying spare columns costs nothing on read, only a little disk.
# -----------------------------------------------------------------------------
apc_cols <- c("STUDY_ID","ADMIDATE","DISDATE","ADMIMETH","DISMETH","PROCODE3",
              "SITETRET","EPISTART","EPIEND","EPIORDER","EPITYPE","EPISTAT",
              "ADMISORC","DISDEST","CLASSPAT","MAINSPEF","TRETSPEF","STARTAGE",
              "SEX","ETHNOS","LSOA11",
              paste0("OPERTN_", str_pad(1:24, 2, pad = "0")),
              paste0("OPDATE_", str_pad(1:24, 2, pad = "0")),
              paste0("DIAG_4_", str_pad(1:20, 2, pad = "0")),
              paste0("DIAG_3_", str_pad(1:20, 2, pad = "0")))

op_cols  <- c("STUDY_ID","APPTDATE","ATTENDED","PROCODET","TRETSPEF","MAINSPEF",
              "STAFFTYP","FIRSTATT","OUTCOME","PRIORITY","SEX","ETHNOS","LSOA11",
              paste0("OPERTN_", str_pad(1:24, 2, pad = "0")),
              paste0("OPDATE_", str_pad(1:24, 2, pad = "0")),
              paste0("DIAG_4_", str_pad(1:12, 2, pad = "0")))

# helper: write a frame into year=YYYY/ partitions, ~ROWS_PER_FILE per file.
# Each call uses a unique basename tag, so "overwrite" only ever replaces a file
# of the same name (a re-run of the same source/chunk), never another partition.
write_year_partitioned <- function(df, out_dir, basename_tag) {
  if (nrow(df) == 0) return(invisible())
  df %>%
    group_by(year) %>%
    write_dataset(
      path = out_dir, format = "parquet", partitioning = "year",
      existing_data_behavior = "overwrite",
      basename_template = paste0(basename_tag, "-{i}.parquet"),
      max_rows_per_file = ROWS_PER_FILE, compression = COMPRESSION)
}

# =============================================================================
# HES APC  -  already Parquet; re-read per source file (each is one year-ish),
# project columns, derive the year, and write into the partitioned store.
# Reading one source file at a time keeps memory bounded.
# =============================================================================
message("Converting HES APC ...")
apc_out   <- file.path(hes_out, "APC_parquet")
apc_files <- list.files(apc_dir, pattern = "FILE", full.names = TRUE)
stopifnot(length(apc_files) > 0)

for (f in apc_files) {
  tag <- str_extract(f, "(?<=HES_APC_)\\d{4}")          # leading year in the file
  if (is.na(tag) || as.integer(tag) < MIN_YEAR) {
    message("  skip APC ", basename(f), " (before ", MIN_YEAR, ")"); next
  }
  message("  APC ", basename(f))
  ds <- open_dataset(f) %>% select(any_of(apc_cols)) %>% collect()
  ds <- ds %>%
    mutate(across(any_of(c("ADMIDATE","DISDATE","EPISTART","EPIEND")), as.Date),
           year = as.integer(coalesce(str_sub(as.character(EPISTART), 1, 4), tag))) %>%
    filter(year >= MIN_YEAR)                              # guard mixed-year files
  write_year_partitioned(ds, apc_out, paste0("apc-", tag))
  rm(ds); gc()
}
message("  APC done -> ", apc_out)

# =============================================================================
# HES OP  -  pipe-delimited text, 5-12 GB per file. STREAM in row chunks so a
# whole year is never resident. read_delim_chunked applies a callback per chunk;
# each chunk is projected, dated, and written straight into the store.
# =============================================================================
message("Converting HES OP ...")
op_out   <- file.path(hes_out, "OP_parquet")
op_files <- list.files(op_dir, pattern = "\\.txt$", full.names = TRUE)
stopifnot(length(op_files) > 0)

for (f in op_files) {
  tag <- str_extract(f, "(?<=HES_OP_)\\d{4}")
  if (is.na(tag) || as.integer(tag) < MIN_YEAR) {
    message("  skip OP ", basename(f), " (before ", MIN_YEAR, ")"); next
  }
  message("  OP ", basename(f), " (streaming)")
  chunk_i <- 0L
  callback <- function(chunk, pos) {
    chunk_i <<- chunk_i + 1L
    ch <- chunk %>% select(any_of(op_cols))
    if (!"APPTDATE" %in% names(ch)) return(invisible())
    ch <- ch %>%
      mutate(appt_date = as.Date(APPTDATE),
             year = as.integer(coalesce(str_sub(as.character(appt_date), 1, 4), tag))) %>%
      filter(year >= MIN_YEAR)
    write_year_partitioned(ch, op_out, paste0("op-", tag, "-c", chunk_i))
  }
  read_delim_chunked(
    f, callback = SideEffectChunkCallback$new(callback),
    delim = "|", chunk_size = OP_CHUNK_ROWS,
    col_types = cols(.default = col_character()), progress = TRUE)
  message("    ", basename(f), ": ", chunk_i, " chunks written")
}
message("  OP done -> ", op_out)

# -----------------------------------------------------------------------------
# Sanity: row counts per year, both stores (cheap - metadata only)
# -----------------------------------------------------------------------------
report <- function(path, label) {
  if (!dir.exists(path)) { message(label, ": not built"); return(invisible()) }
  ds <- open_dataset(path)
  message(label, " total rows: ", format(nrow(ds), big.mark = ","))
  ds %>% count(year) %>% collect() %>% arrange(year) %>% print()
}
report(apc_out, "APC_parquet")
report(op_out,  "OP_parquet")

message("\nDone. Now point 01_define_parameters.R at:\n  ", apc_out, "\n  ", op_out)