# =============================================================================
# Build the partitioned HES Parquet store  (one-off, programme-wide)
# -----------------------------------------------------------------------------
# Converts the raw HES APC and HES OP text extracts into year-partitioned
# Parquet stores under Extracts/#2045_ICON_TACTIC/HES/, sized for fast cohort
# reads. Run once for the whole research programme, not per project. After it,
# 01_define_parameters.R points at the new stores and the cohort extracts in 03
# read partitioned Parquet (column-projected, predicate-pushed) in seconds.
#
# Reads straight from the authoritative pipe-delimited .txt (not the pre-made
# APC parquet), so both stores are built the same way and every column is stored
# as character. HES code fields carry non-numeric values (ADMIMETH = "2D"); all-
# character storage keeps the schema consistent across partitions and avoids the
# int-vs-char unification error a later read would otherwise hit.
#
# Design:
#   - partition by YEAR only (derived per record from the date; the cohort filter
#     is STUDY_ID %in% ids, which no semantic sub-partition speeds up).
#   - target ~256 MB files within a year (max_rows_per_file).
#   - STREAM both APC and OP in row chunks: a multi-GB file is never resident.
#   - keep ALL patients and a column superset; per-project slicing stays in 03.
#
# Takes a while, but only once. Safe to re-run: writes fresh *_parquet folders.
# =============================================================================

suppressPackageStartupMessages({
  library(arrow); library(dplyr); library(stringr); library(readr); library(purrr)
})

dir_raw <- "E:/Data_PHE"
hes_in  <- file.path(dir_raw, "Raw data files received from PHE READ ONLY/HES")
hes_out <- file.path(dir_raw, "Extracts/#2045_ICON_TACTIC/HES")
dir.create(hes_out, recursive = TRUE, showWarnings = FALSE)

ROWS_PER_FILE <- 2e6      # rows per Parquet file within a year partition
CHUNK_ROWS    <- 2e6      # rows read per streaming chunk from the .txt
COMPRESSION   <- "zstd"
MIN_YEAR      <- 2012L    # skip records before this (cohort 2015-2022; lookback 2014)

apc_in <- file.path(hes_in, "APC")
op_in  <- file.path(hes_in, "OP")

# only the real data files: HES_<type>_<year>99.txt, never *_rowcount / .dta / .mdmp
data_files <- function(dir, type)
  list.files(dir, pattern = paste0("HES_", type, "_\\d{6}\\.txt$"), full.names = TRUE)

# column supersets (generous, for programme reuse; any_of skips absent ones)
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

# write a frame into year=YYYY/ partitions, ~ROWS_PER_FILE per file. Each call
# uses a unique tag so "overwrite" only replaces a same-named file (a re-run of
# the same chunk), never another partition.
write_year_partitioned <- function(df, out_dir, tag) {
  if (nrow(df) == 0) return(invisible())
  df %>% group_by(year) %>%
    write_dataset(path = out_dir, format = "parquet", partitioning = "year",
                  existing_data_behavior = "overwrite",
                  basename_template = paste0(tag, "-{i}.parquet"),
                  max_rows_per_file = ROWS_PER_FILE, compression = COMPRESSION)
}

# stream one pipe-delimited HES .txt into the store: read in chunks, keep the
# wanted columns as character, derive year from the date column, drop old years.
convert_txt <- function(f, out_dir, keep_cols, date_col, tag) {
  message("  ", basename(f), " (streaming)")
  chunk_i <- 0L
  cb <- function(chunk, pos) {
    chunk_i <<- chunk_i + 1L
    ch <- chunk[, intersect(keep_cols, names(chunk)), drop = FALSE]
    if (!date_col %in% names(ch)) return(invisible())
    d <- as.Date(ch[[date_col]])
    ch$year <- as.integer(coalesce(str_sub(as.character(d), 1, 4), tag))
    ch <- ch[!is.na(ch$year) & ch$year >= MIN_YEAR, , drop = FALSE]
    write_year_partitioned(ch, out_dir, paste0(tag, "-c", chunk_i))
  }
  read_delim_chunked(f, callback = SideEffectChunkCallback$new(cb),
                     delim = "|", chunk_size = CHUNK_ROWS,
                     col_types = cols(.default = col_character()), progress = TRUE)
  message("    ", chunk_i, " chunks")
}

# =============================================================================
# APC
# =============================================================================
message("Converting HES APC from text ...")
apc_out <- file.path(hes_out, "APC_parquet")
apc_files <- data_files(apc_in, "APC")
stopifnot(length(apc_files) > 0)
for (f in apc_files) {
  tag <- str_extract(f, "(?<=HES_APC_)\\d{4}")
  if (is.na(tag) || as.integer(tag) < MIN_YEAR) {
    message("  skip ", basename(f), " (before ", MIN_YEAR, ")"); next }
  convert_txt(f, apc_out, apc_cols, "EPISTART", paste0("apc-", tag))
}
message("  APC done -> ", apc_out)

# =============================================================================
# OP
# =============================================================================
message("Converting HES OP from text ...")
op_out <- file.path(hes_out, "OP_parquet")
op_files <- data_files(op_in, "OP")
stopifnot(length(op_files) > 0)
for (f in op_files) {
  tag <- str_extract(f, "(?<=HES_OP_)\\d{4}")
  if (is.na(tag) || as.integer(tag) < MIN_YEAR) {
    message("  skip ", basename(f), " (before ", MIN_YEAR, ")"); next }
  convert_txt(f, op_out, op_cols, "APPTDATE", paste0("op-", tag))
}
message("  OP done -> ", op_out)

# -----------------------------------------------------------------------------
# Sanity: row counts per year (metadata only, cheap)
# -----------------------------------------------------------------------------
report <- function(path, label) {
  if (!dir.exists(path)) { message(label, ": not built"); return(invisible()) }
  ds <- open_dataset(path)
  message(label, " total rows: ", format(nrow(ds), big.mark = ","))
  ds %>% count(year) %>% collect() %>% arrange(year) %>% print()
}
report(apc_out, "APC_parquet")
report(op_out,  "OP_parquet")
message("\nDone. 01_define_parameters.R already points at these stores.")