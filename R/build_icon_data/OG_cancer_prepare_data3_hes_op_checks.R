library(arrow)
library(haven)
library(dplyr)
library(tidyr)
library(stringr)
library(forcats)
library(purrr)
library(lubridate)
library(readr)


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

ncras_og     <- readRDS("E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/ncras_og_2015_2022.rds")
ncras_og_ids <- ncras_og %>% distinct(pseudo_patientid) %>% pull()
#View(hes_op_raw )
hes_op_raw_test <- map_dfr(
  hes_op_file_list,
  ~{
    read_delim(
      .x,
      delim          = "|",
      col_select     = any_of(op_cols_select),
      col_types      = cols(.default = col_character()),
      show_col_types = FALSE,
      n_max = 1000000
    ) %>%
      filter(
        STUDY_ID %in% ncras_og_ids#,
   #     ATTENDED %in% c("5", "6")
      ) %>%
      mutate(
        STUDY_ID  = as.character(STUDY_ID),
        appt_date = as.Date(APPTDATE)
      )
  },
  .progress = TRUE
)

cat("HES-OP rows (OG cohort, attended):", nrow(hes_op_raw_test), "\n")
cat("Patients:                         ", n_distinct(hes_op_raw_test$STUDY_ID), "\n")

# Helper: normalise OPCS codes for consistent matching
normalise_opcs <- function(x) str_replace_all(str_to_upper(as.character(x)), "\\.", "")


hes_op_raw_test%>%
  pivot_longer(
    cols      = starts_with("OPERTN_"),
    names_to  = "op_position", 
    values_to = "opcs_code"
  ) %>%
  filter(!is.na(opcs_code), opcs_code != "-") %>%
  mutate(opcs4 = normalise_opcs(opcs_code)) %>%
  filter(opcs4 %in% opcs_diagnostic_endoscopy) %>%
  count(ATTENDED, sort = TRUE)

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

