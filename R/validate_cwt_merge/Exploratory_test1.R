#######

library(haven)
library(dplyr)

base_dir <- "E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/"

test_data <- read_dta("D:/Projects/#2045_ICON_TACTIC/Project1_interim_bowel/Code_for_OG/synthetic_data_OG_dx_tx.dta")
names(test_data)
print(test_data, n = 10)

og_data <- readRDS("E:/Data_PHE/Extracts/#2045_ICON_TACTIC/Derived/og_cohort_cwt_2015_2022.rds")
write_dta(og_data, paste0(base_dir,  "og_cohort_cwt_2015_2022.dta"))

og_data
names(og_data)
sort(summary(as.factor(og_cohort$diag_hosp[og_cohort$received_curative_tx == T])))
length(unique(og_cohort$diag_hosp))
length(unique(og_cohort$diag_trust))
length(unique(og_cohort$diag_hosp[og_cohort$received_curative_tx == T]))
sort(summary(as.factor(og_cohort$diag_hosp[og_cohort$received_curative_tx == T])))

syn_cohort <- readRDS(paste0(base_dir, "og_cohort_precwt_SYNTH.rds"))
syn_cwt <-   readRDS(paste0(base_dir, "cwt_records_SYNTH.rds"))
syn_cohort
syn_cwt

syn_cohort_out <- readRDS(paste0(base_dir, "og_cohort_precwt_SYNTH.rds"))
syn_cohort_out <- syn_cohort_out %>% mutate(across(where(is.factor), as.character))
write_dta(syn_cohort_out, paste0(base_dir, "og_cohort_precwt_SYNTH.dta"))
write_dta(syn_cwt, paste0(base_dir,  "cwt_records_SYNTH.dta"))

#syn_cohort_out









