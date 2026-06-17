library(haven)

test_registry <- read_dta(paste0(base_dir, "og_ncras_treatment_synthetic.dta"))
test_cwt <- read_dta(paste0(base_dir, "og_cwt_records_synthetic.dta"))
View(test_registry)
View(test_cwt)
test_all <- og_cwt_merge(A_tableA, syn_cwt)
View(test_all)

og_cwt_merge

########################################















