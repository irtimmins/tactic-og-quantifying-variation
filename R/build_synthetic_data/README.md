# Synthetic OG cancer waiting-times data

Code that builds a fully synthetic version of the OG (oesophago-gastric) cancer
waiting-times datasets, driven by a disclosure-safe aggregate profile of the
real ICON cohort. The synthetic data lets the downstream merge / provider
analysis be developed and tested off-server, then run unchanged on the real data.

## Layout

```
Data/
  ICON/                                  real ICON inputs (read-only source)
    og_cohort_cwt_2015_2022.rds          post-merge cohort  (needed by script 02)
    ncras_og_2015_2022.rds               }
    OG_endoscopy_anchor_combined.rds     }
    OG_emresd_anchor.rds                 } anchor objects, only needed by
    og_surgery_anchor_2015_2022.rds      } script 01's optional cohort rebuild
    og_sact_anchor_2015_2022.rds         }
    rt_anchor_og.rds                     }
  synthetic/                             everything the build produces
    og_pipeline_spec.rds                 spec object (script 01)
    og_profile_for_synthetic.rds         aggregate profile (script 02)
    og_cohort_precwt_SYNTH.rds / .dta    synthetic Table A (script 03)
    cwt_records_SYNTH.rds / .dta         synthetic Table B (script 03)
    provider_level/
      OG_pairwise_distance_matrix_SYNTH.dta   (script 05)
      NHSHospitals_valid_sites_SYNTH.dta      (script 05)
      NHSHospitals_services_SYNTH.xlsx        (script 06)
R/
  build_synthetic_data/
    01_pre_cwt_specification.R
    02_profile_for_synthetic_data.R
    03_generate_synthetic_data.R
    04_validate_synthetic_dataset.R
    05_create_synthetic_distance_matrix.R
    06_create_synthetic_provider_excel.R
Stata/
  01_data_format_checks.do
```

Paths are resolved with `here::here()`, so scripts run from anywhere inside the
repository (the project root is found via the `.git` / `.Rproj` marker). No paths
need editing between machines.

## Run order

| # | Script | Reads | Writes |
|---|--------|-------|--------|
| 1 | 01_pre_cwt_specification.R | ICON anchors (optional) | og_pipeline_spec.rds |
| 2 | 02_profile_for_synthetic_data.R | spec + og_cohort_cwt_2015_2022.rds | og_profile_for_synthetic.rds |
| 3 | 03_generate_synthetic_data.R | profile | og_cohort_precwt_SYNTH, cwt_records_SYNTH (.rds + .dta) |
| 4 | 04_validate_synthetic_dataset.R | spec + synthetic cohort | (console report) |
| 5 | 05_create_synthetic_distance_matrix.R | synthetic cohort | distance matrix + valid sites |
| 6 | 06_create_synthetic_provider_excel.R | valid sites | provider xlsx |

Scripts 01 and 02 read the real cohort and so run on-server; only the two
artefacts they emit (the spec and the aggregate profile) need to leave the
server. Scripts 03-06 use only those artefacts and run anywhere.

Script 01 always writes the spec (it has no data dependency). It only rebuilds
the real reference cohort if the full ICON anchor set is present in `Data/ICON`,
so the pipeline still runs with just `og_cohort_cwt_2015_2022.rds` in place.

## Stata format check

`Stata/01_data_format_checks.do` confirms the two synthetic `.dta` files are in
the shape the OG merge and provider scripts expect (variable presence, storage
types, parseable CWT dates, C15x/C16x site codes, merge-key overlap). Run it
after script 03:

```stata
global syn "Data/synthetic"
do "Stata/01_data_format_checks.do"
```

## Dependencies

R packages: tidyverse, lubridate, haven, writexl, arrow, here.
