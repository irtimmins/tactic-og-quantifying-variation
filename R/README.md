# OG Cancer Waiting Times - synthetic data workflow

This bundle lets a partner institution develop and test the **CWT merge step** of the
oesophago-gastric (OG) cancer waiting-times pipeline against a realistic synthetic
dataset, then apply the identical merge code to their own real cohort.

The partner has already built the first stage (cancer registry + HES-APC + SACT + RTDS
                                               resolved to one row per patient). What they have not yet done is merge the Cancer
Waiting Times (CWT) data and derive the decision-to-treat (DTT) interval metrics. These
scripts provide a known-good target format, a realistic synthetic stand-in for both
tables, and validation that the synthetic data is a drop-in for the real thing.


## Two environments

The work spans a secure environment and an open one. Only one aggregate file crosses
the boundary.

**Inside the Trusted Research Environment (real patient data):**
  the canonical cohort is built, the real CWT merge is run, and a non-disclosive
statistical profile is extracted.

**Outside (synthetic, shareable):**
  the profile drives a generator that produces a synthetic cohort + CWT pair, which is
validated and used to develop the merge.

The **only** artefact that leaves the secure environment is
`og_profile_for_synthetic.rds` - aggregate distributions with counts rounded to 5 and
cells under 10 suppressed. It still must pass the local output-checking / DAR process
before export.


## Artefacts and run order

Run in this sequence. "TRE" = inside the secure environment; "Open" = anywhere.

1. **(TRE) Existing scripts 1-2** - build `ncras_og` and the treatment anchor objects
(endoscopy, EMR/ESD, surgery, SACT, RTDS). Pre-existing; not part of this bundle.

2. **(TRE) `OG_build_pre_cwt_spec.R`** - assembles the canonical pre-CWT cohort in the
agreed spec format and defines the authoritative data dictionary.
- In: `ncras_og_*.rds` + the anchor objects.
- Out: `og_cohort_precwt_spec_2015_2022.rds` (the cohort) and
`og_pipeline_spec.rds` (the `pre_cwt_spec` / `cwt_spec` manifests, pathway levels,
                        globals). The manifest is non-disclosive and travels with this bundle.
- The single curative clock-stop is `first_tx_date`; the treating trust is `tx_trust`.

3. **(TRE) CWT merge - `OG_cancer_prepare_data4_cwt_merge.R`** - joins the CWT DTT anchor
onto the cohort and derives `wt_endo_to_dtt`, `wt_dtt_to_tx`, `wt_dx_to_dtt`,
`dtt_valid`, plus the validation tables. This is the step the partner is building.
- In: `og_cohort_precwt_spec_*.rds`, the partitioned CWT dataset, `og_pipeline_spec.rds`.
- Out: `og_cohort_cwt_2015_2022.rds` - the analysis-ready cohort.

4. **(TRE) `OG_profile_for_synthetic.R`** - extracts the distributional parameters needed
to generate realistic synthetic data, disclosure-safe.
- In: `og_cohort_cwt_*.rds`, the partitioned CWT dataset.
- Out: `og_profile_for_synthetic.rds` ??? **output-check ??? export**.

5. **(Open) `OG_synthetic_generators.R`** - generates the synthetic pair from the profile.
- In: `og_profile_for_synthetic.rds`.
- Out: `og_cohort_precwt_SYNTH.rds` (Table A, spec-conformant) and
`cwt_records_SYNTH.rds` (Table B, raw `dd/mm/yyyy` dates). Prints a merge-QC block
comparing synthetic behaviour to the real targets.

6. **(Open) `OG_synthetic_validation.R`** - confirms the synthetic cohort conforms to the
spec and is internally consistent (see Validation below).

7. **(Open) Partner develops the merge** - runs the CWT merge (step 3) against the
synthetic pair, iterates on their code, then applies the same merge unchanged to their
own real pre-CWT cohort inside their TRE.


## The merge contract

Two tables join on `pseudo_patientid`, plus two globals (`og_icd10`, `tx_window_days`).

**Table A - pre-CWT cohort** (`og_cohort_precwt_*`): one row per patient. Merge-critical
columns are `pseudo_patientid`, `diagmdy`, `endoscopy_date`, `surgery_date`, `sact_date`,
`rt_date`, `first_tx_date`, `tx_pathway`, `route_combined`,
`NHSE_reversed_imd_quintile_lsoas`, `wt_endo_to_tx`, `wt_dx_to_tx`. The rest is needed to
regenerate pathways/trust/PIs and for realism. Full dictionary in `pre_cwt_spec`.

**Table B - CWT records** (`cwt_records_*`): one or more rows per patient, dates stored as
character `dd/mm/yyyy` so the merge's `as.Date()` parsing runs. Fields: `pseudo_patientid`,
`site_icd10`, `modality`, `crtp_date`, `date_first_seen`, `mdt_date`, `treat_period_start`
(??? DTT), `treat_start`.

The merge keeps the earliest in-window DTT per patient, joins it to Table A, and flags
`dtt_valid` (DTT after diagnosis and at/before treatment; `NA` for EMR/ESD pathways where
DTT is less meaningful).


## How the synthetic data is built (and why the merge behaves realistically)

The two tables are not generated independently. The cohort is generated first -
`tx_pathway` from its marginal, then `stage`, `subtype`, and `year` conditional on
pathway, then treatment dates constructed per pathway so the real build logic re-derives
the same pathway and `first_tx_date`. The CWT record is then **anchored** to the cohort:
`cwt_treat_date = first_tx_date + ??` (?? from the per-pathway `cwt_vs_first_tx`
distribution), and `cwt_dtt_date = cwt_treat_date ??? offset`. That is what reproduces the
~71% exact treatment-date agreement and the EMR/ESD-then-surgery negative-interval quirk
rather than leaving them to chance.


## Validation

Run `OG_synthetic_validation.R`. Expected results:

- **Conformance:** no missing required columns, no type mismatches, all `tx_pathway`
  values within the 14 canonical levels, `stage_clean` in {1,2,3}, zero duplicate IDs.
- **Internal consistency:** `tx_pathway` re-derivation match and `first_tx_date` match
  should both be ~100% - confirming the saved cohort round-trips through the build logic
  on its own.
- **Generator merge-QC** (printed by `OG_synthetic_generators.R`, block 5), latest run:

  | Metric | Real (ICON) | Synthetic |
  |---|---|---|
  | CWT anchor patients | 36,197 | ~35,200 |
  | `days_dx_to_dtt` median | 39 | 39 |
  | `dtt_to_cwt_treat` median | 11 | 11 |
  | `cwt_vs_first_tx` exact match | 71.1% | ~67-70% |
  | within 14 days | 85.6% | ~90% |
  | negative `dtt_to_tx` | 5.3% | ~7% |
  | EMR/ESD-then-surgery negative | ~50% | ~43% |

  Medians match exactly and every rate is within a few points. The residual gap is a
  single artefact: the real `cwt_vs_first_tx` has a sharp spike at 0 plus a fat tail, and
  any quantile-interpolation sampler smears that spike slightly.


## Known limitations / tuning levers

- **`cwt_vs_first_tx` spike.** The exact-match / within-14 rates are the most sensitive
  knobs. A `pct_exact` point mass at 0 pins exact-match; drawing the non-zero branch from
  the outer quantiles only would tighten within-14 further. Optional polish.
- **Rare pathways** (e.g. Surgery + neoadjuvant RT, n???80) fall back to overall interval
  distributions, so their synthetic spread is approximate.
- **stage × subtype** correlation is not jointly preserved - pathway is conditioned on
  each separately. Add a `pathway × subtype × stage` cube to the profiler if that
  secondary association matters.
- **Counts** are rounded to 5 in the profile, so synthetic marginals match the real ones
  to within rounding, not exactly.


## File index

| File | Environment | Purpose |
|---|---|---|
| `OG_build_pre_cwt_spec.R` | TRE | Build canonical pre-CWT cohort + spec manifests |
| `OG_cancer_prepare_data4_cwt_merge.R` | TRE | CWT merge + DTT intervals + validation |
| `OG_profile_for_synthetic.R` | TRE | Disclosure-safe distributional profile |
| `OG_synthetic_generators.R` | Open | Generate synthetic cohort + CWT pair |
| `OG_synthetic_validation.R` | Open | Conformance + internal-consistency checks |
| `og_pipeline_spec.rds` | both | Spec manifests, pathway levels, globals |
| `og_profile_for_synthetic.rds` | exported | Aggregate profile (after output-checking) |
| `og_cohort_precwt_SYNTH.rds` | Open | Synthetic Table A |
| `cwt_records_SYNTH.rds` | Open | Synthetic Table B |
