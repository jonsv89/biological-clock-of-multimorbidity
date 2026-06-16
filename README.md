# Follow-up duration determines multimorbidity detection in primary care

**Sánchez-Valle J, Zambrana C, Navarro-Martínez A, Costa FX, Rocha L, Cirillo D, Violán C, Valencia A.**

Barcelona Supercomputing Center (BSC)

---

## Overview

This repository contains the analysis pipeline used to study how observation window length determines comorbidity detectability in primary care, using 11 years of EHR data from 5.8 million individuals in Catalonia (SIDIAP, 2008–2018). The pipeline builds matched case-control cohorts for each index disease, computes directed comorbidity associations across nine temporal windows (five cumulative: 0–1 to 0–5 years; four conditional: 1–2 to 4–5 years), applies empirical Bayes shrinkage, and generates all manuscript figures and tables.

---

## Pipeline

Scripts must be run in order. Scripts 01–02 are run per disease pair; scripts 03–05 aggregate results.

### `01_build_pool_A.R` — Build control candidate pool for disease A

For a given index disease A, identifies eligible cases (patients with at least 5 years of follow-up after their first diagnosis of A) and builds a pool of candidate controls matched on sex, age group, and index date (±2 months). Controls must also have ≥5 years of follow-up and cannot have a prior diagnosis of A. The candidate pool is written to disk in FST format (for high-prevalence diseases) or CSV (for low-prevalence diseases), along with a row index for efficient downstream access.

```bash
Rscript 01_build_pool_A.R --disease_a E11 --data_dir ./CaseControlStudy/ --out_dir ./pools/
```

---

### `02_analyze_AB.R` — Match cases to controls and record disease B outcome

For a given disease pair A→B, reads the candidate pool produced by script 01, filters out controls that have a prior diagnosis of B or whose index diagnosis is B, and assigns up to 5 controls per case without replacement (prioritising cases with fewer available candidates). Records whether each case and control develops disease B within the follow-up period, and writes complete and partial matched sets to disk.

```bash
Rscript 02_analyze_AB.R --disease_a E11 --disease_b I10 \
  --data_dir ./CaseControlStudy/ --pool_dir ./pools/ --out_dir ./matched/
```

---

### `03_compute_OR_clogit.R` — Odds ratio via conditional logistic regression

Fits a conditional logistic regression model (using `survival::clogit`) for each combination of temporal window and population stratum (all, women, men). Both cumulative windows (0–k years) and conditional windows (year k−1 to k, restricted to individuals free of B before that interval) are evaluated. Returns OR, 95% CI, SE, and p-value per stratum and window.

```bash
Rscript 03_compute_OR_clogit.R --disease_a E11 --disease_b I10 \
  --matched_dir ./matched/ --out_dir ./or_clogit/
```

---

### `03_compute_RR_contingency.R` — Relative risk and OR via contingency tables

Computes RR and OR from 2×2 contingency tables for the same window × population combinations. RR is estimated following Morris & Gardner (1988); OR uses a continuity-corrected formula. Uses the fully matched (complete, 1:5) dataset. Outputs event counts alongside effect estimates, which are needed for downstream filtering.

```bash
Rscript 03_compute_RR_contingency.R --disease_a E11 --disease_b I10 \
  --matched_dir ./matched/ --counts_dir ./counts/ --out_dir ./rr_ct/
```

---

### `04_shrinkage.R` — Empirical Bayes shrinkage with `ashr`

Aggregates per-pair RR or OR estimates across all disease pairs for a given window and population stratum, and applies adaptive shrinkage (`ashr`) to obtain posterior effect sizes, posterior SEs, and local false sign rates (lfsr). Two pre-screening strategies are supported: `events` (at least one event in cases or controls) and `westergaard` (replicating the criteria from Westergaard et al. 2019). Outputs one shrunk results file per window × population combination.

```bash
Rscript 04_shrinkage.R --results_dir ./Results --method rr_contingency \
  --prescreening events --out_dir ./Results/shrinkage_events
```

---

### `05_manuscript_analyses.R` — All manuscript figures and analyses

Modular script driven by a command-line argument (`args[1]`) that selects the analysis to run. Each module reads from the shrinkage output files and produces figures (PDF) and results tables. Key modules include:

| Argument | Description |
|---|---|
| `prepare_networks` | Filters dagger-asterisk pairs, computes directionality (θ, binomial test, FDR) following Jensen et al. 2014 / Westergaard et al. 2019 |
| `compare_prevalences` | Sex-biased disease prevalence in Catalonia and Denmark; Fisher tests, enrichment by ICD-10 category, UpSet plots |
| `compare_age_of_diagnoses` | Sex differences in age at diagnosis; Welch t-tests, WLS regression by category, comparison with Danish data |
| `calculate_correlations_between_diseaseprevalence_and_numberofcomorbidities` | Correlation between disease prevalence and number of comorbidity partners across time windows |
| `biological_clock` | Two-dimensional biological clock of multimorbidity: cumulative RR slopes, conditional risk persistence, 4-quadrant classification, sex-stratified analysis, paradigmatic examples, directionality within the clock |
| `BiologicalClock_level2` | Biological clock at ICD-10 subcategory level using a WHO 2016 subcategory mapping |
| `Validate_Biological_Clock` | Robustness analyses: pair-level bootstrap (1,000 replicates), leave-one-index-disease-out sensitivity, k-means quadrant comparison |
| `Sex_window_differences` | Sex-specific network overlap and enrichment by disease category pair; directional reversals between women and men |
| `Temporal_robustness_of_comorbidities` | Nine-window detection taxonomy (omnipresent, persistent, late-emerging, transient, conditional-only, mixed); stacked bar and heatmap figures |
| `Temporal_robustness_of_comorbidities_by_pairs` | Same taxonomy but at the catA × catB pair level |
| `compare_westergaard` | Compares the Catalonia network (window 0–5) against the Westergaard et al. 2019 Danish hospital network. Computes network overlap (Jaccard, hypergeometric test), Spearman concordance of RR for shared pairs, enrichment of catA×catB pairs among shared vs population-specific associations, and directional concordance/reversals between the two datasets. Outputs UpSet plots, concordance scatter plots, and enrichment heatmaps |

```bash
Rscript 05_manuscript_analyses.R biological_clock
Rscript 05_manuscript_analyses.R compare_time_windows incremental
```

---

## Dependencies

R packages: `data.table`, `fst`, `survival`, `ashr`, `igraph`, `UpSetR`, `ggplot2`, `ggrepel`, `ggtext`, `ggnewscale`, `patchwork`, `plotly`, `ggalluvial`, `forestplot`, `EbayesThresh`, `VennDiagram`, `dendextend`, `gplots`, `gridExtra`, `MASS`, `flextable`, `officer`, `dplyr`, `tidyr`, `scales`

---

## Input data

The pipeline expects anonymised EHR data from SIDIAP (not publicly available). Required input files:

- `cohort.rds` — one row per patient-diagnosis with fields `idp`, `cod`, `dat`, `sexe`, `rangos`, `followup_end`, `followup_years`
- `valid_diseases_any.rds` — ICD-10 codes eligible as control index diagnoses
- `valid_diseases_A.rds` — ICD-10 codes eligible as index disease A (prevalence < 20%)
- `Data/ICD10_prevalence_Catalonia.txt` — disease prevalence by sex
- `Data/ICD10_three_digits_names.txt` — ICD-10 code → disease name mapping
- `icd10_3digitos_categoria_subcategoria_who2016.csv` — ICD-10 subcategory mapping (for script 05 `BiologicalClock_level2`)
- `Epidemiology/41467_2019_8475_MOESM4_ESM.txt` — Danish prevalence data (Westergaard et al. 2019, Supplementary Data 4)
- `Epidemiology/41467_2019_8475_MOESM5_ESM.txt` — Danish age-of-diagnosis data (Westergaard et al. 2019, Supplementary Data 5)

---
