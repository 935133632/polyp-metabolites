# polyp-metabolites

NMR metabolomic associations with incident colon polyp (ICD-10 K63.5) and stomach/duodenal polyp (K31.7) in the UK Biobank.
This repeats the earlier Olink plasma-proteomics polyp project with 251 Nightingale NMR metabolites.

Only code and aggregate results are tracked. UK Biobank individual-level data (`processed/`) are excluded and cannot be shared.
`ANALYSIS_LOG.md` records each analysis step, the decisions taken and the headline numbers.

## Pipeline (run from the repository root, in order)

| Script | Step |
|---|---|
| `code/00_build_cohort.py` | cohort, covariates, outcome dates, metabolite matrix (rank-based inverse-normal transform) |
| `code/01_primary_cox.py` | per-metabolite Cox, Model 1, Bonferroni |
| `code/02_secondary_cox.py` | sensitivity (exclude events <2 y), subgroups, split-sample validation, bootstrap (all 251 metabolites; `main`, `boot`, `merge` modes) |
| `code/03_prevalent_logistic.py` | baseline prevalent companion analysis |
| `code/04_dual_site.py` | overlap of the two outcomes |
| `code/05_prediction.py`, `code/05b_coxnet_tuned.py` | binary and survival prediction, feature selection inside training data only |
| `code/06_figures_tables.py`, `code/07_prediction_summary.py` | main figures (Fig 1-7) and tables |
| `code/08_bonferroni_summary.py` | class breakdown and consistency counts for Bonferroni-significant metabolites |
| `code/09_supplement.py` | supplementary tables (Excel + CSV) and supplementary figures |
| `code/10_discussion_numbers.py` | family-level effect summaries quoted in the Discussion |
| `code/11_extra_covariates.py`, `code/12_model_adjustment.py` | extended adjustment (waist, body fat %, HbA1c, CRP) in the same complete-case sample |
| `code/13_specificity.py` | Model 1 scan for 23 other digestive outcomes and colorectal cancer |
| `code/14_score_dose_response.py` | metabolite score (train/test) and spline dose-response |
| `code/15_extended_figures.py` | figures and tables for steps 11-13 |
| `code/run_queue.ps1` | throttled launcher for the extended analyses (Windows) |

`code/common.py` holds paths, outcome definitions and the Cox scan; `code/metab_annot.py` the metabolite class annotation; `code/circular_barplot.py` the circular bar and heat-ring plots.

Multiple testing uses Bonferroni (P < 0.05/251) throughout.

## Outputs
- `figures/`: Fig1 volcano; Fig2 class summary; Fig3A/B circular bar plots of all 251 metabolites; Fig4 cross-outcome and prevalent-vs-incident; Fig5A/B circular heat rings by sex/age stratum, Fig5C sensitivity scatter; Fig6A dual-site counts, Fig6B circular OR plot; Fig7 ROC curves; Fig8A/B circular heat rings after extended adjustment; Fig9 specificity against other digestive diseases; Fig10 metabolite score (cumulative incidence and quartile HRs).
- `figures/supplement/`: S1 sensitivity (circular), S2 correlation structure and effective tests, S3 original vs polyp-specific definition, S4 internal validation, S5 age cut-off 65, S6 prevalent vs incident, S7 prediction bars, S8 dual-site timing, S9 dose-response curves.
- `tables/`: main tables; `tables/supplement/Supplementary_Tables.xlsx` (index + 21 tables, also as CSV).
- `results/`: per-step aggregate results; `manuscript/`: Introduction, Methods, Results and Discussion draft (references to be verified).

## Model 1 covariates
Age, sex, Townsend deprivation index, BMI, systolic and diastolic blood pressure, smoking status, drinking status.

## Differences from the Olink project
- Outcomes use polyp-specific codes (K63.5, K31.7). The Olink project used the 3-character codes K63 and K31; K31 mostly captures non-polyp stomach disease. 3-character versions are run for the primary scan as `K63any` and `K31any`.
- Follow-up is censored at death.
- Bonferroni instead of BH-FDR; subgroup age cut-off 60 (65 supplementary) instead of 58; all analyses use all 251 metabolites rather than top-N subsets.
- Pathway enrichment, druggability and protein-network steps are not applicable to metabolites and are not repeated.
- scikit-survival and PyTorch were unavailable, so random survival forest and DeepSurv are replaced by a penalised Cox model and XGBoost survival.
- The metabolome and the Olink project use different participant-ID spaces; all inputs here come from the `WH-*` files.
