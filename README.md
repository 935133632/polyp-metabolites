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
| `code/08_bonferroni_summary.py` | class breakdown and consistency counts for Bonferroni-significant metabolites |
| `code/06_figures_tables.py`, `code/07_prediction_summary.py` | figures and tables |

`code/common.py` holds paths, outcome definitions and the Cox scan; `code/metab_annot.py` the metabolite class annotation.

Multiple testing uses Bonferroni (P < 0.05/251) throughout.

## Model 1 covariates
Age, sex, Townsend deprivation index, BMI, systolic and diastolic blood pressure, smoking status, drinking status.

## Differences from the Olink project
- Outcomes use polyp-specific codes (K63.5, K31.7). The Olink project used the 3-character codes K63 and K31; K31 mostly captures non-polyp stomach disease. 3-character versions are run for the primary scan as `K63any` and `K31any`.
- Follow-up is censored at death.
- Pathway enrichment, druggability and protein-network steps are not applicable to metabolites and are not repeated.
- scikit-survival and PyTorch were unavailable, so random survival forest and DeepSurv are replaced by a penalised Cox model and XGBoost survival.
- The metabolome and the Olink project use different participant-ID spaces; all inputs here come from the `WH-*` files.
