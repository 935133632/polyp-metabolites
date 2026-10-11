# Analysis log: polyp x NMR metabolites (UK Biobank)

Replicates the Olink-proteomics polyp project (K31 / K63) with 251 NMR metabolites. Inputs are read-only
(`D:\research\UKB_olink\Data`); the original project folder is not modified. Only aggregate results are tracked here.
**Multiple testing: Bonferroni throughout** (P < 0.05 / number of tests; 251 metabolites per scan, so P < 1.99e-4). Reported `p_bonf` = min(P x number of tests, 1).

## Step 0 - data audit and cohort construction (`code/00_build_cohort.py`)
- Metabolome (`09.UKB Metabolome Clean Data.Rdata`) uses the **WH eid space**; the old project's ICD-10/smoking/alcohol files (`HW-*`) use a different, non-linkable eid space (only chance overlap, birth-year agreement 3.5%). All inputs here are therefore the `WH-*` files.
- Participants: 274,260 with NMR metabolites; 251 metabolites.
- Metabolites: rank-based inverse-normal transform, then per-metabolite z-score within complete cases (HR per 1 SD).
- Outcomes (first ICD-10 date): **K63.5 colon polyp**, **K31.7 stomach/duodenal polyp**. Sensitivity definitions `K63any`/`K31any` = original 3-character codes used in the Olink project (K31 is mostly non-polyp stomach disease, so the old "K31 polyp" label was broader than polyps).
- Prevalent cases (disease on/before baseline) excluded; follow-up to first event, death, or 2022-10-31.
- Covariates (Model 1): age, sex, TDI, BMI, SBP, DBP (automated, manual fallback), smoking, drinking; "Prefer not to answer" -> missing.
- Counts: K63.5 N=272,310 / 15,677 events; K31.7 N=273,634 / 7,213 events; K63any 18,443; K31any 10,130.

## Step 1 - primary Cox, Model 1 (`code/01_primary_cox.py`, `results/step01_primary/`)
Per-metabolite Cox: `Surv(days, event) ~ metabolite_z + age + sex + TDI + BMI + SBP + DBP + smoking + drinking`; Bonferroni across 251 metabolites.
- K63.5 colon polyp: 157/251 significant (63 higher, 94 lower risk). Top: Linopt (linoleic acid %) HR 0.89; CE_VSVLDL_pt HR 0.90; TG_VSVLDL_pt HR 1.11.
- K31.7 stomach/duodenal polyp: 174/251 significant (68 higher, 106 lower). Top: CE_VLHDL HR 0.84; CE_LHDL HR 0.83; GlycoAcet HR 1.18.
- 127 metabolites are significant in both outcomes, all in the same direction. Effect correlation between outcomes r = 0.81.
- The significant metabolites are mostly lipoprotein-subclass lipid content and composition variables (highly correlated). The metabolite correlation matrix implies about 19-38 effectively independent tests, so Bonferroni over 251 is conservative.

## Step 2 - sensitivity, excluding events in first 2 years (`code/02_secondary_cox.py main`, `results/step02_sensitivity_exclude_first2y/`)
K63.5: 14,044 events, 142 significant; 141 of the 157 main hits remain, all same direction. K31.7: 6,609 events, 166 significant; 166 of 174 remain, all same direction.

## Step 3 - subgroups (`results/step03_subgroup/`)
All 251 metabolites re-analysed within female / male / age<60 / age>=60 (sex or age removed from covariates in its own stratum). The age cut-off of 58 inherited from the Olink project (the cohort median) was replaced by 60 (conventional); 65 is reported as a supplementary cut-off (`*_subgroups_age65.csv`).
Bonferroni-significant per stratum (female / male / <60 / >=60): K63.5 102 / 120 / 114 / 121; K31.7 151 / 132 / 148 / 151. Supplementary (<65 / >=65): K63.5 143 / 65; K31.7 164 / 101.
Direction agreement with the main result among main hits (female / male / <60 / >=60): K63.5 0.98 / 1.00 / 0.98 / 1.00; K31.7 1.00 / 1.00 / 0.91 / 1.00; at the 65 cut-off 0.99-1.00 for both. Significant in all four sex and age-60 strata: 51 (K63.5), 75 (K31.7).

## Step 3b - baseline prevalent companion (`code/03_prevalent_logistic.py`, `results/step03b_prevalent/`)
Adjusted logistic regression per metabolite for prevalent disease at baseline (same covariates). Significant: K63.5 75/251, K31.7 46/251.

## Step 4 - Bonferroni summary (`code/08_bonferroni_summary.py`)
Class breakdown and consistency counts for the Bonferroni-significant metabolites.

## Step 5 - prediction (`code/05_prediction.py`, `code/05b_coxnet_tuned.py`, `code/07_prediction_summary.py`; `results/step11_13_prediction/`, `tables/Table_*prediction*.csv`, `figures/Fig7_prediction.*`)
Three feature sets: covariates only; covariates + metabolites passing Bonferroni in a univariate screen run **inside the training data** (on average 214 metabolites for K63.5, 167 for K31.7); covariates + all 251. Binary: 5-fold CV (logistic, elastic net, LightGBM, XGBoost; AUC). Survival: 70/30 split (penalised Cox with penalty chosen on an inner validation split, XGBoost Cox); C-index and 5/10-year AUC.
- Colon polyp: AUC 0.626-0.627 (covariates) -> 0.630-0.631 (+metabolites); C-index 0.630 -> 0.636; time-dependent AUC (tuned Cox) 5 y 0.676 -> 0.682, 10 y 0.648 -> 0.656.
- Stomach/duodenal polyp: AUC 0.627-0.632 -> 0.636-0.647; C-index 0.634-0.635 -> 0.647-0.651; time-dependent AUC 5 y 0.656 -> 0.681, 10 y 0.650 -> 0.667.
- Fig 7 shows ROC curves (binary: out-of-fold predictions of the best model; survival: time-dependent ROC at 5 and 10 years). Bar charts of all models are in supplementary Fig S7.
- Bonferroni-selected and all-251 sets perform alike. Absolute discrimination is low (<= 0.65).
- Random forest and extra trees were dropped for run time; random survival forest and DeepSurv were not run (packages unavailable).

## Step 6 - internal validation (`results/step06_internal_validation/`)
Split-sample (70/30, all 251 in both halves): of metabolites significant in the discovery half (K63.5 134, K31.7 149), 100% have the same direction in the validation half, 95% (127, 142) are nominally significant (P<0.05) there and 41% (55) / 77% (114) also pass Bonferroni there; the validation half has only 30% of the sample. Bootstrap (30 resamples, all 251): among the main hits, P<0.05 in >=95% of resamples for 156/157 (K63.5) and 168/174 (K31.7). These show stability within one cohort, not independent replication.

## Step 7 - dual-site overlap (`code/04_dual_site.py`, `results/step07_dual_site/`)
Participants eligible for both incident analyses: none 250,137; K63 only 14,497; K31 only 5,975; both 1,110. All 251 metabolites compared between dual-site and single-site cases (logistic, Bonferroni across 251): one metabolite is significant, Pho_LHDL_pt (OR 1.14 per SD, 95% CI 1.06-1.22). Treat as hypothesis-generating (1,110 cases).

## Step 8 - figures and tables (`code/06_figures_tables.py`, `code/07_prediction_summary.py`, `code/metab_annot.py`, `figures/`, `tables/`)
Fig1 volcano (labels repelled; 5 leading metabolites per direction); Fig2 counts of Bonferroni-significant metabolites by class (rule-based annotation in `metab_annot.py`); Fig3A/B circular bar plot of all 251 metabolites, one figure per outcome (`code/circular_barplot.py`): bars start at HR = 1 (red higher, blue lower risk), black lines 95% CI, stars Bonferroni-adjusted P, grey sectors = groups (lipoprotein particles, total lipids, phospholipids, cholesterol, esterified and free cholesterol, triglycerides, composition %, sizes/apolipoproteins/other, fatty acids, amino acids/glycolysis/ketone bodies/fluid balance); Fig4 cross-outcome and prevalent-vs-incident effect scatter; Fig5 sensitivity + subgroup (sex, age 60) heatmap (20 leading metabolites shown for display; all 251 analysed); Fig6A dual-site group sizes and Fig6B circular bar plot of all 251 metabolites (OR per 1 SD, both polyps vs one; red OR>1, blue OR<1); Fig7 ROC curves. Main tables: cohort characteristics, top-20 per outcome, full primary results (Table S1), class summary. Circular bar plots use the same grouping for all 251 metabolites (Fig3, Fig6B, FigS1); linear top-N forest plots were dropped.

## Step 9 - manuscript draft (`manuscript/polyp_metabolites_draft.md`)
Methods, Results and limitations drafted from the tables above; Introduction and Discussion not drafted.

## Step 10 - supplementary material (`code/09_supplement.py`, `tables/supplement/`, `figures/supplement/`)
`Supplementary_Tables.xlsx` (index sheet + 13 tables, also as CSV): S1 primary Cox; S2 sensitivity; S3 subgroups (sex, age 60); S4 age-65 subgroups; S5 prevalent logistic; S6 split-sample; S7 bootstrap; S8 dual-site; S9 original 3-character definitions; S10 binary prediction; S11 survival prediction; S12 metabolite annotation; S13 cohort characteristics.
Supplementary figures: S1 sensitivity analysis (circular); S2 correlation structure and effective number of tests; S3 polyp-specific vs original 3-character definitions; S4 split-sample concordance and bootstrap stability; S5 age cut-offs 60 and 65; S6 prevalent vs incident (both outcomes); S7 prediction bar charts; S8 dual-site groups and timing; S1A/B circular bar plots of the sensitivity analysis excluding events in the first 2 years (all 251 metabolites).

Note on Fig3 grouping: metabolite abbreviations were assigned to groups by rule. `TC_Chy_VLDL` is treated as triglycerides in the chylomicron/XXL-VLDL subclass (it is the seventh item of the `*_Chy_VLDL` series); `TPR`, `GluLac` and `SpeAla` are placed under 'Sizes, apolipoproteins & other'. Check these against the data dictionary if the figure is used in a paper.
