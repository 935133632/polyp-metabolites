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
All 251 metabolites re-analysed within female / male / age<60 / age>=60 (sex or age removed from covariates in its own stratum). The age cut-off of 58 inherited from the Olink project (the cohort median) was replaced by 60 (conventional; balanced events: colon 7,107 vs 8,377, stomach/duodenal 2,964 vs 4,157). A cut-off of 65 was also run and gave the same conclusions (direction agreement 99-100%) but was dropped to keep one cut-off.
Bonferroni-significant per stratum (female / male / <60 / >=60): K63.5 102 / 120 / 114 / 121; K31.7 151 / 132 / 148 / 151.
Direction agreement with the main result among main hits (female / male / <60 / >=60): K63.5 0.98 / 1.00 / 0.98 / 1.00; K31.7 1.00 / 1.00 / 0.91 / 1.00. Significant in all four strata: 51 (K63.5), 75 (K31.7).

## Step 3b - baseline prevalent companion (`code/03_prevalent_logistic.py`, `results/step03b_prevalent/`)
Adjusted logistic regression per metabolite for prevalent disease at baseline (same covariates). Significant: K63.5 75/251, K31.7 46/251.

## Step 4 - Bonferroni summary (`code/08_bonferroni_summary.py`)
Class breakdown and consistency counts for the Bonferroni-significant metabolites.

## Step 5 - prediction (`code/05_prediction.py`, `code/05b_coxnet_tuned.py`, `code/07_prediction_summary.py`; `results/step11_13_prediction/`, `tables/Table_*prediction*.csv`, `figures/Fig7_prediction_ROC.*`)
Three feature sets: covariates only; covariates + metabolites passing Bonferroni in a univariate screen run **inside the training data** (on average 214 metabolites for K63.5, 167 for K31.7); covariates + all 251. Binary: 5-fold CV (logistic, elastic net, LightGBM, XGBoost; AUC). Survival: 70/30 split (penalised Cox with penalty chosen on an inner validation split, XGBoost Cox); C-index and 5/10-year AUC.
- Colon polyp: AUC 0.626-0.627 (covariates) -> 0.630-0.631 (+metabolites); C-index 0.630 -> 0.636; time-dependent AUC (tuned Cox) 5 y 0.676 -> 0.682, 10 y 0.648 -> 0.656.
- Stomach/duodenal polyp: AUC 0.627-0.632 -> 0.636-0.647; C-index 0.634-0.635 -> 0.647-0.651; time-dependent AUC 5 y 0.656 -> 0.681, 10 y 0.650 -> 0.667.
- Fig 7 shows ROC curves (binary: out-of-fold predictions of the best model; survival: time-dependent ROC at 5 and 10 years). Bar charts of all models are in supplementary Fig S7.
- Bonferroni-selected and all-251 sets perform alike. Absolute discrimination is low (<= 0.65).
- Random forest and extra trees were dropped for run time; random survival forest and DeepSurv were not run (packages unavailable).

## Step 6 - internal validation (`results/step06_internal_validation/`)
Split-sample (70/30, all 251 in both halves): of metabolites significant in the discovery half (K63.5 134, K31.7 149), 100% have the same direction in the validation half, 95% (127, 142) are nominally significant (P<0.05) there and 41% (55) / 77% (114) also pass Bonferroni there; the validation half has only 30% of the sample. Bootstrap (30 resamples, all 251; summary table also holds the 2.5th/97.5th percentile of the bootstrap coefficient, used in FigS4): among the main hits, P<0.05 in >=95% of resamples for 156/157 (K63.5) and 168/174 (K31.7). These show stability within one cohort, not independent replication.

## Step 7 - dual-site overlap (`code/04_dual_site.py`, `results/step07_dual_site/`)
Participants eligible for both incident analyses: none 250,137; K63 only 14,497; K31 only 5,975; both 1,110. All 251 metabolites compared between dual-site and single-site cases (logistic, Bonferroni across 251): one metabolite is significant, Pho_LHDL_pt (OR 1.14 per SD, 95% CI 1.06-1.22). Treat as hypothesis-generating (1,110 cases).

## Step 11 - extended adjustment (`code/11_extra_covariates.py`, `code/12_model_adjustment.py`, `results/step12_model_adjustment/`, `results/step15_summary/model_adjustment_summary.csv`)
Added waist circumference, body fat %, HbA1c and log CRP (baseline). All models are fitted in the same complete-case sample (K63.5 240,720 participants / 13,822 events; K31.7 241,891 / 6,370). M1 = Model 1; M2a = M1 + waist + body fat %; M2b = M1 + HbA1c + log CRP; M2c = all four.
Main hits (Bonferroni-significant in M1 in this sample: 145 colon, 172 stomach/duodenal) still significant: M2a 117 / 155, M2b 122 / 148, M2c 88 / 139. All keep their direction. Median attenuation of ln HR (colon / stomach/duodenal): M2a 17% / 12%, M2b 12% / 13%, M2c 26% / 22%. Lipid-lowering, aspirin and PPI use are not available in the WH exports and are not adjusted. Fig8A/B circular heat rings.

## Step 12 - specificity (`code/13_specificity.py`, `results/step13_specificity/`, `results/step15_summary/specificity_profile_correlation.csv`)
Model 1 Cox scan of all 251 metabolites for 23 other digestive diseases with >= 2,000 incident events plus colorectal cancer (C18-C20 combined). The vector of 251 ln(HR) values of each outcome was correlated with the colon and stomach/duodenal polyp profiles.
Result: the polyp profile is **not specific**. Correlation with the colon polyp profile is 0.88-0.93 for 12 of 24 outcomes (e.g. diverticular disease 0.93, anorectal diseases 0.93, gastritis 0.92, liver disease incl. fatty liver 0.92, haemorrhoids 0.91, irritable bowel syndrome 0.90); with the stomach/duodenal polyp profile it is 0.93-0.97 for reflux, hiatus hernia, gallstones, peritoneal disorders, diverticular disease and irritable bowel syndrome. Colorectal cancer correlates 0.76 / 0.57 but only 4 metabolites pass Bonferroni. Inguinal hernia is inversely related (-0.88 / -0.69) and umbilical hernia weakly related (0.43 / 0.56), so the profile is not simply shared by every outcome. By metabolite class, glycoprotein acetyls (GlycA) have a positive mean ln(HR) in all outcomes. Fig9.

## Step 13 - metabolite score and dose-response (`code/14_score_dose_response.py`, `results/step14_score_dose_response/`)
Score = metabolite part of a penalised Cox model (penalty 0.0005, l1 ratio 0.1, Model 1 covariates + 251 metabolites) fitted in the 70% training set; evaluated in the independent 30% test set (same split as the prediction analysis; colon 81,693 participants / 4,703 events; stomach/duodenal 82,091 / 2,164).
Colon polyp: HR per SD 1.19 (1.16-1.23); Q4 vs Q1 1.52 (1.39-1.66); P-trend 2e-22; 14-year cumulative incidence about 4.1% (Q1) vs 8.8% (Q4). Stomach/duodenal polyp: HR per SD 1.34 (1.27-1.40); Q4 vs Q1 2.02 (1.76-2.31); P-trend 7e-25; about 1.9% vs 3.9%. Fig10. (A first run failed at the cumulative-incidence step because a column named `T` collided with the DataFrame transpose; the cumulative incidence was recomputed from the saved weights with `... cuminc`, no refit.)
Dose-response: restricted cubic splines (4 knots) for the 3 strongest positive and negative metabolites plus GlycA, HDL-C and total triglycerides. Curves are mostly monotonic; nominal P for non-linearity < 0.05 in 7 of 9 (colon) and 5 of 8 (stomach/duodenal) metabolites (curvature at the extremes, e.g. HDL-C higher risk only at low values); not corrected for multiple testing. FigS8.

## Step 8 - figures and tables (`code/06_figures_tables.py`, `code/07_prediction_summary.py`, `code/metab_annot.py`, `figures/`, `tables/`)
Fig1 volcano (labels repelled; 5 leading metabolites per direction); Fig2 counts of Bonferroni-significant metabolites by class (rule-based annotation in `metab_annot.py`); Fig3A/B circular bar plot of all 251 metabolites, one figure per outcome (`code/circular_barplot.py`): bars start at HR = 1 (red higher, blue lower risk), black lines 95% CI, stars Bonferroni-adjusted P, grey sectors = groups (lipoprotein particles, total lipids, phospholipids, cholesterol, esterified and free cholesterol, triglycerides, composition %, sizes/apolipoproteins/other, fatty acids, amino acids/glycolysis/ketone bodies/fluid balance); Fig4 cross-outcome and prevalent-vs-incident effect scatter; Fig5A/B circular heat rings of all 251 metabolites by stratum (overall, female, male, age <60, >=60), Fig5C sensitivity scatter (excluding events <2 y); Fig6A dual-site group sizes and Fig6B circular bar plot of all 251 metabolites (OR per 1 SD, both polyps vs one; red OR>1, blue OR<1); Fig7 ROC curves. Main tables: cohort characteristics, top-20 per outcome, full primary results (Table S1), class summary. Circular bar plots use the same grouping for all 251 metabolites (Fig3, Fig6B, FigS1); linear top-N forest plots were dropped.

## Step 9 - manuscript draft (`manuscript/polyp_metabolites_draft.md`)
Introduction, Methods, Results and Discussion (with limitations and conclusion) drafted from the tables above. Numbers in the Discussion come from `code/10_discussion_numbers.py` (family-level mean HRs, prevalent-vs-incident agreement). The reference list at the end of the draft marks every entry as found or to be verified: the literature was identified from search-result summaries, not read in full, so each citation and estimate must be checked before submission; two statements still carry `[ref needed]` (adenoma-carcinoma background, and a literature check that no earlier metabolome-wide polyp analysis exists).

## Step 10 - supplementary material (`code/09_supplement.py`, `tables/supplement/`, `figures/supplement/`)
`Supplementary_Tables.xlsx` (index sheet + 5 tables, also as CSV): **S1** metabolite-level results, one row per outcome x metabolite (502 rows, 73 columns): primary Cox, sensitivity excluding events <2 y, prevalent odds ratio, sex and age-60 subgroups, split-sample, bootstrap, extended adjustment (M1, M2a-M2c in the same complete-case sample), original 3-character definition, metabolite-score weight, dual-site odds ratio; **S2** metabolite annotation (class, lipoprotein family); **S3** prediction performance (binary and survival, three feature sets); **S4** specificity profile correlations; **S5** specificity effect estimates (23 digestive outcomes + colorectal cancer). Small tables whose content is given in the text were dropped (cohort characteristics = main Table 1; adjustment summary, score quartiles, spline curves and non-linearity P values are in the Results text and Fig S8).
Supplementary figures: S1A/B sensitivity analysis (circular); S2A circular heat rings of the correlation of each metabolite with principal components 1-8, S2B cumulative variance (effective number of tests); S3 polyp-specific vs original 3-character definitions; S4 split-sample concordance and bootstrap caterpillar plot (main estimate with bootstrap 95% interval, all 251 metabolites); S5 prevalent vs incident; S6 prediction bar charts; S7 dual-site groups and timing; S8 dose-response curves (restricted cubic splines).

Note on Fig3 grouping: metabolite abbreviations were assigned to groups by rule. `TC_Chy_VLDL` is treated as triglycerides in the chylomicron/XXL-VLDL subclass (it is the seventh item of the `*_Chy_VLDL` series); `TPR`, `GluLac` and `SpeAla` are placed under 'Sizes, apolipoproteins & other'. Check these against the data dictionary if the figure is used in a paper.

## Decision log (what changed from the first version and why)
- **Outcome definition**: polyp-specific ICD-10 codes (K63.5, K31.7) instead of the 3-character K63/K31 of the Olink project; 3-character versions kept as `K63any`/`K31any` (Step 0).
- **Multiple testing**: BH-FDR was used first and replaced by Bonferroni everywhere (all result tables carry `p_bonf`, no FDR columns remain).
- **No arbitrary top-N**: subgroup, split-sample, bootstrap, dual-site and prediction analyses use all 251 metabolites; figures show all 251 in circular plots (linear top-N forest plots dropped).
- **Prediction**: features are chosen inside the training data (earlier top-30 selection on the full data leaked information); penalised Cox penalty tuned on an inner validation split.
- **Age cut-off**: 60 instead of 58 (the Olink project's hard-coded cohort median); 65 was checked and dropped.
- **Smoking and drinking** were added as covariates once the WH-eid exports were provided.
- **Figures**: ROC curves for prediction (not bar charts); circular bar/heat plots for all-metabolite results; volcano labels repelled.

## Repository notes
- Git history was rewritten (four times, on request) to remove superseded versions of code, results, figures and the analysis log that used top-N subsets or BH-FDR; the history now contains one commit per analysis step plus later revisions. Commit hashes before 2026-10-11 therefore no longer exist.
- Run logs (`results/*.log`, `*.err`) and bootstrap replicate files (`results/step06_internal_validation/bootstrap_chunks/`) are not tracked; individual-level data (`processed/`) are never tracked.
- **Extended analyses (Steps 11-13)**: added after the main results to address adiposity/glycaemia/inflammation confounding, specificity and effect size at the score level. The specificity result changes the interpretation: the metabolic profile is largely shared with other digestive diseases, so it should be described as a general metabolic-inflammatory risk profile rather than a polyp-specific signature.

## Consistency review of the manuscript draft
Every number in the Results and Discussion was compared with the result files (cohort characteristics, HRs and CIs of the leading metabolites, subgroup agreement, split-sample and bootstrap counts, adjustment, specificity, score, prediction); no numerical errors were found. Wording was aligned after the specificity result: the Introduction now lists all analyses; the Discussion summary, the inflammation paragraph and the site-difference paragraph no longer imply polyp-specific or strongly site-specific signatures; the prediction summary now states the same limits in Results and Discussion (AUC and C-index up to 0.65, time-dependent AUC up to 0.68). Every figure and table cited in the text exists, and every supplementary figure and table (S1-S9, S1-S21) is cited at least once. Circular heat rings: group names were moved closer to the metabolite labels.
- **Supplementary tables consolidated** from 21 to 5: all metabolite-level results in one wide table, small summary tables replaced by text.
