# Plasma NMR metabolomics and incident colon and stomach/duodenal polyps in the UK Biobank (draft: Methods and Results)

Draft generated from the analysis record in `ANALYSIS_LOG.md`. Numbers come from `tables/` and `results/`. Introduction and Discussion are not drafted.

## Methods

### Study population and outcomes
We used UK Biobank participants with baseline NMR metabolomic profiling (n = 274,260; 251 metabolite measures). Incident disease was defined from the first hospital inpatient ICD-10 diagnosis date: colon polyp (K63.5) and polyp of the stomach and duodenum (K31.7). Participants with the outcome on or before the baseline date were excluded from the corresponding incident analysis. Follow-up ended at the first event, death or 31 October 2022. The analytic cohorts comprised 272,310 participants with 15,677 incident colon polyps (median follow-up 13.6 years) and 273,634 participants with 7,213 incident stomach/duodenal polyps (median follow-up 13.7 years). As a secondary definition, the original 3-character codes K63 and K31 were also analysed.

### Metabolites and covariates
Each metabolite was rank-based inverse-normal transformed and then standardised within the complete cases of each model, so hazard ratios (HR) are per 1 SD. Models were adjusted for age, sex, Townsend deprivation index, body mass index, systolic and diastolic blood pressure, smoking status and alcohol drinking status; "prefer not to answer" was set to missing.

### Metabolome-wide association analysis
Cox proportional hazards models were fitted for each metabolite and outcome separately. Significance was assessed with Bonferroni correction across the 251 metabolites (P < 0.05/251 = 1.99 × 10⁻⁴) in all analyses. Because the metabolites are highly correlated (about 19–38 effectively independent tests by principal-component and Li-Ji estimates), this correction is conservative, and the number of significant metabolites does not equal the number of independent signals.

### Additional analyses
Sensitivity analysis excluded events within the first 2 years. Subgroup analyses of all metabolites were run by sex and by age (<60 vs ≥60 years; <65 vs ≥65 as a supplementary cut-off). Internal validation used a 70/30 split-sample replication and 30 bootstrap resamples, both over all 251 metabolites. A baseline companion analysis used adjusted logistic regression for prevalent disease. For participants eligible for both outcomes, groups were classified as colon only, stomach/duodenal only or both, and all 251 metabolites were compared between dual-site and single-site cases by logistic regression.

### Prediction
Models used covariates only, covariates plus metabolites passing a Bonferroni-corrected univariate screen run within the training data, or covariates plus all 251 metabolites. Binary prediction used 5-fold stratified cross-validation (logistic regression, elastic net, LightGBM, XGBoost; AUC). Survival prediction used a 70/30 split with a penalised Cox model (penalty chosen on an inner validation split) and an XGBoost Cox model, evaluated by C-index.

## Results

### Cohort
Participants who developed colon polyps were older (58.6 vs 56.4 years), more often male (55.5% vs 45.4%) and had higher BMI (28.4 vs 27.4 kg/m²) than those who did not (Table 1). Those who developed stomach/duodenal polyps were older (59.5 vs 56.5 years) and less often male (39.3% vs 46.2%).

### Metabolome-wide associations
Of 251 metabolites, 157 were associated with incident colon polyp (63 with higher and 94 with lower risk) and 174 with incident stomach/duodenal polyp (68 higher, 106 lower) after Bonferroni correction (Fig. 1, Table S1). Most were lipoprotein-subclass lipid content or composition measures, with both higher and lower risk directions (Fig. 2). For colon polyp, the strongest associations were with the proportion of linoleic acid in total fatty acids (HR 0.89 per SD, 95% CI 0.88–0.91; lower risk), the cholesteryl-ester share of very small VLDL (HR 0.90, 0.88–0.92; lower risk) and the triglyceride share of very small VLDL (HR 1.11, 1.09–1.13; higher risk). For stomach/duodenal polyp, the strongest were cholesteryl esters in very large and large HDL (HR 0.84, 0.81–0.86 and HR 0.83, 0.81–0.85; lower risk) and glycoprotein acetyls (HR 1.18, 1.15–1.21; higher risk) (Fig. 3 shows all 251 metabolites as a circular bar plot for each outcome; Table 2 and Fig. S1 show the 20 leading).

### Consistency across outcomes and analyses
Effect estimates for the two outcomes were correlated (r = 0.81); 127 metabolites were significant for both, all in the same direction (Fig. 4). Excluding events in the first 2 years left 141 of the 157 (colon) and 166 of the 174 (stomach/duodenal) main hits significant, all with unchanged direction (Fig. 5). Among the main hits, the direction of effect was the same in women, men and participants aged 60 or older in 98–100% of cases for both outcomes; in participants under 60 years it was 98% for colon polyp and 91% for stomach/duodenal polyp, the lowest agreement of any stratum (65-year cut-off: 99–100%). Fifty-one (colon) and 75 (stomach/duodenal) main hits were also significant in all four sex and age strata. In split-sample validation, of the metabolites significant in the discovery half (134 and 149), all had the same direction in the validation half, 95% were nominally significant (P < 0.05) in it, and 41% and 77% also passed Bonferroni correction there; the validation half contained 30% of the sample. In bootstrap resamples, 156 of 157 and 168 of 174 main hits had P < 0.05 in at least 95% of resamples. These checks show stability within one cohort, not independent replication. In baseline logistic models, 75 (colon) and 46 (stomach/duodenal) metabolites were associated with prevalent disease.

### Dual-site polyps
Among eligible participants, 14,497 developed colon polyps only, 5,975 stomach/duodenal polyps only and 1,110 both (Fig. 6). Comparing participants with both polyps against those with one, one of 251 metabolites remained significant after Bonferroni correction: the phospholipid share of large HDL (OR 1.14 per SD, 95% CI 1.06–1.22). This comparison rests on 1,110 dual-site cases and is hypothesis-generating.

### Prediction
Adding metabolites to the covariates changed discrimination little (Fig. 7, Table S10–S11). For colon polyp, AUC rose from 0.626–0.627 to 0.630–0.631, the C-index from 0.630 to 0.636, and the 5- and 10-year time-dependent AUC from 0.676 to 0.682 and from 0.648 to 0.656. For stomach/duodenal polyp, AUC rose from 0.627–0.632 to 0.636–0.647, the C-index from 0.634–0.635 to 0.647–0.651, and the 5- and 10-year AUC from 0.656 to 0.681 and from 0.650 to 0.667. Overall discrimination remained low (AUC ≤ 0.68), and the Bonferroni-selected set performed like the full set of 251 metabolites.

## Limitations to state
- One cohort with no external validation.
- Hospital-recorded polyps only; polyps found at screening without admission are missed, and the outcome depends on diagnostic activity. Colonoscopy and screening history were not adjusted for.
- Model 1 does not adjust for diet, physical activity, lipid-lowering medication or family history.
- Metabolite measures are highly collinear; effect-size comparisons across metabolites are descriptive.
- The Cox penalty grid in the prediction step bottomed out at its lowest value, so gains from metabolites may be slightly underestimated.
