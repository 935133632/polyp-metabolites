# Before submission (working notes, not part of the manuscript)

## 1. References
All 43 references have DOIs resolved and bibliographic details taken from Crossref, are cited at least once, and are stored in the Zotero library (My Library > collection `polyp_metabolites`, tag `polyp_metabolites`). `references_verified.csv` lists them in citation order.
- Five items already existed elsewhere in the Zotero library (Bray 2024, Sudlow 2015, Yuan 2025, Harrell 2015, Chen 2016) and now appear twice; merge them with Zotero's Duplicate Items view so the collection keeps one copy each.
- Claims were checked against the abstracts of: Tian 2015 (triglyceride RR 1.06 and LDL-C RR 1.11 for adenoma; no HDL-C association), Yao 2014 (triglycerides RR 1.18; HDL-C 0.84, 0.69-1.02), Wu 2021 (OR 1.39), El Herrag 2017 (RR 1.43), Lu 2022 (no consistent association of the metabolic syndrome with metachronous neoplasia), Li 2025 (UK Biobank plasma PUFA, HR per SD 0.90-0.96), Atashi 2025 (dietary linoleic acid RR 1.15; tissue levels null), Chandler 2016 (GlycA and colorectal cancer), Yang 2025 (HDL subclasses positively associated with GI cancers), Yuan 2025, Ouyang 2021, Willer 2013 (157 lipid loci also associated with type 2 diabetes, waist-hip ratio and BMI), Voight 2012 (genetic HDL-C score not associated with myocardial infarction). Read the full texts of these before submission and confirm the sentences that cite them.
- Statements that rest on the title or general knowledge of the paper and need a full-text check: Bray 2024 (colorectal cancer burden), Winawer 1993 and Zauber 2012 (polypectomy and cancer incidence and mortality), Nordestgaard 2016 (triglyceride-rich lipoproteins as cardiovascular risk factors), Buergel 2022 and Julkunen 2023 (multi-disease prediction and association atlas), Kettunen 2018 and Connelly 2017 (GlycA), Schaeffer 2006 (FADS1/FADS2 and phospholipid fatty-acid composition; title-level), Lu 2023 (fatty acids and colorectal cancer inconsistent), Fry 2017 (healthy-volunteer selection).

## 2. Statements to confirm
- Introduction, paragraph 3: the statement that the metabolic profile preceding polyps has received less attention. Confirm with a systematic literature search and soften or remove it if earlier analyses exist.
- Introduction, paragraph 1: "adenomas can progress to cancer" is supported by the gastric polyp reviews (refs 5, 6) only for the stomach; add a colorectal reference if wanted.

## 3. Data dictionary (resolved by arithmetic checks)
Metabolite definitions were verified from exact identities in the data (see ANALYSIS_LOG, Step 14): `TC_HDL` is non-HDL cholesterol, `TC_Chy_VLDL` is the triglyceride measure of the chylomicron/XXL-VLDL series, `TPR` is the triglyceride-to-phosphoglyceride ratio, `Phospho` is phosphoglycerides, and `GluLac` and `SpeAla` are corrected glucose and alanine. Compare these with the official Nightingale/UK Biobank field list once if the supplementary annotation (Table S2) is published.

## 4. Lipoprotein results
The audit changed the interpretation of LDL cholesterol (confounding by treatment and health-care contact) and of HDL cholesterol (stomach/duodenal polyp only). Prior hospital diagnoses are an incomplete proxy for medication use; if medication (lipid-lowering drugs, aspirin, proton-pump inhibitors) and endoscopy/screening records can be exported later, repeat M3/M4 with them.
