"""Main-text tables in publication format.
Table 1: baseline characteristics by incident polyp status (two outcomes).
Table 2: hazard ratios of key metabolites in Model 1, Model 3 (+ prior diagnoses) and Model 4 (fully adjusted), same complete-case sample."""
import sys; sys.path.insert(0, r"D:\research\Xingang Li\polyp_metabolites\code")
import warnings; warnings.filterwarnings("ignore")
from scipy import stats
from common import *
from docx_helpers import *

OUT = TAB; OUT.mkdir(exist_ok=True)
LAB = {"K63": "Colon polyp (K63.5)", "K31": "Stomach/duodenal polyp (K31.7)"}


def pfmt(p): return "<0.001" if p < 0.001 else f"{p:.3f}"
def ms(x, d=1): x = x.dropna(); return f"{x.mean():.{d}f} ({x.std():.{d}f})"
def mi(x, d=1):
    x = x.dropna(); q = x.quantile([.25, .5, .75]); sep = " to " if q[.25] < 0 else "–"
    return f"{q[.5]:.{d}f} ({q[.25]:.{d}f}{sep}{q[.75]:.{d}f})"
def cnt(x, v): x = x.dropna(); n = int((x == v).sum()); return f"{n:,} ({100 * n / len(x):.1f})"


def table1():
    extra = pd.read_pickle(PROC / "extra_covars.pkl"); extra["CRP_mgL"] = np.exp(extra.CRP)
    prior = pd.read_pickle(PROC / "prior_dx.pkl"); raw = load_metabolome(False)[["eid", "MetaTG", "HDLC", "LDLC", "GlycoAcet"]]
    VARS = [("Age at recruitment, years", "AgeRecr", "ms", 1), ("Female sex, n (%)", "Sex", "cnt0", None), ("Townsend deprivation index", "TDI", "mi", 2), ("Body mass index, kg/m²", "BMI", "ms", 1),
            ("Waist circumference, cm", "Waist", "ms", 1), ("Systolic blood pressure, mmHg", "SBPmean", "ms", 0), ("Diastolic blood pressure, mmHg", "DBPmean", "ms", 0),
            ("Smoking status, n (%)", None, "head", None), ("Never", "SmokeStatus", "cntNever", None), ("Previous", "SmokeStatus", "cntPrevious", None), ("Current", "SmokeStatus", "cntCurrent", None),
            ("Alcohol drinking, n (%)", None, "head", None), ("Never ", "DrinkStatus", "cntNever", None), ("Previous ", "DrinkStatus", "cntPrevious", None), ("Current ", "DrinkStatus", "cntCurrent", None),
            ("HbA1c, mmol/mol", "HbA1c", "mi", 1), ("C-reactive protein, mg/L", "CRP_mgL", "mi", 2),
            ("Prior hospital diagnoses, number of ICD-10 codes", "n_prior_dx", "mi", 0), ("Disorders of lipoprotein metabolism (E78), n (%)", "prior_E78", "cnt1", None),
            ("Ischaemic heart disease (I20–I25), n (%)", "prior_IHD", "cnt1", None), ("Hypertension (I10–I15), n (%)", "prior_HTN", "cnt1", None), ("Diabetes (E10–E14), n (%)", "prior_DM", "cnt1", None),
            ("Total triglycerides, mmol/L", "MetaTG", "mi", 2), ("HDL cholesterol, mmol/L", "HDLC", "ms", 2), ("LDL cholesterol, mmol/L", "LDLC", "ms", 2), ("Glycoprotein acetyls, mmol/L", "GlycoAcet", "mi", 2),
            ("Follow-up, years", "follow_up_years", "mi", 1)]
    cols, rows = {}, {}
    for k in LAB:
        c = build_cohort(k).merge(extra, on="eid", how="left").merge(prior, on="eid", how="left").merge(raw, on="eid", how="left")
        g = {"All": c, "Case": c[c.event == 1], "None": c[c.event == 0]}; cols[k] = {n: len(d) for n, d in g.items()}
        for lab, col, kind, dec in VARS:
            if kind == "head": rows.setdefault(lab, {})[k] = ["", "", "", ""]; continue
            out = []
            for n in ("All", "Case", "None"):
                x = g[n][col]
                out.append(ms(x, dec) if kind == "ms" else mi(x, dec) if kind == "mi" else cnt(x, 0.0 if kind == "cnt0" else 1) if kind in ("cnt0", "cnt1") else cnt(x, kind[3:]))
            a, b = g["Case"][col].dropna(), g["None"][col].dropna()
            if kind == "ms": p = stats.ttest_ind(a, b, equal_var=False).pvalue
            elif kind == "mi": p = stats.mannwhitneyu(a, b).pvalue
            else:
                lv = pd.concat([g["Case"][col], g["None"][col]]); isv = (lv == (0.0 if kind == "cnt0" else 1 if kind == "cnt1" else kind[3:])).astype(int)
                grp = np.r_[np.ones(len(g["Case"][col])), np.zeros(len(g["None"][col]))]; ok = lv.notna().values
                p = stats.chi2_contingency(pd.crosstab(isv.values[ok], grp[ok]))[1]
            rows.setdefault(lab, {})[k] = out + [pfmt(p)]
    return VARS, cols, rows


def table2():
    MET = [("Lipids and lipoproteins", None), ("Total triglycerides", "MetaTG"), ("HDL cholesterol", "HDLC"), ("Apolipoprotein A1", "ApoA1"), ("LDL cholesterol", "LDLC"), ("Apolipoprotein B", "ApoB"),
           ("Total cholesterol", "MetaTC"), ("Triglycerides in very small VLDL, % of lipids", "TG_VSVLDL_pt"), ("Inflammation", None), ("Glycoprotein acetyls", "GlycoAcet"), ("Albumin", "Albumin"),
           ("Fatty acids, % of total fatty acids", None), ("Linoleic acid", "Linopt"), ("Omega-6 fatty acids", "O6pt"), ("Monounsaturated fatty acids", "Monopt"), ("Saturated fatty acids", "Satupt")]
    R = {(k, m): pd.read_csv(RES / f"step12_model_adjustment/{k}_{m}.csv").query("status == 'ok'").set_index("metabolite") for k in LAB for m in ("M1", "M3", "M4")}
    body, bold, groups = [], [], []
    for name, key in MET:
        if key is None: groups.append(len(body)); body.append([name] + [""] * 6); continue
        row = [name]
        for j, k in enumerate(LAB):
            for mj, m in enumerate(("M1", "M3", "M4")):
                r = R[(k, m)].loc[key]; row.append(f"{r.hr:.2f} ({r.hr_lower95:.2f}–{r.hr_upper95:.2f})")
                if r.p_bonf < 0.05: bold.append((len(body), 1 + j * 3 + mj))
        body.append(row)
    ns = {k: (int(R[(k, 'M1')].n_complete.max()), int(R[(k, 'M1')].n_events.max())) for k in LAB}
    return body, bold, groups, ns


if __name__ == "__main__":
    VARS, cols, rows = table1()
    head1 = ["Characteristic"] + [LAB["K63"], "", "", "", LAB["K31"], "", "", ""]
    head2 = [""] + sum([[f"All participants\n(n = {cols[k]['All']:,})", f"Incident polyp\n(n = {cols[k]['Case']:,})", f"No polyp\n(n = {cols[k]['None']:,})", "P"] for k in LAB], [])
    body = [[lab] + rows[lab]["K63"] + rows[lab]["K31"] for lab, *_ in VARS]
    indent = [i for i, (lab, col, kind, d) in enumerate(VARS) if kind.startswith("cntNever") or kind.startswith("cntPrevious") or kind.startswith("cntCurrent")]
    heads = [i for i, v in enumerate(VARS) if v[2] == "head"]
    pd.DataFrame([head2] + body, columns=head1).to_csv(OUT / "Table1_baseline_characteristics.csv", index=False, encoding="utf-8-sig")
    d = new_document(size=9, line_spacing=1.0); landscape(d.sections[0])
    caption(d, "**Table 1.** Baseline characteristics of UK Biobank participants by incident colon polyp and incident stomach/duodenal polyp status", size=10)
    three_line_table(d, [head1, head2], body, col_widths_cm=[6.2] + [2.2, 2.2, 2.2, 1.3] * 2, font_size=7.5, merges=[(0, 1, 4), (0, 5, 8)], indent_rows=indent, group_rows=heads)
    caption(d, "Values are mean (standard deviation), median (interquartile range) or n (%); percentages exclude missing values. Participants with the outcome before baseline were excluded from the corresponding cohort. "
               "P values compare participants with and without incident polyp (Welch t test, Mann–Whitney test or chi-squared test). Prior hospital diagnoses are those recorded before baseline. Metabolite concentrations are NMR measurements.", size=8)
    d.save(OUT / "Table1_baseline_characteristics.docx")

    body2, bold2, groups2, ns = table2()
    h1 = ["Metabolite", LAB["K63"], "", "", LAB["K31"], "", ""]; h2 = [""] + ["Model 1", "Model 3", "Model 4"] * 2
    pd.DataFrame([h2] + body2, columns=h1).to_csv(OUT / "Table2_key_metabolites.csv", index=False, encoding="utf-8-sig")
    d = new_document(size=9, line_spacing=1.0); landscape(d.sections[0])
    caption(d, "**Table 2.** Associations of key metabolites with incident colon polyp and stomach/duodenal polyp across adjustment models (hazard ratio per 1 SD, 95% CI)", size=10)
    three_line_table(d, [h1, h2], body2, col_widths_cm=[7.0] + [3.0] * 6, font_size=8, merges=[(0, 1, 3), (0, 4, 6)], bold_cells=bold2, indent_rows=[i for i in range(len(body2)) if i not in groups2], group_rows=groups2)
    caption(d, f"Same complete-case sample for all models (colon polyp: {ns['K63'][0]:,} participants, {ns['K63'][1]:,} events; stomach/duodenal polyp: {ns['K31'][0]:,}, {ns['K31'][1]:,}). "
               "Model 1: age, sex, Townsend deprivation index, body mass index, systolic and diastolic blood pressure, smoking and alcohol drinking. Model 3: Model 1 plus prior hospital diagnoses (log number of ICD-10 codes; disorders of lipoprotein metabolism, ischaemic heart disease, hypertension, diabetes, reflux disease, gastritis or ulcer). "
               "Model 4: Model 3 plus waist circumference, body fat percentage, HbA1c and C-reactive protein. Metabolites were rank-normalised and standardised. Bold: P < 0.05 after Bonferroni correction across 251 metabolites. Results for all 251 metabolites are in Table S1 and Fig. 3.", size=8)
    d.save(OUT / "Table2_key_metabolites.docx")
    with pd.ExcelWriter(OUT / "Main_Tables.xlsx", engine="xlsxwriter") as xw:
        pd.DataFrame([head2] + body, columns=head1).to_excel(xw, "Table 1", index=False); pd.DataFrame([h2] + body2, columns=h1).to_excel(xw, "Table 2", index=False)
    print("Table 1 rows:", len(body), "| Table 2 rows:", len(body2), "| complete-case n:", ns)
    print(pd.DataFrame([h2] + body2, columns=h1).to_string(index=False))
