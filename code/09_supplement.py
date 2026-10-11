"""Supplementary tables (Excel workbook + CSVs) and supplementary figures. Items whose inputs are missing are skipped with a message."""
import sys; sys.path.insert(0, r"D:\research\Xingang Li\polyp_metabolites\code")
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import *
from metab_annot import metab_class, lipo_family

plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"], "font.size": 8, "axes.spines.top": False,
                     "axes.spines.right": False, "pdf.fonttype": 42, "savefig.dpi": 300})
UP, DN, GR = "#D55E00", "#0072B2", "#BBBBBB"
LAB = {"K63": "Colon polyp (K63.5)", "K31": "Stomach/duodenal polyp (K31.7)"}
ST, SF = TAB / "supplement", FIG / "supplement"; ST.mkdir(parents=True, exist_ok=True); SF.mkdir(parents=True, exist_ok=True)
sheets = {}


def rd(path):
    f = RES / path
    if not f.exists(): print("skip (missing):", path); return None
    return pd.read_csv(f)


def both(pattern, **kw):
    out = []
    for k in LAB:
        d = rd(pattern.format(k=k))
        if d is not None: out.append(d.assign(outcome=k))
    return pd.concat(out) if out else None


def add(name, df, title):
    if df is None: return
    df = df.copy()
    if "metabolite" in df.columns and "class" not in df.columns: df.insert(df.columns.get_loc("metabolite") + 1, "class", df.metabolite.map(metab_class))
    sheets[name] = (title, df); df.round(6).to_csv(ST / f"{name}.csv", index=False)


def save(fig, name): fig.savefig(SF / f"{name}.png", bbox_inches="tight"); fig.savefig(SF / f"{name}.pdf", bbox_inches="tight"); plt.close(fig)


# ------------------------------------------------------------------ tables
add("TableS1_primary_cox", both("step01_primary/model1_{k}_cox_all_metabolites.csv"), "Primary Cox, Model 1, all 251 metabolites (HR per 1 SD; p_bonf = Bonferroni-adjusted P)")
add("TableS2_sensitivity_excl_2y", both("step02_sensitivity_exclude_first2y/{k}_sensitivity.csv"), "Sensitivity analysis excluding events in the first 2 years")
add("TableS3_subgroups_sex_age60", both("step03_subgroup/{k}_subgroups.csv"), "Subgroup Cox: female, male, age <60, age >=60 (all 251 metabolites)")
add("TableS4_subgroups_age65", both("step03_subgroup/{k}_subgroups_age65.csv"), "Supplementary age cut-off: age <65 and >=65")
add("TableS5_prevalent_logistic", both("step03b_prevalent/{k}_prevalent_logistic.csv"), "Baseline prevalent disease, adjusted logistic regression (OR per 1 SD)")
add("TableS6_split_sample", both("step06_internal_validation/{k}_split_sample.csv"), "Split-sample validation (70% discovery / 30% validation), all 251 metabolites")
add("TableS7_bootstrap", both("step06_internal_validation/{k}_bootstrap_all_metabolites.csv"), "Bootstrap stability (30 resamples), all 251 metabolites")
add("TableS8_dual_site", rd("step07_dual_site/dual_vs_single_logistic.csv"), "Dual-site (both polyps, n=1,110) vs single-site cases, all 251 metabolites")
orig = [rd(f"step01_primary/model1_{k}any_cox_all_metabolites.csv") for k in LAB]
add("TableS9_original_3char_definition", pd.concat([d.assign(outcome=f"{k}any") for k, d in zip(LAB, orig) if d is not None]) if any(d is not None for d in orig) else None,
    "Primary Cox with the original 3-character codes (K63, K31) used in the Olink project")
bp = [pd.read_csv(RES / f"step11_13_prediction/binary_performance_{k}.csv") for k in LAB if (RES / f"step11_13_prediction/binary_performance_{k}.csv").exists()]
sp = [pd.read_csv(RES / f"step11_13_prediction/survival_performance_{k}.csv") for k in LAB if (RES / f"step11_13_prediction/survival_performance_{k}.csv").exists()]
st = [pd.read_csv(RES / f"step11_13_prediction/survival_coxnet_tuned_{k}.csv") for k in LAB if (RES / f"step11_13_prediction/survival_coxnet_tuned_{k}.csv").exists()]
if bp: add("TableS10_binary_prediction", pd.concat(bp), "Binary prediction, 5-fold cross-validation, three feature sets")
if sp: add("TableS11_survival_prediction", pd.concat([pd.concat(sp)[lambda d: d.model != "CoxElasticNet"]] + st), "Survival prediction, 70/30 split, three feature sets")
M = load_metabolome(False); names = [c for c in M.columns if c != "eid"]
ann = pd.DataFrame({"metabolite": names, "class": [metab_class(n) for n in names], "lipoprotein_family": [lipo_family(n) for n in names],
                    "n_nonmissing": [int(M[n].notna().sum()) for n in names]})
ann = pd.concat([ann.drop(columns="class"), ann["class"]], axis=1); sheets["TableS12_metabolite_annotation"] = ("Metabolite abbreviations (as supplied) with rule-based class", ann); ann.to_csv(ST / "TableS12_metabolite_annotation.csv", index=False)
t1 = pd.read_csv(TAB / "Table1_cohort_characteristics.csv"); sheets["TableS13_cohort"] = ("Cohort characteristics", t1); t1.to_csv(ST / "TableS13_cohort.csv", index=False)

# extended analyses (steps 12-14)
def cat(pattern, names, label):
    parts = [rd(pattern.format(n=n)) for n in names]; parts = [p.assign(**{label: n}) for p, n in zip(parts, names) if p is not None]
    return pd.concat(parts) if parts else None
adj = [rd(f"step12_model_adjustment/{k}_{m}.csv") for k in LAB for m in ("M1", "M2a", "M2b", "M2c")]
add("TableS14_model_adjustment", pd.concat([a.assign(outcome=k) for a, k in zip(adj, [k for k in LAB for _ in range(4)]) if a is not None]) if all(a is not None for a in adj) else None,
    "Primary Cox in the same complete-case sample: M1; M2a (+waist, body fat %); M2b (+HbA1c, log CRP); M2c (+all four)")
add("TableS15_model_adjustment_summary", rd("step15_summary/model_adjustment_summary.csv"), "Attenuation of main hits after extended adjustment")
codes = ["K57", "K21", "K29", "K44", "K59", "K62", "K64", "K92", "K80", "K40", "K22", "K52", "K76", "K58", "K20", "K30", "K56", "K66", "K25", "K42", "K83", "K43", "K26", "CRC"]
add("TableS16_specificity_effects", cat("step13_specificity/{n}.csv", codes, "outcome_code"), "Model 1 Cox of all 251 metabolites for other digestive outcomes (K = 3-character ICD-10; CRC = C18-C20)")
add("TableS17_specificity_summary", rd("step15_summary/specificity_profile_correlation.csv"), "Correlation of each outcome's 251 ln(HR) values with the colon and stomach/duodenal polyp profiles")
sq = cat("step14_score_dose_response/{n}_score_quartiles.csv", list(LAB), "outcome"); add("TableS18_score_quartiles", sq, "Metabolite score in the 30% test set: HR per SD and by quartile")
add("TableS19_score_weights", cat("step14_score_dose_response/{n}_score_weights.csv", list(LAB), "outcome"), "Score weights (penalised Cox, training set)")
add("TableS20_dose_response", cat("step14_score_dose_response/{n}_dose_response.csv", list(LAB), "outcome"), "Restricted cubic spline curves (HR vs median)")
add("TableS21_nonlinearity", cat("step14_score_dose_response/{n}_nonlinearity.csv", list(LAB), "outcome"), "P for non-linearity (likelihood ratio, 4-knot spline vs linear)")

with pd.ExcelWriter(ST / "Supplementary_Tables.xlsx", engine="xlsxwriter") as xw:
    idx = pd.DataFrame({"sheet": list(sheets), "description": [v[0] for v in sheets.values()]}); idx.to_excel(xw, "Index", index=False)
    for name, (title, df) in sheets.items(): df.to_excel(xw, name[:31], index=False)
print("tables:", len(sheets))

# ------------------------------------------------------------------ figures
# S2 correlation structure: (A) circular heat rings of the correlation of each metabolite with principal components 1-8; (B) cumulative variance
from circular_barplot import circular_heatmap
Mt = load_metabolome().drop(columns="eid").sample(60000, random_state=1).dropna()
Z = (Mt - Mt.mean()) / Mt.std(); C = Z.corr(); evals, evecs = np.linalg.eigh(C.values); o = np.argsort(evals)[::-1]; evals, evecs = evals[o], evecs[:, o]
cum = np.cumsum(evals) / evals.sum(); n95 = int(np.searchsorted(cum, .95) + 1)
rings = {f"PC{j + 1} ({100 * evals[j] / evals.sum():.0f}%)": pd.DataFrame({"metabolite": C.columns, "coef": evecs[:, j] * np.sqrt(evals[j]), "p_bonf": 0.0}) for j in range(8)}
circular_heatmap(rings, "Correlation of each metabolite with principal components", str(SF / "FigS2A_correlation_structure_circular"), vmax=1.0, note="rings: PC1 (inner) to PC8; % = variance explained", cbar_label="Correlation with PC")
fig, a = plt.subplots(figsize=(3.8, 3.2))
a.plot(np.arange(1, 61), cum[:60], c="k"); a.axhline(.95, c=GR, ls="--"); a.axvline(n95, c=GR, ls="--")
a.set_xlabel("Number of principal components"); a.set_ylabel("Cumulative variance explained"); a.set_title(f"{n95} components explain 95%", fontsize=8); fig.tight_layout(); save(fig, "FigS2B_effective_tests_scree")

# S3 original 3-character definition vs polyp-specific codes
fig, ax = plt.subplots(1, 2, figsize=(7.1, 3.3)); ok = True
for a, k in zip(ax, LAB):
    p = rd(f"step01_primary/model1_{k}_cox_all_metabolites.csv"); q = rd(f"step01_primary/model1_{k}any_cox_all_metabolites.csv")
    if p is None or q is None: ok = False; break
    m = p.merge(q, on="metabolite", suffixes=("_s", "_a")); a.scatter(m.coef_s, m.coef_a, c=GR, s=8, lw=0); lim = [min(m.coef_s.min(), m.coef_a.min()), max(m.coef_s.max(), m.coef_a.max())]
    a.plot(lim, lim, "k--", lw=.5); a.set_xlabel(f"ln(HR), {k}.{'5' if k == 'K63' else '7'} (polyp-specific)"); a.set_ylabel(f"ln(HR), {k} (3-character code)")
    a.set_title(f"{LAB[k].split(' (')[0]}: r={m.coef_s.corr(m.coef_a):.2f}\n{(m.p_bonf_s < .05).sum()} vs {(m.p_bonf_a < .05).sum()} Bonferroni-significant", fontsize=8)
if ok: fig.tight_layout(); save(fig, "FigS3_original_vs_polyp_definition")

# S4 internal validation: split-sample concordance (top) and bootstrap caterpillar plot of all 251 metabolites (bottom)
fig, ax = plt.subplots(2, 2, figsize=(7.1, 6.4)); ok = True
for j, k in enumerate(LAB):
    s_ = rd(f"step06_internal_validation/{k}_split_sample.csv"); b = rd(f"step06_internal_validation/{k}_bootstrap_all_metabolites.csv"); p = rd(f"step01_primary/model1_{k}_cox_all_metabolites.csv")
    if s_ is None or b is None or "coef_p2_5" not in b: ok = False; break
    a = ax[0, j]; sig = s_.p_bonf_disc < .05; a.scatter(s_.coef_disc, s_.coef_valid, c=np.where(sig, UP, GR), s=8, lw=0); lim = [s_.coef_disc.min(), s_.coef_disc.max()]
    a.plot(lim, lim, "k--", lw=.5); a.set_xlabel("ln(HR) discovery (70%)"); a.set_ylabel("ln(HR) validation (30%)"); a.set_title(f"{LAB[k]}" + chr(10) + f"split-sample, r={s_.coef_disc.corr(s_.coef_valid):.2f}", fontsize=8)
    a = ax[1, j]; m = p[p.status == "ok"].merge(b, on="metabolite").sort_values("coef").reset_index(drop=True); x = np.arange(len(m))
    a.fill_between(x, m.coef_p2_5, m.coef_p97_5, color="#9DB7D5", alpha=.8, lw=0, label="Bootstrap 95% interval (30 resamples)")
    col = np.where(m.p_bonf < .05, np.where(m.coef > 0, UP, DN), GR); a.scatter(x, m.coef, c=col, s=5, lw=0, zorder=3, label="Main analysis")
    a.axhline(0, c="k", lw=.5); a.set_xlabel("Metabolites ranked by main-analysis ln(HR)"); a.set_ylabel("ln(HR) per 1 SD"); a.set_title(f"Bootstrap stability, all {len(m)} metabolites", fontsize=8)
    a.legend(frameon=False, fontsize=6, loc="upper left")
if ok: fig.tight_layout(); save(fig, "FigS4_internal_validation")

# S5 supplementary age cut-off 65 as circular heat rings (all 251 metabolites): overall, <65, >=65
for k, tag in (("K63", "A_colon_polyp"), ("K31", "B_stomach_duodenal_polyp")):
    g65 = rd(f"step03_subgroup/{k}_subgroups_age65.csv"); p = rd(f"step01_primary/model1_{k}_cox_all_metabolites.csv")
    if g65 is None or p is None: continue
    g65 = g65[g65.status == "ok"]
    rings = {"Overall": p[p.status == "ok"], "<65 y": g65[g65.subgroup == "Age<65"], "≥65 y": g65[g65.subgroup == "Age>=65"]}
    circular_heatmap(rings, LAB[k] + ", age cut-off 65", str(SF / f"FigS5{tag}_age65_circular_heatmap"))

# S6 prevalent vs incident, both outcomes
fig, ax = plt.subplots(1, 2, figsize=(7.1, 3.3))
for a, k in zip(ax, LAB):
    v = rd(f"step03b_prevalent/{k}_prevalent_logistic.csv"); p = rd(f"step01_primary/model1_{k}_cox_all_metabolites.csv")
    m = v.merge(p, on="metabolite"); x = np.log(m["or"]); a.scatter(x, m.coef, c=GR, s=8, lw=0); a.axhline(0, c="k", lw=.4); a.axvline(0, c="k", lw=.4)
    a.set_xlabel("ln(OR) prevalent disease"); a.set_ylabel("ln(HR) incident disease"); a.set_title(f"{LAB[k]}: r={x.corr(m.coef):.2f}", fontsize=8)
fig.tight_layout(); save(fig, "FigS6_prevalent_vs_incident")

# S7 prediction performance bars
if bp and sp:
    b = pd.concat(bp); s = pd.concat([pd.concat(sp)[lambda d: d.model != "CoxElasticNet"]] + st)
    cols = {"covariates": "#BBBBBB", "cov+Bonf-selected": "#0072B2", "cov+all251": "#D55E00"}
    fig, ax = plt.subplots(2, 2, figsize=(7.1, 5))
    for i, k in enumerate(LAB):
        for j, (df, metric, ttl) in enumerate(((b, "auc_mean", "Binary, 5-fold CV AUC"), (s, "c_index", "Survival, test C-index"))):
            a = ax[i, j]; dd = df[df.outcome == k]; models = list(dd.model.unique()); x = np.arange(len(models))
            for m_, (fs, c) in enumerate(cols.items()):
                a.bar(x + (m_ - 1) * .27, dd[dd.feature_set == fs].set_index("model").reindex(models)[metric], .27, color=c, label=fs)
            a.set_xticks(x, models, rotation=20, ha="right"); a.set_ylim(.5, .70); a.axhline(.5, c="k", lw=.5); a.set_ylabel(metric.replace("_mean", "")); a.set_title(f"{LAB[k]}\n{ttl}", fontsize=8)
    ax[0, 0].legend(frameon=False, fontsize=6.5); fig.tight_layout(); save(fig, "FigS7_prediction_performance_bars")

# S8 dual-site timing and group sizes
gc, tc = rd("step07_dual_site/group_counts.csv"), rd("step07_dual_site/timing_counts.csv")
if gc is not None and tc is not None:
    fig, ax = plt.subplots(1, 2, figsize=(7.1, 3)); g = gc[gc.group != "none"]
    ax[0].bar(g.group, g.n, color=[UP, DN, "#009E73"]); ax[0].set_ylabel("Participants"); ax[0].set_title("Incident polyp groups (no polyp: %s)" % f"{int(gc[gc.group == 'none'].n.iloc[0]):,}", fontsize=8)
    ax[1].bar(tc.timing, tc.n, color="#009E73"); ax[1].set_ylabel("Dual-site cases"); ax[1].set_title("Order of diagnosis in dual-site cases", fontsize=8)
    fig.tight_layout(); save(fig, "FigS8_dual_site_timing")
print("figures done")
