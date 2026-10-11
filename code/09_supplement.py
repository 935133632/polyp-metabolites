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

with pd.ExcelWriter(ST / "Supplementary_Tables.xlsx", engine="xlsxwriter") as xw:
    idx = pd.DataFrame({"sheet": list(sheets), "description": [v[0] for v in sheets.values()]}); idx.to_excel(xw, "Index", index=False)
    for name, (title, df) in sheets.items(): df.to_excel(xw, name[:31], index=False)
print("tables:", len(sheets))

# ------------------------------------------------------------------ figures
# S2 correlation structure and effective number of tests
Mt = load_metabolome().drop(columns="eid").sample(60000, random_state=1).dropna()
C = Mt.corr(); ev = np.linalg.eigvalsh(C.values)[::-1]; ev = ev[ev > 0]; cum = np.cumsum(ev) / ev.sum()
n95 = int(np.searchsorted(cum, .95) + 1)
order = sorted(C.columns, key=lambda n: (metab_class(n), n))
fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.4), gridspec_kw={"width_ratios": [1.2, 1]})
im = ax[0].imshow(C.loc[order, order].values, cmap="RdBu_r", vmin=-1, vmax=1); ax[0].set_xticks([]); ax[0].set_yticks([]); ax[0].set_title("Correlation of 251 metabolites\n(ordered by class)", fontsize=8)
fig.colorbar(im, ax=ax[0], shrink=.75, label="Pearson r")
ax[1].plot(np.arange(1, 61), cum[:60], c="k"); ax[1].axhline(.95, c=GR, ls="--"); ax[1].axvline(n95, c=GR, ls="--")
ax[1].set_xlabel("Number of principal components"); ax[1].set_ylabel("Cumulative variance explained"); ax[1].set_title(f"{n95} components explain 95%; Bonferroni over 251 is conservative", fontsize=7.5)
fig.tight_layout(); save(fig, "FigS2_correlation_effective_tests")

# S3 original 3-character definition vs polyp-specific codes
fig, ax = plt.subplots(1, 2, figsize=(7.1, 3.3)); ok = True
for a, k in zip(ax, LAB):
    p = rd(f"step01_primary/model1_{k}_cox_all_metabolites.csv"); q = rd(f"step01_primary/model1_{k}any_cox_all_metabolites.csv")
    if p is None or q is None: ok = False; break
    m = p.merge(q, on="metabolite", suffixes=("_s", "_a")); a.scatter(m.coef_s, m.coef_a, c=GR, s=8, lw=0); lim = [min(m.coef_s.min(), m.coef_a.min()), max(m.coef_s.max(), m.coef_a.max())]
    a.plot(lim, lim, "k--", lw=.5); a.set_xlabel(f"ln(HR), {k}.{'5' if k == 'K63' else '7'} (polyp-specific)"); a.set_ylabel(f"ln(HR), {k} (3-character code)")
    a.set_title(f"{LAB[k].split(' (')[0]}: r={m.coef_s.corr(m.coef_a):.2f}\n{(m.p_bonf_s < .05).sum()} vs {(m.p_bonf_a < .05).sum()} Bonferroni-significant", fontsize=8)
if ok: fig.tight_layout(); save(fig, "FigS3_original_vs_polyp_definition")

# S4 internal validation: split-sample concordance and bootstrap stability
fig, ax = plt.subplots(2, 2, figsize=(7.1, 5.6)); ok = True
for j, k in enumerate(LAB):
    s = rd(f"step06_internal_validation/{k}_split_sample.csv"); b = rd(f"step06_internal_validation/{k}_bootstrap_all_metabolites.csv"); p = rd(f"step01_primary/model1_{k}_cox_all_metabolites.csv")
    if s is None or b is None: ok = False; break
    a = ax[0, j]; sig = s.p_bonf_disc < .05; a.scatter(s.coef_disc, s.coef_valid, c=np.where(sig, UP, GR), s=8, lw=0); lim = [s.coef_disc.min(), s.coef_disc.max()]
    a.plot(lim, lim, "k--", lw=.5); a.set_xlabel("ln(HR) discovery (70%)"); a.set_ylabel("ln(HR) validation (30%)"); a.set_title(f"{LAB[k]}\nr={s.coef_disc.corr(s.coef_valid):.2f}", fontsize=8)
    a = ax[1, j]; sigm = set(p[p.p_bonf < .05].metabolite); bb = b[b.metabolite.isin(sigm)]
    a.hist(bb.freq_p05, bins=np.linspace(0, 1, 21), color=DN); a.set_xlabel("Fraction of bootstrap resamples with P<0.05"); a.set_ylabel("Main-analysis hits"); a.set_title(f"Bootstrap, {len(bb)} Bonferroni-significant metabolites", fontsize=8)
if ok: fig.tight_layout(); save(fig, "FigS4_internal_validation")

# S5 sex and age strata as circular heat rings (all 251 metabolites): overall, female, male, <60, >=60, <65, >=65
from circular_barplot import circular_heatmap
for k, tag in (("K63", "A_colon_polyp"), ("K31", "B_stomach_duodenal_polyp")):
    g = rd(f"step03_subgroup/{k}_subgroups.csv"); g65 = rd(f"step03_subgroup/{k}_subgroups_age65.csv"); p = rd(f"step01_primary/model1_{k}_cox_all_metabolites.csv")
    if g is None or g65 is None or p is None: continue
    allg = pd.concat([g, g65]); allg = allg[allg.status == "ok"]
    rings = {"Overall": p[p.status == "ok"]}
    for nm, lab in {"Female": "Female", "Male": "Male", "Age<60": "<60 y", "Age>=60": "≥60 y", "Age<65": "<65 y", "Age>=65": "≥65 y"}.items():
        rings[lab] = allg[allg.subgroup == nm]
    circular_heatmap(rings, LAB[k], str(SF / f"FigS5{tag}_strata_circular_heatmap"))

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
