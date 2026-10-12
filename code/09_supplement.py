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


# ------------------------------------------------------------------ tables (5 supplementary tables)
for f in ST.glob("TableS*.csv"): f.unlink()          # regenerate from scratch

def pick(df, cols, prefix):
    d = df[["metabolite"] + [c for c in cols if c in df.columns]].copy()
    return d.rename(columns={c: f"{prefix}_{c}" for c in cols})

wide = []
for k in LAB:
    p = rd(f"step01_primary/model1_{k}_cox_all_metabolites.csv")
    W = pick(p, ["n_events", "hr", "hr_lower95", "hr_upper95", "p", "p_bonf"], "primary")
    s_ = rd(f"step02_sensitivity_exclude_first2y/{k}_sensitivity.csv"); W = W.merge(pick(s_, ["n_events", "hr", "hr_lower95", "hr_upper95", "p_bonf"], "sens_excl2y"), on="metabolite", how="left")
    v = rd(f"step03b_prevalent/{k}_prevalent_logistic.csv"); W = W.merge(pick(v, ["n_cases", "or", "or_lower95", "or_upper95", "p_bonf"], "prevalent"), on="metabolite", how="left")
    g = rd(f"step03_subgroup/{k}_subgroups.csv"); g = g[g.status == "ok"]
    for nm, lab in {"Female": "female", "Male": "male", "Age<60": "age_lt60", "Age>=60": "age_ge60"}.items():
        W = W.merge(pick(g[g.subgroup == nm], ["n_events", "hr", "hr_lower95", "hr_upper95", "p_bonf"], f"sub_{lab}"), on="metabolite", how="left")
    sp_ = rd(f"step06_internal_validation/{k}_split_sample.csv"); W = W.merge(pick(sp_, ["coef_disc", "coef_valid", "p_bonf_disc", "p_bonf_valid", "same_direction"], "split"), on="metabolite", how="left")
    b = rd(f"step06_internal_validation/{k}_bootstrap_all_metabolites.csv"); W = W.merge(pick(b, ["freq_p05", "coef_p2_5", "coef_p97_5"], "boot"), on="metabolite", how="left")
    for m in ("M1", "M2a", "M2b", "M2c", "M3", "M4"):
        a = rd(f"step12_model_adjustment/{k}_{m}.csv"); W = W.merge(pick(a, ["hr", "hr_lower95", "hr_upper95", "p_bonf"], f"adj_{m}"), on="metabolite", how="left")
    o = rd(f"step01_primary/model1_{k}any_cox_all_metabolites.csv"); W = W.merge(pick(o, ["n_events", "hr", "hr_lower95", "hr_upper95", "p_bonf"], "orig3char"), on="metabolite", how="left")
    sw = rd(f"step14_score_dose_response/{k}_score_weights.csv"); W = W.merge(sw.rename(columns={"weight": "score_weight"}), on="metabolite", how="left")
    dl = rd("step07_dual_site/dual_vs_single_logistic.csv"); W = W.merge(pick(dl, ["or", "or_lower95", "or_upper95", "p_bonf"], "dual_site"), on="metabolite", how="left")
    W.insert(0, "outcome", LAB[k]); wide.append(W)
add("TableS1_metabolite_results", pd.concat(wide).sort_values(["outcome", "primary_p"]),
    "One row per outcome and metabolite: primary Cox (Model 1; HR per 1 SD), sensitivity excluding events <2 y, prevalent-disease odds ratio, sex and age-60 subgroups, split-sample (ln HR, 70% discovery / 30% validation), bootstrap (30 resamples; ln HR interval), extended adjustment in the same complete-case sample (M1; M2a +waist, body fat %; M2b +HbA1c, log CRP; M2c +waist, body fat %, HbA1c, log CRP; M3 +prior-diagnosis proxies; M4 all), original 3-character code definition, metabolite-score weight, and dual-site vs single-site OR (identical in both outcome blocks). p_bonf = Bonferroni-adjusted P")

M = load_metabolome(False); names = [c for c in M.columns if c != "eid"]
ann = pd.DataFrame({"metabolite": names, "lipoprotein_family": [lipo_family(n) for n in names], "n_nonmissing": [int(M[n].notna().sum()) for n in names]})
ann.insert(1, "class", [metab_class(n) for n in names]); sheets["TableS2_metabolite_annotation"] = ("Metabolite abbreviations (as supplied) with rule-based class and lipoprotein family", ann); ann.to_csv(ST / "TableS2_metabolite_annotation.csv", index=False)

bp = pd.concat([pd.read_csv(RES / f"step11_13_prediction/binary_performance_{k}.csv").assign(type="binary, 5-fold CV") for k in LAB])
sp = pd.concat([pd.read_csv(RES / f"step11_13_prediction/survival_performance_{k}.csv") for k in LAB]); sp = sp[sp.model != "CoxElasticNet"].assign(type="survival, 70/30 split")
st = pd.concat([pd.read_csv(RES / f"step11_13_prediction/survival_coxnet_tuned_{k}.csv") for k in LAB]).assign(type="survival, 70/30 split")
pred = pd.concat([bp, sp, st])[["type", "outcome", "model", "feature_set", "n_features", "auc_mean", "auc_sd", "pr_auc", "n_selected_mean", "c_index", "auc_5y", "auc_10y", "penalizer"]]
add("TableS3_prediction_performance", pred, "Prediction performance: binary (5-fold cross-validation AUC) and survival (test-set C-index, 5- and 10-year AUC) for three feature sets")

def cat(pattern, names, label):
    parts = [rd(pattern.format(n=n)) for n in names]; parts = [p.assign(**{label: n}) for p, n in zip(parts, names) if p is not None]
    return pd.concat(parts) if parts else None
codes = ["K57", "K21", "K29", "K44", "K59", "K62", "K64", "K92", "K80", "K40", "K22", "K52", "K76", "K58", "K20", "K30", "K56", "K66", "K25", "K42", "K83", "K43", "K26", "CRC"]
add("TableS4_specificity_correlations", rd("step15_summary/specificity_profile_correlation.csv"), "Correlation of each outcome's 251 ln(HR) values with the colon and stomach/duodenal polyp profiles (events and number of Bonferroni-significant metabolites)")
add("TableS5_specificity_effects", cat("step13_specificity/{n}.csv", codes, "outcome_code"), "Model 1 Cox of all 251 metabolites for 23 other digestive outcomes (3-character ICD-10 codes) and colorectal cancer (CRC = C18-C20)")

add("TableS6_lipid_audit", rd("step16_lipid_audit/TableS6_lipid_audit.csv"), "Audit of the lipoprotein results: NMR vs clinical chemistry (alignment and independent hazard ratios), stepwise adjustment for prior-diagnosis proxies, mutual adjustment, crude incidence by HDL tertile")

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

# S5 prevalent vs incident, both outcomes
fig, ax = plt.subplots(1, 2, figsize=(7.1, 3.3))
for a, k in zip(ax, LAB):
    v = rd(f"step03b_prevalent/{k}_prevalent_logistic.csv"); p = rd(f"step01_primary/model1_{k}_cox_all_metabolites.csv")
    m = v.merge(p, on="metabolite"); x = np.log(m["or"]); a.scatter(x, m.coef, c=GR, s=8, lw=0); a.axhline(0, c="k", lw=.4); a.axvline(0, c="k", lw=.4)
    a.set_xlabel("ln(OR) prevalent disease"); a.set_ylabel("ln(HR) incident disease"); a.set_title(f"{LAB[k]}: r={x.corr(m.coef):.2f}", fontsize=8)
fig.tight_layout(); save(fig, "FigS5_prevalent_vs_incident")

# S6 prediction performance bars
if len(bp) and len(sp):
    b = bp; s = pd.concat([sp, st])
    cols = {"covariates": "#BBBBBB", "cov+Bonf-selected": "#0072B2", "cov+all251": "#D55E00"}
    fig, ax = plt.subplots(2, 2, figsize=(7.1, 5))
    for i, k in enumerate(LAB):
        for j, (df, metric, ttl) in enumerate(((b, "auc_mean", "Binary, 5-fold CV AUC"), (s, "c_index", "Survival, test C-index"))):
            a = ax[i, j]; dd = df[df.outcome == k]; models = list(dd.model.unique()); x = np.arange(len(models))
            for m_, (fs, c) in enumerate(cols.items()):
                a.bar(x + (m_ - 1) * .27, dd[dd.feature_set == fs].set_index("model").reindex(models)[metric], .27, color=c, label=fs)
            a.set_xticks(x, models, rotation=20, ha="right"); a.set_ylim(.5, .70); a.axhline(.5, c="k", lw=.5); a.set_ylabel(metric.replace("_mean", "")); a.set_title(f"{LAB[k]}\n{ttl}", fontsize=8)
    ax[0, 0].legend(frameon=False, fontsize=6.5); fig.tight_layout(); save(fig, "FigS6_prediction_performance_bars")

# S7 dual-site timing and group sizes
gc, tc = rd("step07_dual_site/group_counts.csv"), rd("step07_dual_site/timing_counts.csv")
if gc is not None and tc is not None:
    fig, ax = plt.subplots(1, 2, figsize=(7.1, 3)); g = gc[gc.group != "none"]
    ax[0].bar(g.group, g.n, color=[UP, DN, "#009E73"]); ax[0].set_ylabel("Participants"); ax[0].set_title("Incident polyp groups (no polyp: %s)" % f"{int(gc[gc.group == 'none'].n.iloc[0]):,}", fontsize=8)
    ax[1].bar(tc.timing, tc.n, color="#009E73"); ax[1].set_ylabel("Dual-site cases"); ax[1].set_title("Order of diagnosis in dual-site cases", fontsize=8)
    fig.tight_layout(); save(fig, "FigS7_dual_site_timing")
print("figures done")
