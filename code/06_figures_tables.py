import sys; sys.path.insert(0, r"D:\research\Xingang Li\polyp_metabolites\code")
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.patches as mp
from common import *
from metab_annot import metab_class
plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"], "font.size": 8, "axes.spines.top": False,
                     "axes.spines.right": False, "pdf.fonttype": 42, "savefig.dpi": 300})
UP, DN, GR = "#D55E00", "#0072B2", "#BBBBBB"   # Wong palette
LAB = {"K63": "Colon polyp (K63.5)", "K31": "Stomach/duodenal polyp (K31.7)"}
FIG.mkdir(exist_ok=True); TAB.mkdir(exist_ok=True); (FIG / "supplement").mkdir(exist_ok=True)

def load(k, kind="step01_primary/model1_{k}_cox_all_metabolites.csv"):
    r = pd.read_csv(RES / kind.format(k=k)); r = r[r.status == "ok"].copy(); r["cls"] = r.metabolite.map(metab_class); return r

def save(fig, name): fig.savefig(FIG / f"{name}.png", bbox_inches="tight"); fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight"); plt.close(fig)

P = {k: load(k) for k in LAB}
# ---- Fig 1 volcano (labels repelled; 5 leading metabolites per direction)
from adjustText import adjust_text
fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.8))
for a, k in zip(ax, LAB):
    r = P[k]; y = -np.log10(r.p.clip(lower=1e-300)); sig = r.p_bonf < .05
    col = np.where(~sig, GR, np.where(r.coef > 0, UP, DN)); a.scatter(r.coef, y, c=col, s=8, lw=0)
    lab = pd.concat([r[r.coef > 0].nsmallest(5, "p"), r[r.coef < 0].nsmallest(5, "p")])
    texts = [a.text(t.coef, -np.log10(max(t.p, 1e-300)), t.metabolite, fontsize=5.5) for _, t in lab.iterrows()]
    a.axhline(-np.log10(.05 / len(r)), c="k", lw=.4, ls=":"); a.axvline(0, c="k", lw=.5)
    a.set_xlabel("ln(HR) per 1 SD"); a.set_ylabel("-log10(P)"); a.set_title(f"{LAB[k]}" + chr(10) + f"{sig.sum()}/{len(r)} Bonferroni<0.05", fontsize=8)
    a.margins(x=.12, y=.08); fig.canvas.draw()
    adjust_text(texts, ax=a, expand=(1.4, 1.6), arrowprops=dict(arrowstyle="-", color="k", lw=.3), time_limit=5)
fig.tight_layout(); save(fig, "Fig1_volcano")
# ---- Fig 2 class summary
rows = []
for k in LAB:
    for cls, g in P[k].groupby("cls"):
        s = g[g.p_bonf < .05]; rows.append({"outcome": k, "class": cls, "n": len(g), "n_sig_pos": int((s.coef > 0).sum()), "n_sig_neg": int((s.coef < 0).sum())})
cs = pd.DataFrame(rows); cs.to_csv(TAB / "class_summary.csv", index=False)
order = cs.groupby("class").n.first().sort_values().index.tolist()
fig, ax = plt.subplots(1, 2, figsize=(7.1, 3.4), sharey=True)
for a, k in zip(ax, LAB):
    d = cs[cs.outcome == k].set_index("class").reindex(order)
    a.barh(order, d.n_sig_pos, color=UP, label="Higher risk"); a.barh(order, -d.n_sig_neg, color=DN, label="Lower risk")
    a.axvline(0, c="k", lw=.5); a.set_xlabel("No. metabolites Bonferroni<0.05"); a.set_title(LAB[k], fontsize=8)
ax[0].legend(frameon=False, fontsize=7, loc="lower left"); fig.tight_layout(); save(fig, "Fig2_class_summary")
# ---- Fig 3 circular bar plot of all 251 metabolites (one figure per outcome; see circular_barplot.py)
from circular_barplot import circular_barplot
for k, tag in (("K63", "A_colon_polyp"), ("K31", "B_stomach_duodenal_polyp")):
    circular_barplot(P[k], LAB[k], str(FIG / f"Fig3{tag}"))
# supplementary S1: sensitivity analysis excluding events in the first 2 years, all 251 metabolites, one circular plot per outcome
for k, tag in (("K63", "A_colon_polyp"), ("K31", "B_stomach_duodenal_polyp")):
    sens = pd.read_csv(RES / f"step02_sensitivity_exclude_first2y/{k}_sensitivity.csv"); sens = sens[sens.status == "ok"]
    circular_barplot(sens, LAB[k] + ", excl. events <2 y", str(FIG / "supplement" / f"FigS1{tag}_sensitivity_circular"))
# ---- Fig 4 cross-outcome scatter + prevalent vs incident
j = P["K31"].merge(P["K63"], on="metabolite", suffixes=("_k31", "_k63"))
fig, ax = plt.subplots(1, 2, figsize=(7.1, 3.3))
b = (j.p_bonf_k31 < .05) & (j.p_bonf_k63 < .05)
ax[0].scatter(j.coef_k63, j.coef_k31, c=np.where(b, UP, GR), s=8, lw=0); lim = [min(j.coef_k63.min(), j.coef_k31.min()), max(j.coef_k63.max(), j.coef_k31.max())]
ax[0].plot(lim, lim, c="k", lw=.5, ls="--"); ax[0].set_xlabel("ln(HR) colon polyp"); ax[0].set_ylabel("ln(HR) stomach/duodenal polyp")
ax[0].set_title(f"Cox effects, r={j.coef_k63.corr(j.coef_k31):.2f}; {b.sum()} shared Bonferroni<0.05", fontsize=8)
pv = pd.read_csv(RES / "step03b_prevalent/K63_prevalent_logistic.csv").merge(P["K63"], on="metabolite")
ax[1].scatter(np.log(pv["or"]), pv.coef, c=GR, s=8, lw=0); ax[1].axhline(0, c="k", lw=.4); ax[1].axvline(0, c="k", lw=.4)
ax[1].set_xlabel("ln(OR) prevalent, colon polyp"); ax[1].set_ylabel("ln(HR) incident, colon polyp")
ax[1].set_title(f"Prevalent vs incident, r={np.log(pv['or']).corr(pv.coef):.2f}", fontsize=8)
fig.tight_layout(); save(fig, "Fig4_cross_outcome"); j.to_csv(TAB / "cross_outcome_cox.csv", index=False)
# ---- Fig 5: (A,B) circular heat rings of all 251 metabolites by sex and age (60) stratum; (C) sensitivity analysis excluding events <2 y
from circular_barplot import circular_heatmap
for k, tag in (("K63", "A_colon_polyp"), ("K31", "B_stomach_duodenal_polyp")):
    g = pd.read_csv(RES / f"step03_subgroup/{k}_subgroups.csv"); g = g[g.status == "ok"]
    rings = {"Overall": P[k]}
    for nm, lab in {"Female": "Female", "Male": "Male", "Age<60": "<60 y", "Age>=60": "≥60 y"}.items(): rings[lab] = g[g.subgroup == nm]
    circular_heatmap(rings, LAB[k], str(FIG / f"Fig5{tag}_strata_circular_heatmap"))
fig, ax = plt.subplots(1, 2, figsize=(7.1, 3.4))
for k, a in zip(LAB, ax):
    s = load(k, "step02_sensitivity_exclude_first2y/{k}_sensitivity.csv").merge(P[k], on="metabolite", suffixes=("_s", ""))
    a.scatter(s.coef, s.coef_s, c=GR, s=8, lw=0); a.plot([s.coef.min(), s.coef.max()], [s.coef.min(), s.coef.max()], "k--", lw=.5)
    a.set_xlabel("ln(HR) main"); a.set_ylabel("ln(HR) excl. events <2y"); a.set_title(f"{LAB[k]}, r={s.coef.corr(s.coef_s):.3f}", fontsize=8)
fig.tight_layout(); save(fig, "Fig5C_sensitivity_scatter")
# ---- Fig 6 dual site: (A) group sizes; (B) circular bar plot of all 251 metabolites, both polyps vs one (OR per 1 SD)
gc = pd.read_csv(RES / "step07_dual_site/group_counts.csv"); dl = pd.read_csv(RES / "step07_dual_site/dual_vs_single_logistic.csv")
fig, ax = plt.subplots(figsize=(3.6, 3.2))
g2 = gc.set_index("group").loc[["K63_only", "K31_only", "both"]]; ax.bar(["Colon only", "Stomach/duod. only", "Both"], g2.n, color=[UP, DN, "#009E73"])
for i, v in enumerate(g2.n): ax.text(i, v, f"{v:,}", ha="center", va="bottom", fontsize=7)
ax.set_ylabel("Participants with incident polyp"); ax.tick_params(axis="x", rotation=20); fig.tight_layout(); save(fig, "Fig6A_dual_site_counts")
circular_barplot(dl, "Both polyps (n=1,110) vs one", str(FIG / "Fig6B_dual_site_circular"), lo=.70, hi=1.30, est="or", est_lo="or_lower95", est_hi="or_upper95",
                 unit="OR per 1 SD", ticks=(.7, .85, 1.0, 1.15, 1.3))
# ---- tables
T1 = []
for k in LAB:
    r = P[k]; T1.append(r.assign(outcome=k)[["outcome", "metabolite", "cls", "n_complete", "n_events", "hr", "hr_lower95", "hr_upper95", "p", "p_bonf"]])
pd.concat(T1).to_csv(TAB / "TableS1_primary_cox_all_metabolites.csv", index=False)
for k in LAB: P[k].nsmallest(20, "p")[["metabolite", "cls", "n_events", "hr", "hr_lower95", "hr_upper95", "p", "p_bonf"]].to_csv(TAB / f"Table2_top20_{k}.csv", index=False)
cases = {k: build_cohort(k) for k in LAB}
rows = []
for k, c in cases.items():
    for lab, g in (("All", c), ("Incident polyp", c[c.event == 1]), ("No polyp", c[c.event == 0])):
        rows.append({"outcome": k, "group": lab, "N": len(g), "age_mean": g.AgeRecr.mean(), "male_pct": 100 * g.Sex.mean(), "BMI_mean": g.BMI.mean(), "TDI_mean": g.TDI.mean(),
                     "SBP_mean": g.SBPmean.mean(), "current_smoker_pct": 100 * (g.SmokeStatus == "Current").mean(), "median_FU_years": g.follow_up_years.median()})
pd.DataFrame(rows).round(2).to_csv(TAB / "Table1_cohort_characteristics.csv", index=False)
print("done")
