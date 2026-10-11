import sys; sys.path.insert(0, r"D:\research\Xingang Li\polyp_metabolites\code")
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.patches as mp
from common import *
from metab_annot import metab_class
plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"], "font.size": 8, "axes.spines.top": False,
                     "axes.spines.right": False, "pdf.fonttype": 42, "savefig.dpi": 300})
UP, DN, GR = "#D55E00", "#0072B2", "#BBBBBB"   # Wong palette
LAB = {"K63": "Colon polyp (K63.5)", "K31": "Stomach/duodenal polyp (K31.7)"}
FIG.mkdir(exist_ok=True); TAB.mkdir(exist_ok=True)

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
# ---- Fig 3 circular forest, all 251 metabolites (HR and 95% CI on a radial log scale)
CLS_ORDER = ["Lipid content of subclass", "Lipid composition (% within subclass)", "Lipoprotein / lipid totals", "Lipoprotein particle concentration",
             "Lipoprotein particle size", "Apolipoproteins", "Fatty acids", "Amino acids", "Glycolysis / ketone bodies", "Fluid balance", "Inflammation (GlycA)"]
PAL = ["#0072B2", "#56B4E9", "#009E73", "#F0E442", "#CC79A7", "#999999", "#D55E00", "#E69F00", "#117733", "#882255", "#000000"]
CC = dict(zip(CLS_ORDER, PAL))
def circular_forest(ax, r, title, lo=-.24, hi=.24):
    r = r.copy(); r["ci"] = pd.Categorical(r.cls, CLS_ORDER, ordered=True); r = r.sort_values(["ci", "coef"]).reset_index(drop=True)
    n = len(r); gap = 2 * np.pi * .05; th = np.linspace(np.pi / 2 - gap / 2, np.pi / 2 - 2 * np.pi + gap / 2, n)   # clockwise from top
    rad = lambda v: .12 + (np.clip(v, lo, hi) - lo) / (hi - lo) * .74
    g = np.linspace(0, 2 * np.pi, 400)
    for v in (-.2, .2): ax.plot(g, np.full_like(g, rad(v)), c="#E3E3E3", lw=.4, zorder=0)
    ax.plot(g, np.full_like(g, rad(0)), c="k", lw=.6, ls="--", zorder=1)
    for t_, (_, m) in zip(th, r.iterrows()):
        c_ = CC[m.cls] if m.p_bonf < .05 else "#BBBBBB"
        ax.plot([t_, t_], [rad(np.log(m.hr_lower95)), rad(np.log(m.hr_upper95))], c=c_, lw=.9, solid_capstyle="butt", zorder=2)
        ax.plot(t_, rad(m.coef), "o", ms=2.2, c=c_, mew=0, zorder=3)
    for cls_, g_ in r.groupby("ci", observed=True):   # class arc outside the data ring
        i0, i1 = g_.index.min(), g_.index.max(); arc = np.linspace(th[i0], th[i1], 40)
        ax.plot(arc, np.full_like(arc, .93), c=CC[cls_], lw=4, solid_capstyle="butt")
    top = pd.concat([r[r.coef > 0].nsmallest(2, "p"), r[r.coef < 0].nsmallest(2, "p")]).sort_index()
    for j, (i, m) in enumerate(top.iterrows()):
        deg = np.degrees(th[i]) % 360; left = 90 < deg < 270
        ax.text(th[i] + .04 * (-1) ** j, .99, m.metabolite, fontsize=5, rotation=deg + (180 if left else 0), ha="right" if left else "left", va="center", rotation_mode="anchor")
    ax.set_ylim(0, 1.32); ax.set_yticks([]); ax.set_xticks([]); ax.spines["polar"].set_visible(False); ax.set_title(title, fontsize=8, pad=16)
fig = plt.figure(figsize=(7.4, 4.9))
for i, k in enumerate(LAB):
    circular_forest(fig.add_subplot(1, 2, i + 1, projection="polar"), P[k], f"{LAB[k]} (n=251)")
h = [mp.Patch(color=CC[c_], label=c_) for c_ in CLS_ORDER] + [mp.Patch(color="#CCCCCC", label="Not significant (Bonferroni)")]
fig.legend(handles=h, loc="lower center", ncol=3, frameon=False, fontsize=5.8, bbox_to_anchor=(.5, .05))
fig.text(.5, .21, "Dashed ring: HR = 1; light rings: HR 0.82 and 1.22 per SD. Dot = HR, line = 95% CI; grey = not significant after Bonferroni correction.", fontsize=5.3, ha="center", color="#444444")
fig.subplots_adjust(left=.02, right=.98, top=.93, bottom=.24, wspace=.02); save(fig, "Fig3_circular_forest")
# supplementary: linear forest of the 20 leading metabolites per outcome
fig, ax = plt.subplots(1, 2, figsize=(7.1, 4.6))
for a, k in zip(ax, LAB):
    r = P[k].nsmallest(20, "p").iloc[::-1]; yy = np.arange(len(r))
    a.errorbar(r.hr, yy, xerr=[r.hr - r.hr_lower95, r.hr_upper95 - r.hr], fmt="o", ms=3, c="k", lw=.8)
    a.axvline(1, c=GR, lw=.7); a.set_yticks(yy, r.metabolite); a.set_xlabel("HR per 1 SD (95% CI)"); a.set_title(LAB[k], fontsize=8)
fig.tight_layout(); save(fig, "FigS1_forest_top20")
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
# ---- Fig 5 sensitivity + subgroup heatmap
fig, ax = plt.subplots(1, 3, figsize=(7.1, 4.6), gridspec_kw={"width_ratios": [1, 1.1, 1.1]})
for k, a in zip(["K63"], [ax[0]]):
    s = load(k, "step02_sensitivity_exclude_first2y/{k}_sensitivity.csv").merge(P[k], on="metabolite", suffixes=("_s", ""))
    a.scatter(s.coef, s.coef_s, c=GR, s=8, lw=0); a.plot([s.coef.min(), s.coef.max()], [s.coef.min(), s.coef.max()], "k--", lw=.5)
    a.set_xlabel("ln(HR) main"); a.set_ylabel("ln(HR) excl. events <2y"); a.set_title(f"{LAB[k]}\nr={s.coef.corr(s.coef_s):.3f}", fontsize=8)
for k, a in zip(LAB, ax[1:]):
    g = pd.read_csv(RES / f"step03_subgroup/{k}_subgroups.csv"); g = g[g.status == "ok"]
    H = g.pivot(index="metabolite", columns="subgroup", values="hr")[["Female", "Male", "Age<60", "Age>=60"]]
    H = H.loc[P[k].metabolite[P[k].metabolite.isin(H.index)].head(20).tolist()]   # display only: 20 leading metabolites; all 251 were analysed
    im = a.imshow(np.log(H.values), cmap="RdBu_r", vmin=-.3, vmax=.3, aspect="auto"); a.set_xticks(range(4), H.columns, rotation=45, ha="right"); a.set_yticks(range(len(H)), H.index, fontsize=6)
    a.set_title(LAB[k].split(" (")[0] + " subgroups", fontsize=8)
fig.colorbar(im, ax=ax[2], shrink=.7, label="ln(HR)"); fig.tight_layout(); save(fig, "Fig5_sensitivity_subgroup")
# ---- Fig 6 dual site
gc = pd.read_csv(RES / "step07_dual_site/group_counts.csv"); dl = pd.read_csv(RES / "step07_dual_site/dual_vs_single_logistic.csv").nsmallest(12, "p").iloc[::-1]
fig, ax = plt.subplots(1, 2, figsize=(7.1, 3.2))
g2 = gc.set_index("group").loc[["K63_only", "K31_only", "both"]]; ax[0].bar(["Colon only", "Stomach/duod. only", "Both"], g2.n, color=[UP, DN, "#009E73"])
for i, v in enumerate(g2.n): ax[0].text(i, v, f"{v:,}", ha="center", va="bottom", fontsize=7)
ax[0].set_ylabel("Participants with incident polyp"); ax[0].tick_params(axis="x", rotation=20)
yy = np.arange(len(dl)); ax[1].errorbar(dl["or"], yy, xerr=[dl["or"] - dl.or_lower95, dl.or_upper95 - dl["or"]], fmt="o", ms=3, c="k", lw=.8); ax[1].axvline(1, c=GR, lw=.7)
ax[1].set_yticks(yy, dl.metabolite); ax[1].set_xlabel("OR per 1 SD, both sites vs single site"); fig.tight_layout(); save(fig, "Fig6_dual_site")
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
