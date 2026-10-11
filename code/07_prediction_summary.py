import sys; sys.path.insert(0, r"D:\research\Xingang Li\polyp_metabolites\code")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from common import *
plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"], "font.size": 8, "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42, "savefig.dpi": 300})
d = RES / "step11_13_prediction"
b = pd.concat([pd.read_csv(d / f"binary_performance_{k}.csv") for k in ("K63", "K31")]); s = pd.concat([pd.read_csv(d / f"survival_performance_{k}.csv") for k in ("K63", "K31")]); t = pd.concat([pd.read_csv(d / f"survival_coxnet_tuned_{k}.csv") for k in ("K63", "K31")])
s = pd.concat([s[s.model != "CoxElasticNet"], t])          # drop the untuned (penalty too strong) Cox fit
b.round(4).to_csv(TAB / "Table_binary_prediction.csv", index=False); s.round(4).to_csv(TAB / "Table_survival_prediction.csv", index=False)
cols = {"covariates": "#999999", "cov+Bonf-selected": "#0072B2", "cov+all251": "#D55E00"}; LAB = {"K63": "Colon polyp (K63.5)", "K31": "Stomach/duodenal polyp (K31.7)"}
auc_ = lambda d: float(np.trapz(d.tpr, d.fpr))
# ---- Fig 7: ROC curves. Binary = out-of-fold predictions of the best cov+all251 model; survival = time-dependent ROC (cases: event by h; controls: follow-up > h), tuned penalised Cox
fig, ax = plt.subplots(2, 3, figsize=(7.4, 4.9))
for i, k in enumerate(LAB):
    br = pd.read_csv(d / f"binary_roc_{k}.csv"); best = b[(b.outcome == k) & (b.feature_set == "cov+all251")].sort_values("auc_mean").model.iloc[-1]
    sr = pd.read_csv(d / f"survival_roc_coxnet_{k}.csv")
    panels = [(br[br.model == best], f"{LAB[k]}" + chr(10) + f"Binary ({best}), 5-fold CV", None), (sr, "Survival, 5-year ROC", 5), (sr, "Survival, 10-year ROC", 10)]
    for j, (df, ttl, h) in enumerate(panels):
        a = ax[i, j]
        for fs, c in cols.items():
            dd = df[df.feature_set == fs] if h is None else df[(df.feature_set == fs) & (df.horizon_years == h)]
            a.plot(dd.fpr, dd.tpr, c=c, lw=1.1, label=f"{fs} ({auc_(dd):.3f})")
        a.plot([0, 1], [0, 1], c="k", lw=.5, ls="--"); a.set_xlabel("1 - specificity"); a.set_ylabel("Sensitivity"); a.set_title(ttl, fontsize=7); a.legend(frameon=False, fontsize=5.3, loc="lower right"); a.set_aspect("equal")
fig.tight_layout(); fig.savefig(FIG / "Fig7_prediction_ROC.png", bbox_inches="tight"); fig.savefig(FIG / "Fig7_prediction_ROC.pdf", bbox_inches="tight")
inc = []
for k in LAB:
    for nm, df, met in (("binary", b, "auc_mean"), ("survival", s, "c_index")):
        p = df[df.outcome == k].pivot_table(index="model", columns="feature_set", values=met)
        inc.append(pd.DataFrame({"outcome": k, "type": nm, "model": p.index, "covariates": p["covariates"].values, "cov+all251": p["cov+all251"].values, "gain": (p["cov+all251"] - p["covariates"]).values}))
g = pd.concat(inc); g.round(4).to_csv(TAB / "Table_prediction_gain.csv", index=False); print(g.round(3).to_string(index=False))
