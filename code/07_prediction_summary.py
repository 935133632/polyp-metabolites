import sys; sys.path.insert(0, r"D:\research\Xingang Li\polyp_metabolites\code")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from common import *
plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"], "font.size": 8, "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42, "savefig.dpi": 300})
d = RES / "step11_13_prediction"
b = pd.concat([pd.read_csv(d / f"binary_performance_{k}.csv") for k in ("K63", "K31")]); s = pd.concat([pd.read_csv(d / f"survival_performance_{k}.csv") for k in ("K63", "K31")]); t = pd.concat([pd.read_csv(d / f"survival_coxnet_tuned_{k}.csv") for k in ("K63", "K31")])
s = pd.concat([s[s.model != "CoxElasticNet"], t])          # drop the untuned (penalty too strong) Cox fit
b.round(4).to_csv(TAB / "Table_binary_prediction.csv", index=False); s.round(4).to_csv(TAB / "Table_survival_prediction.csv", index=False)
cols = {"covariates": "#BBBBBB", "cov+Bonf-selected": "#0072B2", "cov+all251": "#D55E00"}; LAB = {"K63": "Colon polyp (K63.5)", "K31": "Stomach/duodenal polyp (K31.7)"}
fig, ax = plt.subplots(2, 2, figsize=(7.1, 5))
for i, k in enumerate(LAB):
    for j, (df, metric, ttl) in enumerate(((b, "auc_mean", "Binary, 5-fold CV AUC"), (s, "c_index", "Survival, test C-index"))):
        a = ax[i, j]; dd = df[df.outcome == k]; models = list(dd.model.unique()); x = np.arange(len(models))
        for m_, (fs, c) in enumerate(cols.items()):
            v = dd[dd.feature_set == fs].set_index("model").reindex(models)[metric]; a.bar(x + (m_ - 1) * .27, v, .27, color=c, label=fs)
        a.set_xticks(x, models, rotation=20, ha="right"); a.set_ylim(.5, .70); a.axhline(.5, c="k", lw=.5); a.set_ylabel(metric.replace("_mean", "")); a.set_title(f"{LAB[k]}\n{ttl}", fontsize=8)
ax[0, 0].legend(frameon=False, fontsize=6.5); fig.tight_layout(); fig.savefig(FIG / "Fig7_prediction.png", bbox_inches="tight"); fig.savefig(FIG / "Fig7_prediction.pdf", bbox_inches="tight")
inc = []
for k in LAB:
    for nm, df, met in (("binary", b, "auc_mean"), ("survival", s, "c_index")):
        p = df[df.outcome == k].pivot_table(index="model", columns="feature_set", values=met)
        inc.append(pd.DataFrame({"outcome": k, "type": nm, "model": p.index, "covariates": p["covariates"].values, "cov+all251": p["cov+all251"].values, "gain": (p["cov+all251"] - p["covariates"]).values}))
g = pd.concat(inc); g.round(4).to_csv(TAB / "Table_prediction_gain.csv", index=False); print(g.round(3).to_string(index=False))
