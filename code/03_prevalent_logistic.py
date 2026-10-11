"""Baseline prevalent companion analysis: adjusted logistic regression per metabolite."""
import sys; sys.path.insert(0, r"D:\research\Xingang Li\polyp_metabolites\code")
import statsmodels.api as sm, warnings; warnings.filterwarnings("ignore")
from common import *
d = RES / "step03b_prevalent"; d.mkdir(exist_ok=True)
M = load_metabolome().set_index("eid")
for k in sys.argv[1:] or ["K63", "K31"]:
    c = prevalent_cohort(k); X = design_covars(c); Mc = M.reindex(c.eid.values); rows = []
    for m in M.columns:
        df = X.copy(); df["x"] = Mc[m].values; df["y"] = c.prevalent.values; df = df.dropna()
        if df.y.sum() < MIN_EVENTS or df.x.std() == 0: continue
        df["x"] = (df.x - df.x.mean()) / df.x.std()
        try:
            f = sm.Logit(df.y, sm.add_constant(df.drop(columns="y"))).fit(disp=0)
            ci = f.conf_int().loc["x"]
            rows.append({"metabolite": m, "n": len(df), "n_cases": int(df.y.sum()), "or": np.exp(f.params["x"]), "or_lower95": np.exp(ci[0]),
                         "or_upper95": np.exp(ci[1]), "p": f.pvalues["x"]})
        except Exception: pass
    r = pd.DataFrame(rows); r["p_bonf"] = bonferroni(r.p.values); r.sort_values("p").to_csv(d / f"{k}_prevalent_logistic.csv", index=False)
    print(k, len(r), "Bonferroni<.05:", (r.p_bonf < .05).sum(), flush=True)
