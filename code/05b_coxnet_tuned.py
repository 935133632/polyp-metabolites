"""Re-fit the penalised Cox model with the penalty chosen on an inner validation split of the training data
(the fixed penalty used in 05_prediction.py degraded C-index when metabolites were added). Same 70/30 split as 05."""
import sys; sys.path.insert(0, r"D:\research\Xingang Li\polyp_metabolites\code")
import warnings; warnings.filterwarnings("ignore")
from scipy.stats import ttest_ind
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from common import *
M = load_metabolome().set_index("eid"); MET = list(M.columns)
GRID = [0.0005, 0.002, 0.01, 0.05]
rows, sroc = [], []
for k in sys.argv[1:] or ["K63", "K31"]:
    c = build_cohort(k); Xc = design_covars(c).reset_index(drop=True); Xm = M.reindex(c.eid.values).reset_index(drop=True); X = pd.concat([Xc, Xm], axis=1)
    T, E = c.follow_up_days.values / 365.25, c.event.values
    tr, te = train_test_split(np.arange(len(c)), test_size=.3, stratify=E, random_state=1)
    t = ttest_ind(Xm.iloc[tr].values[E[tr] == 1], Xm.iloc[tr].values[E[tr] == 0], nan_policy="omit", equal_var=False)
    sel = [m for m, q in zip(MET, bonferroni(t.pvalue)) if q < .05]
    itr, iva = train_test_split(tr, test_size=.2, stratify=E[tr], random_state=2)
    for s, cols in (("covariates", list(Xc.columns)), ("cov+Bonf-selected", list(Xc.columns) + sel), ("cov+all251", list(X.columns))):
        Xf = X[cols].fillna(X[cols].iloc[tr].median()); Z = pd.DataFrame(StandardScaler().fit(Xf.iloc[tr]).transform(Xf), columns=cols)
        def fit(idx, pen):
            return CoxPHFitter(penalizer=pen, l1_ratio=.1).fit(pd.concat([Z.iloc[idx].reset_index(drop=True), pd.DataFrame({"T": T[idx], "E": E[idx]})], axis=1), "T", "E")
        score = {p: concordance_index(T[iva], -fit(itr, p).predict_partial_hazard(Z.iloc[iva]).values, E[iva]) for p in GRID}
        best = max(score, key=score.get); r = fit(tr, best).predict_partial_hazard(Z.iloc[te]).values
        row = {"outcome": k, "model": "CoxElasticNet_tuned", "feature_set": s, "n_features": len(cols), "penalizer": best, "c_index": concordance_index(T[te], -r, E[te])}
        for h in (5, 10):
            case = (E[te] == 1) & (T[te] <= h); ctrl = T[te] > h
            row[f"auc_{h}y"] = roc_auc_score(np.r_[np.ones(case.sum()), np.zeros(ctrl.sum())], np.r_[r[case], r[ctrl]])
        for h in (5, 10):
            case = (E[te] == 1) & (T[te] <= h); ctrl = T[te] > h
            sroc.append(roc_points(np.r_[np.ones(case.sum()), np.zeros(ctrl.sum())], np.r_[r[case], r[ctrl]]).assign(outcome=k, model="CoxElasticNet_tuned", feature_set=s, horizon_years=h))
        pd.concat(sroc).to_csv(RES / "step11_13_prediction" / f"survival_roc_coxnet_{k}.csv", index=False)
        rows.append(row); print(k, s, best, round(row["c_index"], 3), flush=True)
        pd.DataFrame(rows).to_csv(RES / "step11_13_prediction" / f"survival_coxnet_tuned_{k}.csv", index=False)
