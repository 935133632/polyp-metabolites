"""Prediction with three feature sets, all selection done inside the training data (no leakage):
   (A) covariates only; (B) covariates + metabolites passing Bonferroni (P<0.05/251) in a univariate screen of the training fold; (C) covariates + all 251 metabolites.
   Binary: 5-fold stratified CV. Survival: 70/30 split."""
import sys; sys.path.insert(0, r"D:\research\Xingang Li\polyp_metabolites\code")
import warnings; warnings.filterwarnings("ignore")
from scipy.stats import ttest_ind
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score
import lightgbm as lgb, xgboost as xgb
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from common import *
d = RES / "step11_13_prediction"; d.mkdir(exist_ok=True)
M = load_metabolome().set_index("eid"); MET = list(M.columns)


def screen(Xm, y):
    """Bonferroni univariate screen (Welch t, cases vs non-cases) on training data only."""
    t = ttest_ind(Xm[y == 1], Xm[y == 0], nan_policy="omit", equal_var=False)
    return [m for m, q in zip(MET, bonferroni(t.pvalue)) if q < .05]


def mk_models(spw):
    return {
        "Logistic": make_pipeline(SimpleImputer(strategy="median"), StandardScaler(), LogisticRegression(max_iter=500, class_weight="balanced")),
        "ElasticNet": make_pipeline(SimpleImputer(strategy="median"), StandardScaler(), LogisticRegression(penalty="elasticnet", solver="saga", l1_ratio=.5, C=.1, max_iter=300, class_weight="balanced")),
        "LightGBM": lgb.LGBMClassifier(n_estimators=300, learning_rate=.03, num_leaves=15, scale_pos_weight=spw, subsample=.8, subsample_freq=1, colsample_bytree=.8, verbose=-1, n_jobs=16),
        "XGBoost": xgb.XGBClassifier(n_estimators=300, learning_rate=.03, max_depth=4, scale_pos_weight=spw, subsample=.8, colsample_bytree=.8, n_jobs=16, tree_method="hist"),
    }


for k in sys.argv[1:] or ["K63", "K31"]:
    bin_rows, surv_rows, roc_rows = [], [], []
    c = build_cohort(k); y = c.event.values
    Xc = design_covars(c).reset_index(drop=True); Xm = M.reindex(c.eid.values).reset_index(drop=True)
    X = pd.concat([Xc, Xm], axis=1); spw = (y == 0).sum() / y.sum()
    # ---------- binary, 5-fold CV
    skf = StratifiedKFold(5, shuffle=True, random_state=1); nsel = []
    for name in mk_models(spw):
        res = {s: {"auc": [], "oof": np.zeros(len(y))} for s in "ABC"}
        for tr, te in skf.split(X, y):
            sel = screen(Xm.iloc[tr].values, y[tr]); nsel.append(len(sel))
            for s, cols in (("A", list(Xc.columns)), ("B", list(Xc.columns) + sel), ("C", list(X.columns))):
                m_ = mk_models(spw)[name]; m_.fit(X.iloc[tr][cols], y[tr]); p = m_.predict_proba(X.iloc[te][cols])[:, 1]
                res[s]["oof"][te] = p; res[s]["auc"].append(roc_auc_score(y[te], p))
        for s, lab in (("A", "covariates"), ("B", "cov+Bonf-selected"), ("C", "cov+all251")):
            bin_rows.append({"outcome": k, "model": name, "feature_set": lab, "auc_mean": np.mean(res[s]["auc"]), "auc_sd": np.std(res[s]["auc"]),
                             "pr_auc": average_precision_score(y, res[s]["oof"]), "n_selected_mean": np.mean(nsel) if s == "B" else np.nan, "prevalence": y.mean()})
        for s, lab in (("A", "covariates"), ("B", "cov+Bonf-selected"), ("C", "cov+all251")):
            roc_rows.append(roc_points(y, res[s]["oof"]).assign(outcome=k, model=name, feature_set=lab, auc_oof=roc_auc_score(y, res[s]["oof"])))
        pd.concat(roc_rows).to_csv(d / f"binary_roc_{k}.csv", index=False)
        print(k, name, {s: round(np.mean(res[s]["auc"]), 3) for s in "ABC"}, flush=True)
        pd.DataFrame(bin_rows).to_csv(d / f"binary_performance_{k}.csv", index=False)
    # ---------- survival, 70/30
    sroc = []
    T, E = c.follow_up_days.values / 365.25, c.event.values
    tr, te = train_test_split(np.arange(len(c)), test_size=.3, stratify=E, random_state=1)
    sel = screen(Xm.iloc[tr].values, E[tr])
    for s, cols in (("covariates", list(Xc.columns)), ("cov+Bonf-selected", list(Xc.columns) + sel), ("cov+all251", list(X.columns))):
        Xf = X[cols].fillna(X[cols].iloc[tr].median()); sc = StandardScaler().fit(Xf.iloc[tr]); Z = pd.DataFrame(sc.transform(Xf), columns=cols)
        risk = {}
        cph = CoxPHFitter(penalizer=.05, l1_ratio=.5).fit(pd.concat([Z.iloc[tr].reset_index(drop=True), pd.DataFrame({"T": T[tr], "E": E[tr]})], axis=1), "T", "E")
        risk["CoxElasticNet"] = cph.predict_partial_hazard(Z.iloc[te]).values
        b = xgb.XGBRegressor(objective="survival:cox", n_jobs=16, tree_method="hist", n_estimators=300, learning_rate=.03, max_depth=3, subsample=.8, colsample_bytree=.8).fit(Z.iloc[tr], np.where(E[tr] == 1, T[tr], -T[tr]))
        risk["XGBoostCox"] = b.predict(Z.iloc[te])
        for name, r in risk.items():
            row = {"outcome": k, "model": name, "feature_set": s, "n_features": len(cols), "c_index": concordance_index(T[te], -r, E[te])}
            for h in (5, 10):  # cumulative/dynamic AUC without IPCW: cases = event by h, controls = follow-up > h
                case = (E[te] == 1) & (T[te] <= h); ctrl = T[te] > h
                row[f"auc_{h}y"] = roc_auc_score(np.r_[np.ones(case.sum()), np.zeros(ctrl.sum())], np.r_[r[case], r[ctrl]])
            surv_rows.append(row)
            if name == "XGBoostCox":
                for h in (5, 10):   # time-dependent ROC: cases = event by h years, controls = follow-up beyond h
                    case = (E[te] == 1) & (T[te] <= h); ctrl = T[te] > h
                    sroc.append(roc_points(np.r_[np.ones(case.sum()), np.zeros(ctrl.sum())], np.r_[r[case], r[ctrl]]).assign(outcome=k, model=name, feature_set=s, horizon_years=h))
                pd.concat(sroc).to_csv(d / f"survival_roc_{k}.csv", index=False)
        print(k, s, [(r["model"], round(r["c_index"], 3)) for r in surv_rows[-2:]], flush=True)
        pd.DataFrame(surv_rows).to_csv(d / f"survival_performance_{k}.csv", index=False)
