"""Metabolite risk score (penalised Cox fitted on 70% training data, evaluated in the 30% test set) and dose-response (restricted cubic splines).
usage: python 14_score_dose_response.py <K63|K31> <score|spline>"""
import sys; sys.path.insert(0, r"D:\research\Xingang Li\polyp_metabolites\code")
import warnings; warnings.filterwarnings("ignore")
from scipy.stats import chi2, norm
from sklearn.model_selection import train_test_split
from lifelines import CoxPHFitter, KaplanMeierFitter
from common import *

k, mode = sys.argv[1], sys.argv[2]
d = RES / "step14_score_dose_response"; d.mkdir(exist_ok=True)
c = build_cohort(k); M = load_metabolome().set_index("eid"); MET = list(M.columns)
Xc = design_covars(c).reset_index(drop=True); Xm = M.reindex(c.eid.values).reset_index(drop=True)
T, E = c.follow_up_days.values / 365.25, c.event.values

if mode in ("score", "cuminc"):
    tr, te = train_test_split(np.arange(len(c)), test_size=.3, stratify=E, random_state=1)       # same split as the prediction analysis
    med_c = Xc.iloc[tr].median(); med_m = Xm.iloc[tr].median()
    Xtr = pd.concat([Xc.iloc[tr].fillna(med_c), Xm.iloc[tr].fillna(med_m)], axis=1)
    mu, sd = Xtr.mean(), Xtr.std().replace(0, 1)
    Z = lambda X: ((X - mu) / sd)
    fit_df = pd.concat([Z(Xtr).reset_index(drop=True), pd.DataFrame({"T": T[tr], "E": E[tr]})], axis=1)
    if mode == "score":
        cph = CoxPHFitter(penalizer=.0005, l1_ratio=.1).fit(fit_df, "T", "E")
        beta = cph.params_[MET]
        pd.DataFrame({"metabolite": MET, "weight": beta.values}).sort_values("weight", key=abs, ascending=False).to_csv(d / f"{k}_score_weights.csv", index=False)
    else:                                                      # "cuminc": reuse the saved weights (no refit)
        beta = pd.read_csv(d / f"{k}_score_weights.csv").set_index("metabolite").weight.reindex(MET)
    # score = metabolite part of the linear predictor; evaluated in the test set only
    Xte_m = Xm.iloc[te].fillna(med_m); score = ((Xte_m - mu[MET]) / sd[MET]).values @ beta.values
    t = pd.concat([Xc.iloc[te].reset_index(drop=True), pd.DataFrame({"T": T[te], "E": E[te], "score": score})], axis=1)
    t = t.fillna(t.median())
    t["score_z"] = (t.score - t.score.mean()) / t.score.std()
    cov = [c_ for c_ in Xc.columns]
    r = CoxPHFitter().fit(t[cov + ["score_z", "T", "E"]], "T", "E").summary.loc["score_z"]
    rows = [{"term": "per 1 SD", "hr": r["exp(coef)"], "lower95": r["exp(coef) lower 95%"], "upper95": r["exp(coef) upper 95%"], "p": r["p"], "n": len(t), "events": int(t.E.sum())}]
    t["q"] = pd.qcut(t.score_z, 4, labels=False)
    for q in (1, 2, 3): t[f"Q{q + 1}"] = (t.q == q).astype(float)
    s = CoxPHFitter().fit(t[cov + ["Q2", "Q3", "Q4", "T", "E"]], "T", "E").summary
    for q in ("Q2", "Q3", "Q4"):
        rows.append({"term": f"{q} vs Q1", "hr": s.loc[q, "exp(coef)"], "lower95": s.loc[q, "exp(coef) lower 95%"], "upper95": s.loc[q, "exp(coef) upper 95%"], "p": s.loc[q, "p"]})
    t["qmed"] = t.groupby("q").score_z.transform("median")                      # trend across quartile medians
    tr_ = CoxPHFitter().fit(t[cov + ["qmed", "T", "E"]], "T", "E").summary.loc["qmed"]
    rows.append({"term": "P-trend (quartile medians)", "hr": tr_["exp(coef)"], "lower95": tr_["exp(coef) lower 95%"], "upper95": tr_["exp(coef) upper 95%"], "p": tr_["p"]})
    pd.DataFrame(rows).to_csv(d / f"{k}_score_quartiles.csv", index=False)
    grid = np.linspace(0, 14, 57); out = []
    for q in range(4):
        km = KaplanMeierFitter().fit(t["T"][t.q == q], t["E"][t.q == q]); out.append(pd.DataFrame({"quartile": f"Q{q + 1}", "years": grid, "cum_incidence": 1 - km.predict(grid).values}))
    pd.concat(out).to_csv(d / f"{k}_score_cumincidence.csv", index=False)
    print(k, pd.DataFrame(rows).round(3).to_string(index=False))

else:   # restricted cubic splines (4 knots at the 5th, 35th, 65th and 95th percentiles) for key metabolites
    p = pd.read_csv(RES / f"step01_primary/model1_{k}_cox_all_metabolites.csv").query("status == 'ok'").sort_values("p")
    pos, neg = p[p.coef > 0].metabolite.head(3).tolist(), p[p.coef < 0].metabolite.head(3).tolist()
    keys = list(dict.fromkeys(pos + neg + ["GlycoAcet", "HDLC", "MetaTG"]))

    def rcs(x, knots):
        t1, t2, t3, t4 = knots; pp = lambda v: np.maximum(v, 0) ** 3
        h = lambda tj: (pp(x - tj) - pp(x - t3) * (t4 - tj) / (t4 - t3) + pp(x - t4) * (t3 - tj) / (t4 - t3)) / (t4 - t1) ** 2
        return np.column_stack([x, h(t1), h(t2)])
    rows, nl = [], []
    for m in keys:
        df = pd.concat([Xc, pd.DataFrame({"x": Xm[m].values, "T": T, "E": E})], axis=1).dropna().reset_index(drop=True)
        xv = df.x.values; cv = df[list(Xc.columns)]; te_ = df[["T", "E"]]
        knots = np.percentile(xv, [5, 35, 65, 95]); B = rcs(xv, knots)
        lin = CoxPHFitter().fit(pd.concat([cv, df[["x"]], te_], axis=1), "T", "E")
        sp = CoxPHFitter().fit(pd.concat([cv, pd.DataFrame(B, columns=["s1", "s2", "s3"]), te_], axis=1), "T", "E")
        lr = 2 * (sp.log_likelihood_ - lin.log_likelihood_); pn = chi2.sf(lr, 2); nl.append({"metabolite": m, "p_nonlinear": pn, "n": len(df), "events": int(df.E.sum())})
        V = sp.variance_matrix_.loc[["s1", "s2", "s3"], ["s1", "s2", "s3"]].values; b = sp.params_[["s1", "s2", "s3"]].values
        g = np.linspace(-2.5, 2.5, 51); Bg = rcs(g, knots); B0 = rcs(np.array([0.0]), knots); D = Bg - B0
        lp = D @ b; se = np.sqrt(np.einsum("ij,jk,ik->i", D, V, D))
        rows.append(pd.DataFrame({"metabolite": m, "z": g, "hr": np.exp(lp), "lower95": np.exp(lp - 1.96 * se), "upper95": np.exp(lp + 1.96 * se)}))
        print(k, m, "P-nonlinear", round(pn, 4), flush=True)
    pd.concat(rows).to_csv(d / f"{k}_dose_response.csv", index=False); pd.DataFrame(nl).to_csv(d / f"{k}_nonlinearity.csv", index=False)
