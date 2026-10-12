"""Adjustment for prior-diagnosis proxies (lipid-lowering/aspirin indication, health-care contact, reflux/gastritis).
usage: python 17_proxy_adjustment.py scan <K63|K31>   -> all 251 metabolites, Model 1 + proxies (model "M3")
       python 17_proxy_adjustment.py decomp          -> which proxies matter, for key lipid metabolites
       python 17_proxy_adjustment.py defs            -> arithmetic checks of metabolite definitions"""
import sys; sys.path.insert(0, r"D:\research\Xingang Li\polyp_metabolites\code")
import warnings; warnings.filterwarnings("ignore")
from lifelines import CoxPHFitter
from common import *

OUT = RES / "step17_proxy_adjustment"; OUT.mkdir(exist_ok=True)
prior = pd.read_pickle(PROC / "prior_dx.pkl"); prior["log_prior_dx"] = np.log1p(prior.n_prior_dx)
PROX = ["log_prior_dx", "prior_E78", "prior_IHD", "prior_HTN", "prior_DM", "prior_GERD", "prior_gastritis_ulcer"]
mode = sys.argv[1]

if mode == "scan":
    k = sys.argv[2]; c = build_cohort(k).merge(prior, on="eid", how="left")
    r = cox_scan(c, extra=PROX); r["model"] = "M3"; r.to_csv(OUT / f"{k}_M3_all_metabolites.csv", index=False)
    print(k, "M3: Bonferroni-significant", int((r.p_bonf < .05).sum()), "of", int((r.status == "ok").sum()))

elif mode == "decomp":
    M = load_metabolome().set_index("eid"); rows = []
    SETS = {"Model 1": [], "+ hyperlipidaemia (E78)": ["prior_E78"], "+ E78, IHD, hypertension, diabetes": ["prior_E78", "prior_IHD", "prior_HTN", "prior_DM"],
            "+ number of prior diagnoses": ["log_prior_dx"], "+ reflux, gastritis/ulcer": ["prior_GERD", "prior_gastritis_ulcer"], "+ all proxies": PROX}
    for k in ("K63", "K31"):
        c = build_cohort(k).merge(prior, on="eid", how="left")
        for m in ("HDLC", "LDLC", "ApoB", "MetaTG", "GlycoAcet"):
            c[m] = M.reindex(c.eid.values)[m].values
            for name, ex in SETS.items():
                d = pd.concat([design_covars(c, extra=ex), c[[m, "follow_up_days", "event"]]], axis=1).dropna(); d[m] = (d[m] - d[m].mean()) / d[m].std()
                s = CoxPHFitter().fit(d, "follow_up_days", "event").summary.loc[m]
                rows.append({"outcome": k, "metabolite": m, "model": name, "hr": s["exp(coef)"], "lo": s["exp(coef) lower 95%"], "hi": s["exp(coef) upper 95%"], "p": s["p"]})
    D = pd.DataFrame(rows); D.round(4).to_csv(OUT / "decomposition_key_lipids.csv", index=False)
    print(D.pivot_table(index=["outcome", "metabolite"], columns="model", values="hr").round(3)[list(SETS)].to_string())
    # association of the proxies with the exposures
    base = load_base().merge(prior, on="eid"); X = M.reindex(base.eid.values).reset_index(drop=True)
    print("\nmean LDLC / HDLC (RINT z) by prior hyperlipidaemia diagnosis:", {f: round(float(X.loc[base.prior_E78.values == f, "LDLC"].mean()), 3) for f in (0, 1)}, {f: round(float(X.loc[base.prior_E78.values == f, "HDLC"].mean()), 3) for f in (0, 1)})

else:   # defs
    R = load_metabolome(False).set_index("eid"); co = lambda a, b: round(float(pd.concat([a, b], axis=1).dropna().corr().iloc[0, 1]), 4)
    ad = lambda a, b: round(float((a - b).abs().mean()), 4)
    print("MetaTC = TEC + TFC:", co(R.MetaTC, R.TEC + R.TFC), ad(R.MetaTC, R.TEC + R.TFC))
    print("NonHDLnonLDL = MetaTC - HDLC - LDLC:", co(R.NonHDLnonLDL, R.MetaTC - R.HDLC - R.LDLC), ad(R.NonHDLnonLDL, R.MetaTC - R.HDLC - R.LDLC))
    print("TC_HDL = MetaTC - HDLC:", co(R.TC_HDL, R.MetaTC - R.HDLC), ad(R.TC_HDL, R.MetaTC - R.HDLC))
    print("VLDL(-C) + LDLC + HDLC vs MetaTC:", co(R.VLDL + R.LDLC + R.HDLC, R.MetaTC), ad(R.VLDL + R.LDLC + R.HDLC, R.MetaTC))
    for sub in ("VVLDL", "LVLDL", "SLDL", "LHDL"):
        print(f"{sub}: Cho = CE + FC:", ad(R[f"Cho_{sub}"], R[f"CE_{sub}"] + R[f"FC_{sub}"]), "| Lip = Pho + CE + FC + TG:", ad(R[f"Lip_{sub}" if f"Lip_{sub}" in R else f"Lip__{sub}"], R[f"Pho_{sub}"] + R[f"CE_{sub}"] + R[f"FC_{sub}"] + R[f"TG_{sub}"]))
    print("Chy_VLDL series: Cho = CE + FC:", ad(R.Cho_Chy_VLDL, R.CE_Chy_VLDL + R.FC_Chy_VLDL), "| Lip = Pho + CE + FC + X(TC_Chy_VLDL):", ad(R.Lip_Chy_VLDL, R.Pho_Chy_VLDL + R.CE_Chy_VLDL + R.FC_Chy_VLDL + R.TC_Chy_VLDL), "| mean TC_Chy_VLDL", round(float(R.TC_Chy_VLDL.mean()), 4), "mean TG_VVLDL", round(float(R.TG_VVLDL.mean()), 4))
    print("MetaTG vs TG_VLDL+TG_LDL+TG_HDL:", co(R.MetaTG, R.TG_VLDL + R.TG_LDL + R.TG_HDL), ad(R.MetaTG, R.TG_VLDL + R.TG_LDL + R.TG_HDL))
