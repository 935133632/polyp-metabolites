"""Audit of the lipoprotein results (HDL, LDL, triglycerides): data alignment, variable definitions, independent clinical-chemistry replication,
crude incidence rates, mutual adjustment and proxies for lipid-lowering treatment / health-care contact.
usage: python 16_lipid_audit.py            (runs everything; prior-diagnosis proxies are cached in processed/prior_dx.pkl)"""
import sys; sys.path.insert(0, r"D:\research\Xingang Li\polyp_metabolites\code")
import re, warnings; warnings.filterwarnings("ignore")
from scipy.stats import spearmanr
from lifelines import CoxPHFitter
from common import *

OUT = RES / "step16_lipid_audit"; OUT.mkdir(exist_ok=True)
CLIN = {"HDL": "p30760_i0", "LDL": "p30780_i0", "ApoA": "p30630_i0", "ApoB": "p30640_i0", "TG": "p30870_i0", "Cholesterol": "p30690_i0",
        "Glucose": "p30740_i0", "Albumin": "p30600_i0", "Creatinine": "p30700_i0", "CRP": "p30710_i0"}
NMR_MATCH = {"HDL": "HDLC", "LDL": "LDLC", "ApoA": "ApoA1", "ApoB": "ApoB", "TG": "MetaTG", "Cholesterol": "MetaTC", "Glucose": "Glucose", "Albumin": "Albumin", "Creatinine": "Creatinine"}

# ---------------------------------------------------------------- prior diagnoses (proxies for treatment indication and health-care contact)
f = PROC / "prior_dx.pkl"
if not f.exists():
    base = load_base().set_index("eid"); eids = set(base.index); bd = base.baseline_date
    dcols_a = [f"p41280_a{i}" for i in range(198)]; dcols_b = [f"p41280_a{i}" for i in range(198, 259)]
    PB = pd.read_csv(DATA / "WH-Outcome-Diagnosis-ICD10-PartB-20251128.csv", usecols=["eid", *dcols_b], low_memory=False).set_index("eid")
    FLAGS = {"prior_E78": r"^E78", "prior_IHD": r"^I2[0-5]", "prior_HTN": r"^I1[0-5]", "prior_DM": r"^E1[0-4]", "prior_GERD": r"^K2[0-2]", "prior_gastritis_ulcer": r"^K(2[5-9])"}
    rows = {}
    for ca in pd.read_csv(DATA / "WH-Outcome-Diagnosis-ICD10-PartA-20251128.csv", usecols=["eid", "p41270", *dcols_a], chunksize=20000, low_memory=False):
        ca = ca[ca.eid.isin(eids)].reset_index(drop=True); cb = PB.reindex(ca.eid.values).reset_index(drop=True)
        D = pd.concat([ca[dcols_a], cb[dcols_b]], axis=1).apply(pd.to_datetime, errors="coerce").to_numpy()
        for i, (e, s) in enumerate(zip(ca.eid.values, ca.p41270.values)):
            r = {"n_prior_dx": 0, **{k: 0 for k in FLAGS}}
            if isinstance(s, str):
                seen = set(); b = bd.loc[e]
                for j, code in enumerate(s.split("|")):
                    dt = D[i, j] if j < D.shape[1] else None
                    if dt is None or pd.isna(dt) or dt > b: continue
                    c3 = code.strip()[:3]; seen.add(c3)
                    for k, pat in FLAGS.items():
                        if re.match(pat, c3): r[k] = 1
                r["n_prior_dx"] = len(seen)
            rows[int(e)] = r
    pd.DataFrame.from_dict(rows, orient="index").rename_axis("eid").reset_index().to_pickle(f)
prior = pd.read_pickle(f); prior["log_prior_dx"] = np.log1p(prior.n_prior_dx)
print("prior diagnoses:", prior.drop(columns="eid").mean().round(3).to_dict())

clin = pd.read_csv(DATA / "FRS-Biochemistry-Inst0-20251127.csv", usecols=["eid", *CLIN.values()]).rename(columns={v: k for k, v in CLIN.items()})
M = load_metabolome().set_index("eid"); Mraw = load_metabolome(False).set_index("eid")
base = load_base().merge(clin, on="eid", how="left").merge(prior, on="eid", how="left")

# ---------------------------------------------------------------- A. alignment of NMR and clinical chemistry (same participants, same blood draw)
rows = []
for cn, nn in NMR_MATCH.items():
    d = pd.DataFrame({"clin": base.set_index("eid")[cn], "nmr": Mraw[nn]}).dropna()
    r = spearmanr(d.clin, d.nmr)[0]; rows.append({"clinical": cn, "nmr": nn, "n": len(d), "spearman_r": round(r, 3)})
al = pd.DataFrame(rows); al.to_csv(OUT / "A_alignment_nmr_vs_clinical.csv", index=False); print("\nA. NMR vs clinical chemistry\n", al.to_string(index=False))
sx = base.set_index("eid").join(Mraw[["HDLC", "LDLC", "MetaTG"]]).groupby("Sex")[["HDL", "HDLC", "LDL", "LDLC", "Creatinine"]].mean().round(2); print("\nmean by sex (0 female, 1 male):\n", sx)
chk = pd.DataFrame({"TC_HDL": Mraw.TC_HDL, "MetaTC_minus_HDLC": Mraw.MetaTC - Mraw.HDLC}).dropna(); print("\nB. TC_HDL vs MetaTC - HDLC: r =", round(chk.corr().iloc[0, 1], 4), "| mean abs diff", round((chk.TC_HDL - chk.MetaTC_minus_HDLC).abs().mean(), 4))
print("LDLC vs CLDL r =", round(Mraw[["LDLC", "CLDL"]].dropna().corr().iloc[0, 1], 3))

def fit(df, exp, extra=(), outcome_cols=("follow_up_days", "event")):
    X = design_covars(df, extra=extra); d = pd.concat([X, df[[exp, *outcome_cols]]], axis=1).replace([np.inf, -np.inf], np.nan).dropna()
    d[exp] = (d[exp] - d[exp].mean()) / d[exp].std(); d = d.loc[:, d.nunique() > 1]
    s = CoxPHFitter().fit(d, "follow_up_days", "event").summary.loc[exp]
    return {"n": len(d), "events": int(d.event.sum()), "hr": s["exp(coef)"], "lo": s["exp(coef) lower 95%"], "hi": s["exp(coef) upper 95%"], "p": s["p"]}

res = []
for k in ("K63", "K31"):
    c = build_cohort(k).merge(clin, on="eid", how="left").merge(prior, on="eid", how="left")
    for col in ("HDL", "LDL", "ApoA", "ApoB", "TG", "Cholesterol"): c["clin_" + col] = rint(c[col])          # clinical chemistry, rank-normalised like the NMR measures
    c["clin_CRP"] = rint(c["CRP"])
    for m in ("HDLC", "LDLC", "CLDL", "ApoA1", "ApoB", "MetaTG", "MetaTC", "TC_HDL", "NonHDLnonLDL", "GlycoAcet", "Linopt"): c["nmr_" + m] = M.reindex(c.eid.values)[m].values
    c["Sex_z"] = c.Sex
    # C. crude incidence rates by sex-specific tertile of HDL cholesterol (clinical and NMR)
    for src, col in (("clinical", "HDL"), ("NMR", "nmr_HDLC")):
        c["t"] = c.groupby("Sex")[col if src == "clinical" else col].transform(lambda x: pd.qcut(x, 3, labels=["T1 low", "T2", "T3 high"]))
        g = c.groupby(["Sex", "t"], observed=True).apply(lambda d: pd.Series({"n": len(d), "events": int(d.event.sum()), "rate_per_1000py": 1000 * d.event.sum() / d.follow_up_years.sum()})).reset_index()
        g.insert(0, "outcome", k); g.insert(1, "source", src); res.append(("C", g))
    # D. independent clinical-chemistry replication and NMR, Model 1
    rowsD = []
    for e in [f"nmr_{m}" for m in ("HDLC", "LDLC", "CLDL", "ApoA1", "ApoB", "MetaTG", "MetaTC", "TC_HDL", "NonHDLnonLDL", "GlycoAcet", "Linopt")] + [f"clin_{m}" for m in ("HDL", "LDL", "ApoA", "ApoB", "TG", "Cholesterol", "CRP")]:
        rowsD.append({"outcome": k, "exposure": e, "model": "Model 1", **fit(c, e)})
    # E. + proxies for lipid-lowering/aspirin indication, health-care contact, reflux/gastritis
    PROX = ["log_prior_dx", "prior_E78", "prior_IHD", "prior_HTN", "prior_DM", "prior_GERD", "prior_gastritis_ulcer"]
    for e in ["nmr_HDLC", "nmr_LDLC", "nmr_ApoA1", "nmr_ApoB", "nmr_MetaTG", "nmr_GlycoAcet", "clin_HDL", "clin_LDL"]:
        rowsD.append({"outcome": k, "exposure": e, "model": "Model 1 + prior-diagnosis proxies", **fit(c, e, extra=PROX)})
    # F. within-sex fits (the constant sex column is dropped inside fit())
    for e in ["nmr_HDLC", "nmr_LDLC", "clin_HDL", "clin_LDL"]:
        for sx_, nm in ((0, "female"), (1, "male")):
            rowsD.append({"outcome": k, "exposure": e, "model": f"Model 1, {nm} only", **fit(c[c.Sex == sx_], e)})
    D = pd.DataFrame(rowsD); res.append(("D", D))
    # G. mutual adjustment: HDL, LDL and triglycerides together
    X = design_covars(c); d = pd.concat([X, c[["nmr_HDLC", "nmr_LDLC", "nmr_MetaTG", "follow_up_days", "event"]]], axis=1).dropna()
    for e in ("nmr_HDLC", "nmr_LDLC", "nmr_MetaTG"): d[e] = (d[e] - d[e].mean()) / d[e].std()
    s = CoxPHFitter().fit(d, "follow_up_days", "event").summary
    Gd = pd.DataFrame([{"outcome": k, "exposure": e, "model": "HDL-C + LDL-C + TG mutually adjusted", "n": len(d), "events": int(d.event.sum()), "hr": s.loc[e, "exp(coef)"], "lo": s.loc[e, "exp(coef) lower 95%"], "hi": s.loc[e, "exp(coef) upper 95%"], "p": s.loc[e, "p"]} for e in ("nmr_HDLC", "nmr_LDLC", "nmr_MetaTG")]); res.append(("G", Gd))

Dall = pd.concat([x for t, x in res if t == "D"] + [x for t, x in res if t == "G"]); Dall.round(4).to_csv(OUT / "D_G_hazard_ratios.csv", index=False)
pd.concat([x for t, x in res if t == "C"]).round(3).to_csv(OUT / "C_crude_rates_by_hdl_tertile.csv", index=False)
pd.set_option("display.width", 200)
print("\nC. crude incidence per 1000 person-years by sex-specific HDL tertile"); print(pd.concat([x for t, x in res if t == "C"]).round(2).to_string(index=False))
print("\nD/G. hazard ratios per SD"); print(Dall.assign(hr=lambda d: d.hr.round(3), lo=lambda d: d.lo.round(3), hi=lambda d: d.hi.round(3), p=lambda d: d.p.map(lambda x: f"{x:.1e}")).to_string(index=False))
