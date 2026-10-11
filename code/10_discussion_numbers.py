"""Numbers quoted in the Discussion: family-level direction of effects and prevalent-vs-incident agreement."""
import sys; sys.path.insert(0, r"D:\research\Xingang Li\polyp_metabolites\code")
from common import *
from circular_barplot import layout
for k in ["K63", "K31"]:
    p = pd.read_csv(RES / f"step01_primary/model1_{k}_cox_all_metabolites.csv").query("status=='ok'")
    L = layout(list(p.metabolite)).merge(p, on="metabolite")
    sub = L[L.family != ""]
    t = sub.groupby(["group", "family"]).agg(n=("coef", "size"), mean_hr=("hr", "mean"), n_sig=("p_bonf", lambda x: int((x < .05).sum()))).round(3)
    print("\n==", k, "family-level mean HR per SD (subclass measures)"); print(t.loc[["Cholesterol", "Esterified cholesterol", "Free cholesterol", "Triglycerides", "Phospholipids", "Lipoprotein particles", "Triglycerides (%)", "Cholesterol (%)"]].to_string())
    v = pd.read_csv(RES / f"step03b_prevalent/{k}_prevalent_logistic.csv").merge(p, on="metabolite")
    print(k, "prevalent vs incident: r(ln OR, ln HR) =", round(np.log(v["or"]).corr(v.coef), 3), "| same direction among incident hits:", round(((np.sign(np.log(v["or"])) == np.sign(v.coef))[v.p_bonf_y < .05]).mean(), 3))
    fa = L[L.group == "Fatty acids"].set_index("label")
    print(k, "fatty acids:", {l: (round(fa.loc[l, "hr"], 3), "sig" if fa.loc[l, "p_bonf"] < .05 else "ns") for l in ["LA/FA", "PUFA/FA", "FAw6/FA", "SFA/FA", "MUFA/FA", "DHA/FA", "FAw3/FA"]})
    for m in ["GlycoAcet", "HDLC", "LDLC", "ApoB", "ApoA1", "MetaTG", "Glucose", "Creatinine", "Albumin", "Concen_BCAA"]:
        r = p[p.metabolite == m].iloc[0]; print("  ", m, round(r.hr, 3), f"({r.hr_lower95:.3f}-{r.hr_upper95:.3f})", "p_bonf", f"{r.p_bonf:.2g}")
