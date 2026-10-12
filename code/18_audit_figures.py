"""Supplementary Fig S9 and Table S6: audit of the lipoprotein results (NMR vs clinical chemistry, adjustment for prior-diagnosis proxies, crude rates)."""
import sys; sys.path.insert(0, r"D:\research\Xingang Li\polyp_metabolites\code")
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import *

plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"], "font.size": 8, "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42, "savefig.dpi": 300})
UP, DN, GR = "#D55E00", "#0072B2", "#999999"
LAB = {"K63": "Colon polyp (K63.5)", "K31": "Stomach/duodenal polyp (K31.7)"}
A = RES / "step16_lipid_audit"; B = RES / "step17_proxy_adjustment"
al = pd.read_csv(A / "A_alignment_nmr_vs_clinical.csv"); hr = pd.read_csv(A / "D_G_hazard_ratios.csv"); crude = pd.read_csv(A / "C_crude_rates_by_hdl_tertile.csv"); dec = pd.read_csv(B / "decomposition_key_lipids.csv")

# ---------------------------------------------------------------- Table S6 (long format)
rows = []
for _, r in al.iterrows(): rows.append({"analysis": "alignment: NMR vs clinical chemistry", "outcome": "", "exposure": f"{r.nmr} vs {r.clinical}", "model": "", "n": r.n, "estimate": r.spearman_r, "lo": None, "hi": None, "note": "Spearman r"})
for _, r in hr.iterrows(): rows.append({"analysis": "hazard ratio per SD", "outcome": LAB[r.outcome], "exposure": r.exposure, "model": r.model, "n": r.n, "estimate": r.hr, "lo": r.lo, "hi": r.hi, "note": f"events {int(r.events)}; P {r.p:.1e}"})
for _, r in dec.iterrows(): rows.append({"analysis": "hazard ratio per SD, stepwise proxies", "outcome": LAB[r.outcome], "exposure": r.metabolite, "model": r.model, "n": None, "estimate": r.hr, "lo": r.lo, "hi": r.hi, "note": f"P {r.p:.1e}"})
for _, r in crude.iterrows(): rows.append({"analysis": "crude incidence per 1000 person-years", "outcome": LAB[r.outcome], "exposure": f"HDL cholesterol ({r.source})", "model": f"{'female' if r.Sex == 0 else 'male'}, {r.t}", "n": r.n, "estimate": r.rate_per_1000py, "lo": None, "hi": None, "note": f"events {int(r.events)}"})
T = pd.DataFrame(rows); T.round(4).to_csv(A / "TableS6_lipid_audit.csv", index=False)

# ---------------------------------------------------------------- Fig S9
fig, ax = plt.subplots(2, 3, figsize=(7.8, 6.6))
pairs = [("HDL cholesterol", "nmr_HDLC", "clin_HDL"), ("LDL cholesterol", "nmr_LDLC", "clin_LDL"), ("Apolipoprotein A1", "nmr_ApoA1", "clin_ApoA"), ("Apolipoprotein B", "nmr_ApoB", "clin_ApoB"),
         ("Triglycerides", "nmr_MetaTG", "clin_TG"), ("Total cholesterol", "nmr_MetaTC", "clin_Cholesterol")]
for j, k in enumerate(LAB):
    a = ax[0, j]; d = hr[(hr.outcome == k) & (hr.model == "Model 1")].set_index("exposure"); y = np.arange(len(pairs))[::-1]
    for i, (nm, n_, c_) in enumerate(pairs):
        for off, e, col, lab in ((0.16, n_, DN, "NMR"), (-0.16, c_, UP, "Clinical chemistry")):
            a.errorbar(d.loc[e, "hr"], y[i] + off, xerr=[[d.loc[e, "hr"] - d.loc[e, "lo"]], [d.loc[e, "hi"] - d.loc[e, "hr"]]], fmt="o", ms=3.5, c=col, lw=1, label=lab if i == 0 else None)
    a.axvline(1, c=GR, lw=.6); a.set_yticks(y, [p[0] for p in pairs]); a.set_xlabel("HR per SD (95% CI), Model 1"); a.set_title(LAB[k], fontsize=8)
    if j == 0: a.legend(frameon=False, fontsize=6.5, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2)
a = ax[0, 2]; c_ = crude[(crude.outcome == "K31") & (crude.source == "clinical")]
for sx, col, nm in ((0, UP, "Women"), (1, DN, "Men")):
    d = c_[c_.Sex == sx]; a.plot(range(3), d.rate_per_1000py, "o-", c=col, label=nm)
a.set_xticks(range(3), ["T1 (low)", "T2", "T3 (high)"]); a.set_xlabel("Sex-specific tertile of HDL cholesterol"); a.set_ylabel("Incidence per 1000 person-years"); a.set_title("Stomach/duodenal polyp, crude rates", fontsize=8); a.legend(frameon=False, fontsize=6.5)
SETS = ["Model 1", "+ hyperlipidaemia (E78)", "+ E78, IHD, hypertension, diabetes", "+ number of prior diagnoses", "+ reflux, gastritis/ulcer", "+ all proxies"]
SHORT = ["Model 1", "+ hyperlipidaemia", "+ E78, IHD, HTN, DM", "+ prior diagnoses (n)", "+ reflux, gastritis", "+ all proxies"]
for i, (mets, ttl) in enumerate((([("HDLC", "HDL-C")], "HDL cholesterol"), ([("LDLC", "LDL-C"), ("ApoB", "Apo B")], "LDL cholesterol and apolipoprotein B"))):
    a = ax[1, i]; y = np.arange(len(SETS))[::-1]; n_ = 0
    for k, col in (("K63", UP), ("K31", DN)):
        for m, ml in mets:
            d = dec[(dec.outcome == k) & (dec.metabolite == m)].set_index("model").loc[SETS]
            off = (-0.18 if k == "K63" else 0.18) + (0.0 if len(mets) == 1 else (-0.07 if m == "LDLC" else 0.07))
            a.errorbar(d.hr, y + off, xerr=[d.hr - d.lo, d.hi - d.hr], fmt="o" if m != "ApoB" else "s", ms=3, c=col, lw=.9, label=f"{'Colon' if k == 'K63' else 'Stomach/duod.'} {ml}")
    a.axvline(1, c=GR, lw=.6); a.set_yticks(y, SHORT, fontsize=6.5); a.set_xlabel("HR per SD (95% CI)"); a.set_title(ttl, fontsize=8); a.legend(frameon=False, fontsize=6, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2)
a = ax[1, 2]
for m, ml, mk in (("MetaTG", "Total triglycerides", "o"), ("GlycoAcet", "Glycoprotein acetyls", "s")):
    for k, col in (("K63", UP), ("K31", DN)):
        d = dec[(dec.outcome == k) & (dec.metabolite == m)].set_index("model").loc[SETS]
        a.errorbar(d.hr, y + (-0.18 if k == "K63" else 0.18) + (0.0 if m == "MetaTG" else 0.07), xerr=[d.hr - d.lo, d.hi - d.hr], fmt=mk, ms=3, c=col, lw=.9, label=f"{'Colon' if k == 'K63' else 'Stomach/duod.'} {ml.split()[0]}")
a.axvline(1, c=GR, lw=.6); a.set_yticks(y, SHORT, fontsize=6.5); a.set_xlabel("HR per SD (95% CI)"); a.set_title("Triglycerides and glycoprotein acetyls", fontsize=8); a.legend(frameon=False, fontsize=6, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2)
fig.tight_layout(); (FIG / "supplement").mkdir(exist_ok=True)
fig.savefig(FIG / "supplement" / "FigS9_lipid_audit.png", bbox_inches="tight"); fig.savefig(FIG / "supplement" / "FigS9_lipid_audit.pdf", bbox_inches="tight"); plt.close(fig)
print("FigS9 and TableS6 written;", len(T), "rows")
