"""Circular bar plot of per-metabolite hazard ratios (all 251 NMR metabolites), grouped like a Nightingale metabolomics figure.
Bars start at HR = 1 (red = higher risk, blue = lower risk); black lines = 95% CI; stars = Bonferroni-adjusted P."""
import re
import numpy as np, pandas as pd
import matplotlib.pyplot as plt

SUB = {"Chy_VLDL": ("XXL", "VLDL"), "ELVLDL": ("XXL", "VLDL"), "VVLDL": ("XL", "VLDL"), "LVLDL": ("L", "VLDL"), "MVLDL": ("M", "VLDL"), "SVLDL": ("S", "VLDL"),
       "VSVLDL": ("XS", "VLDL"), "IDL": ("IDL", "IDL"), "LLDL": ("L", "LDL"), "MLDL": ("M", "LDL"), "SLDL": ("S", "LDL"),
       "VLHDL": ("XL", "HDL"), "LHDL": ("L", "HDL"), "MHDL": ("M", "HDL"), "SHDL": ("S", "HDL")}
SUB_ORDER = ["Chy_VLDL", "ELVLDL", "VVLDL", "LVLDL", "MVLDL", "SVLDL", "VSVLDL", "IDL", "LLDL", "MLDL", "SLDL", "VLHDL", "LHDL", "MHDL", "SHDL"]
PREFIX_GROUP = {"Concen": "Lipoprotein particles", "Lip": "Total lipids", "Pho": "Phospholipids", "Cho": "Cholesterol", "CE": "Esterified cholesterol",
                "FC": "Free cholesterol", "TG": "Triglycerides"}
PT_GROUP = {"Pho": "Phospholipids (%)", "Cho": "Cholesterol (%)", "CE": "Esterified cholesterol (%)", "FC": "Free cholesterol (%)", "TG": "Triglycerides (%)"}
TOTALS = {  # group -> measures shown after the subclass bars, with display names
    "Lipoprotein particles": {"Concen_Lipo": "Total particles", "Concen_VLDL": "VLDL particles", "Concen_LDL": "LDL particles", "Concen_HDL": "HDL particles"},
    "Total lipids": {"LipLIpo": "Total lipids", "LipVLDL": "Lipids in VLDL", "LipLDL": "Lipids in LDL", "LipHDL": "Lipids in HDL"},
    "Phospholipids": {"Phospho": "Phospholipids", "PhoLipo": "Phospholipids in lipoproteins", "PhoVLDL": "Phospholipids in VLDL", "PhoLDL": "Phospholipids in LDL",
                      "PhoHDL": "Phospholipids in HDL", "Phospha": "Phosphatidylcholines", "Sphingomyelins": "Sphingomyelins", "Cholines": "Total cholines"},
    "Cholesterol": {"MetaTC": "Total cholesterol", "LDLC": "LDL-C", "CLDL": "Clinical LDL-C", "HDLC": "HDL-C", "VLDL": "VLDL-C", "NonHDLnonLDL": "Remnant cholesterol", "TC_HDL": "TC in HDL"},
    "Esterified cholesterol": {"TEC": "Total esterified chol.", "CE_VLDL": "CE in VLDL", "CE_LDL": "CE in LDL", "CE_HDL": "CE in HDL"},
    "Free cholesterol": {"TFC": "Total free chol.", "FC_VLDL": "FC in VLDL", "FC_LDL": "FC in LDL", "FC_HDL": "FC in HDL"},
    "Triglycerides": {"MetaTG": "Total triglycerides", "TG_VLDL": "TG in VLDL", "TG_LDL": "TG in LDL", "TG_HDL": "TG in HDL"},
}
OTHER = {"Diame_VLDL": "VLDL size", "Diame_LDL": "LDL size", "Diame_HDL": "HDL size", "ApoB": "Apo-B", "ApoA1": "Apo-A1", "ApoRatio": "Apo-B/Apo-A1",
         "TPR": "TPR", "GluLac": "Glucose-lactate", "SpeAla": "Spectrometer Ala"}
FA = {"FattyAcids": "Total FA", "Unsaturation": "Unsaturation", "Omega3": "FAw3", "Omega6": "FAw6", "Polyunsaturated": "PUFA", "Monounsaturated": "MUFA", "Saturated": "SFA",
      "Linoleic": "LA", "Docosahexaenoic": "DHA", "O3pt": "FAw3/FA", "O6pt": "FAw6/FA", "Polypt": "PUFA/FA", "Monopt": "MUFA/FA", "Satupt": "SFA/FA", "Linopt": "LA/FA",
      "Docopt": "DHA/FA", "PolyRatio": "PUFA/MUFA", "O6Ratio": "FAw6/FAw3"}
SMALL = {"Alanine": "Alanine", "Glutamine": "Glutamine", "Glycine": "Glycine", "Histidine": "Histidine", "Concen_BCAA": "BCAA", "Isoleucine": "Isoleucine", "Leucine": "Leucine",
         "Valine": "Valine", "Phenylalanine": "Phenylalanine", "Tyrosine": "Tyrosine", "Glucose": "Glucose", "Lactate": "Lactate", "Pyruvate": "Pyruvate", "Citrate": "Citrate",
         "3Hydroxybutyrate": "3-Hydroxybutyrate", "Acetate": "Acetate", "Acetoacetate": "Acetoacetate", "Acetone": "Acetone", "Creatinine": "Creatinine", "Albumin": "Albumin",
         "GlycoAcet": "Glycoprotein acetyls"}
GROUP_ORDER = ["Lipoprotein particles", "Total lipids", "Phospholipids", "Cholesterol", "Esterified cholesterol", "Free cholesterol", "Triglycerides",
               "Phospholipids (%)", "Cholesterol (%)", "Esterified cholesterol (%)", "Free cholesterol (%)", "Triglycerides (%)",
               "Sizes, apolipoproteins & other", "Fatty acids", "Amino acids, glycolysis, ketone bodies & fluid balance"]


def layout(names):
    """Return DataFrame(metabolite, group, label, family, order) covering every name exactly once."""
    tot = {m for d in TOTALS.values() for m in d}
    rows = []
    for n in names:
        m = re.fullmatch(r"(Pho|Cho|CE|FC|TG)_+(\w+?)_pt", n)
        if m and m.group(2) in SUB:
            lab, fam = SUB[m.group(2)]; rows.append((n, PT_GROUP[m.group(1)], lab, fam, SUB_ORDER.index(m.group(2)))); continue
        if n == "Chy_VLDL": rows.append((n, "Lipoprotein particles", "XXL", "VLDL", 0)); continue
        if n == "TC_Chy_VLDL":   # 7th item of the Chy series; assumed to be triglycerides in chylomicrons / XXL-VLDL
            rows.append((n, "Triglycerides", "XXL", "VLDL", 0)); continue
        m = re.fullmatch(r"(Concen|Lip|Pho|Cho|CE|FC|TG)_+(\w+)", n)
        if m and m.group(2) in SUB and n not in tot:
            lab, fam = SUB[m.group(2)]; rows.append((n, PREFIX_GROUP[m.group(1)], lab, fam, SUB_ORDER.index(m.group(2)))); continue
        hit = next(((g, d) for g, d in TOTALS.items() if n in d), None)
        if hit: rows.append((n, hit[0], hit[1][n], "", 100 + list(hit[1]).index(n))); continue
        if n in OTHER: rows.append((n, "Sizes, apolipoproteins & other", OTHER[n], "", list(OTHER).index(n))); continue
        if n in FA: rows.append((n, "Fatty acids", FA[n], "", list(FA).index(n))); continue
        if n in SMALL: rows.append((n, "Amino acids, glycolysis, ketone bodies & fluid balance", SMALL[n], "", list(SMALL).index(n))); continue
        raise KeyError(f"unassigned metabolite: {n}")
    L = pd.DataFrame(rows, columns=["metabolite", "group", "label", "family", "order"])
    assert L.metabolite.is_unique and set(L.metabolite) == set(names), "layout does not cover the metabolites exactly once"
    L["g"] = L.group.map({g: i for i, g in enumerate(GROUP_ORDER)}); return L.sort_values(["g", "order"]).reset_index(drop=True)


def stars(p):
    return "***" if p < .001 else "**" if p < .01 else "*" if p < .05 else ""


def circular_barplot(res, title, outfile_stem, lo=0.80, hi=1.20, size=10.5):
    """res: Cox result table with metabolite, hr, hr_lower95, hr_upper95, p_bonf."""
    L = layout(list(res.metabolite)).merge(res[["metabolite", "hr", "hr_lower95", "hr_upper95", "p_bonf"]], on="metabolite")
    n = len(L); ng = L.group.nunique(); gap = np.radians(1.6); top_gap = np.radians(11)
    w = (2 * np.pi - ng * gap - top_gap) / n
    th = np.zeros(n); t = top_gap / 2; prev = None
    for i, g in enumerate(L.group):
        if prev is not None and g != prev: t += gap
        th[i] = t + w / 2; t += w; prev = g
    L["theta"] = th
    r_in, r_out = .26, .60
    rad = lambda v: r_in + (np.clip(v, lo, hi) - lo) / (hi - lo) * (r_out - r_in)
    fig = plt.figure(figsize=(size, size)); ax = fig.add_axes([.02, .02, .96, .96], projection="polar")
    ax.set_theta_zero_location("N"); ax.set_theta_direction(-1); ax.set_ylim(0, 1.0); ax.axis("off")
    # shaded sectors, alternating
    for k, (g, d) in enumerate(L.groupby("g", sort=True)):
        if k % 2 == 0:
            a0, a1 = d.theta.min() - w / 2, d.theta.max() + w / 2
            ax.bar((a0 + a1) / 2, r_out - r_in + .04, width=a1 - a0, bottom=r_in - .02, color="#EBEBEB", zorder=0, lw=0)
    gg = np.linspace(top_gap / 2, 2 * np.pi - top_gap / 2, 600)   # reference rings stop at the top gap, where the scale labels sit
    for v in (.8, .9, 1.1, 1.2): ax.plot(gg, np.full_like(gg, rad(v)), c="#9A9A9A", lw=.35, zorder=1)
    ax.plot(gg, np.full_like(gg, rad(1.0)), c="k", lw=.7, zorder=4)
    for v in (.8, .9, 1.0, 1.1, 1.2): ax.text(0, rad(v), f"{v:.2f}", fontsize=6.5, ha="center", va="center", zorder=6)
    for _, r in L.iterrows():
        sig = r.p_bonf < .05; col = "#B22A2A" if r.hr > 1 else "#2D5DA8"
        ax.bar(r.theta, rad(r.hr) - rad(1.0), width=w * .78, bottom=rad(1.0), color=col, alpha=1 if sig else .35, lw=0, zorder=3)
        ax.plot([r.theta, r.theta], [rad(r.hr_lower95), rad(r.hr_upper95)], c="k", lw=.35, zorder=5)
    # labels, stars
    fs = 3.9
    for _, r in L.iterrows():
        deg = np.degrees(r.theta); phi = 90 - deg; right = np.cos(np.radians(phi)) >= 0
        rot = phi if right else phi + 180; ha = "left" if right else "right"
        s = stars(r.p_bonf); r_star = .615
        if s: ax.text(r.theta, r_star, s, fontsize=4.5, ha="center", va="center", rotation=rot, rotation_mode="anchor", fontweight="bold", zorder=6)
        ax.text(r.theta, .65, r.label, fontsize=fs, ha=ha, va="center", rotation=rot, rotation_mode="anchor", zorder=6)
    # subclass-family labels (VLDL/IDL/LDL/HDL) beyond the short subclass labels
    for g, d in L[L.family != ""].groupby("group"):
        for fam, f in d.groupby("family", sort=False):
            tm = f.theta.mean(); deg = np.degrees(tm); phi = 90 - deg; right = np.cos(np.radians(phi)) >= 0
            ax.text(tm, .70, fam, fontsize=5, fontweight="bold", ha="left" if right else "right", va="center", rotation=phi if right else phi + 180, rotation_mode="anchor", zorder=6)
    # group names outside
    for g, d in L.groupby("group"):
        tm = (d.theta.min() + d.theta.max()) / 2; deg = np.degrees(tm) % 360
        rot = -deg if (deg < 90 or deg > 270) else -deg + 180
        ax.text(tm, .90 if len(g) < 28 else .90, g.replace("Amino acids, glycolysis, ketone bodies & fluid balance", "Amino acids, glycolysis,\nketone bodies & fluid balance").replace("Sizes, apolipoproteins & other", "Sizes, ApoLP & other"),
                fontsize=8, fontweight="bold", ha="center", va="center", rotation=rot, rotation_mode="anchor", zorder=6)
    ax.text(0, 0, title + "\nHR per 1 SD\n* P$_{Bonf}$<0.05  ** <0.01  *** <0.001\nred: higher risk, blue: lower risk\nfaded bar: not significant", fontsize=7.5, ha="center", va="center", linespacing=1.5)
    fig.savefig(f"{outfile_stem}.png", dpi=300); fig.savefig(f"{outfile_stem}.pdf"); plt.close(fig)
    return L
