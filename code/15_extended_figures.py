"""Figures and tables for the extended analyses: model adjustment (12), specificity (13), metabolite score and dose-response (14)."""
import sys; sys.path.insert(0, r"D:\research\Xingang Li\polyp_metabolites\code")
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, leaves_list
from common import *
from metab_annot import metab_class
from circular_barplot import circular_heatmap

plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"], "font.size": 8, "axes.spines.top": False,
                     "axes.spines.right": False, "pdf.fonttype": 42, "savefig.dpi": 300})
UP, DN, GR = "#D55E00", "#0072B2", "#BBBBBB"
LAB = {"K63": "Colon polyp (K63.5)", "K31": "Stomach/duodenal polyp (K31.7)"}
ST, SF = TAB / "supplement", FIG / "supplement"; ST.mkdir(parents=True, exist_ok=True); SF.mkdir(parents=True, exist_ok=True)
OUT = RES / "step15_summary"; OUT.mkdir(exist_ok=True)


def save(fig, stem): fig.savefig(f"{stem}.png", bbox_inches="tight"); fig.savefig(f"{stem}.pdf", bbox_inches="tight"); plt.close(fig)


# ---------------------------------------------------------------- 12: model adjustment
adj_rows = []
for k, tag in (("K63", "A_colon_polyp"), ("K31", "B_stomach_duodenal_polyp")):
    R = {m: pd.read_csv(RES / f"step12_model_adjustment/{k}_{m}.csv").query("status == 'ok'") for m in ("M1", "M2a", "M2b", "M2c")}
    if not all(len(v) for v in R.values()): continue
    rings = {"Model 1": R["M1"], "+adiposity": R["M2a"], "+HbA1c, CRP": R["M2b"], "+all four": R["M2c"]}
    circular_heatmap(rings, LAB[k] + ", extended adjustment", str(FIG / f"Fig8{tag}_adjustment_circular_heatmap"), note="same participants in every ring", cbar_label="ln(HR) per 1 SD")
    m1 = R["M1"].set_index("metabolite"); hits = m1[m1.p_bonf < .05]
    for m in ("M2a", "M2b", "M2c"):
        x = R[m].set_index("metabolite").loc[hits.index]
        att = (1 - x.coef / hits.coef) * 100                      # % attenuation of ln(HR), positive = smaller after adjustment
        adj_rows.append({"outcome": k, "model": m, "n": int(x.n_complete.max()), "events": int(x.n_events.max()), "M1_hits": len(hits),
                         "still_bonf_sig": int((x.p_bonf < .05).sum()), "same_direction": int((np.sign(x.coef) == np.sign(hits.coef)).sum()),
                         "median_attenuation_pct": float(att.median()), "iqr_lo": float(att.quantile(.25)), "iqr_hi": float(att.quantile(.75)),
                         "bonf_sig_total": int((R[m].p_bonf < .05).sum())})
if adj_rows:
    pd.DataFrame(adj_rows).round(3).to_csv(OUT / "model_adjustment_summary.csv", index=False); print(pd.DataFrame(adj_rows).round(2).to_string(index=False))

# ---------------------------------------------------------------- 13: specificity
NAMES = {"K20": "Oesophagitis", "K21": "Gastro-oesophageal reflux", "K22": "Other oesophageal diseases", "K25": "Gastric ulcer", "K26": "Duodenal ulcer", "K29": "Gastritis / duodenitis",
         "K30": "Dyspepsia", "K40": "Inguinal hernia", "K42": "Umbilical hernia", "K43": "Ventral hernia", "K44": "Hiatus hernia", "K52": "Other non-infective colitis",
         "K56": "Ileus / obstruction", "K57": "Diverticular disease", "K58": "Irritable bowel syndrome", "K59": "Other functional bowel disorders", "K62": "Other anorectal diseases",
         "K64": "Haemorrhoids", "K66": "Other peritoneal disorders", "K76": "Other liver diseases (incl. fatty liver)", "K80": "Gallstones", "K83": "Other biliary diseases",
         "K92": "Other digestive diseases (incl. GI bleeding)", "CRC": "Colorectal cancer (C18-C20)"}
vec, info = {}, []
for k in LAB:
    p = pd.read_csv(RES / f"step01_primary/model1_{k}_cox_all_metabolites.csv").query("status == 'ok'").set_index("metabolite")
    vec[LAB[k]] = p.coef; info.append({"outcome": LAB[k], "events": int(p.n_events.max()), "bonf_sig": int((p.p_bonf < .05).sum())})
for code, nm in NAMES.items():
    f = RES / f"step13_specificity/{code}.csv"
    if not f.exists(): continue
    p = pd.read_csv(f).query("status == 'ok'").set_index("metabolite"); vec[nm] = p.coef; info.append({"outcome": nm, "events": int(p.n_events.max()), "bonf_sig": int((p.p_bonf < .05).sum())})
if len(vec) > 4:
    V = pd.DataFrame(vec); I = pd.DataFrame(info).set_index("outcome")
    ref = list(LAB.values()); tab = pd.DataFrame({r: V.corrwith(V[r]) for r in ref}); tab = tab.join(I); tab.round(3).to_csv(OUT / "specificity_profile_correlation.csv")
    others = [c for c in V.columns if c not in ref]
    order = tab.loc[others].sort_values(ref[0]).index.tolist()
    fig, ax = plt.subplots(1, 2, figsize=(7.4, 5.6), gridspec_kw={"width_ratios": [1, 1.25]})
    y = np.arange(len(order)); ax[0].barh(y + .2, tab.loc[order, ref[0]], height=.4, color=UP, label="vs colon polyp"); ax[0].barh(y - .2, tab.loc[order, ref[1]], height=.4, color=DN, label="vs stomach/duodenal polyp")
    ax[0].set_yticks(y, [f"{o} (n={int(tab.loc[o, 'events']):,})" for o in order], fontsize=6.5); ax[0].axvline(0, c="k", lw=.5); ax[0].set_xlabel("Correlation of 251 ln(HR) values with the polyp profile"); ax[0].legend(frameon=False, fontsize=6.5, loc="lower center", bbox_to_anchor=(.5, 1.0), ncol=2)
    cls = pd.Series({m: metab_class(m) for m in V.index})
    Hc = V.groupby(cls).mean().T                                         # outcomes x classes: mean ln(HR)
    ordo = leaves_list(linkage(Hc.values, "average")); Hc = Hc.iloc[ordo]
    im = ax[1].imshow(Hc.values, cmap="RdBu_r", vmin=-.15, vmax=.15, aspect="auto"); ax[1].set_yticks(range(len(Hc)), Hc.index, fontsize=6.5)
    ax[1].set_xticks(range(Hc.shape[1]), [c if len(c) < 24 else c[:22] + "..." for c in Hc.columns], rotation=60, ha="right", fontsize=6); ax[1].set_title("Mean ln(HR) by metabolite class (outcomes clustered)", fontsize=7.5)
    fig.colorbar(im, ax=ax[1], shrink=.6, label="mean ln(HR) per 1 SD"); fig.tight_layout(); save(fig, str(FIG / "Fig9_specificity"))
    print(tab.round(2).sort_values(ref[0], ascending=False).to_string())

# ---------------------------------------------------------------- 14: score and dose-response
D = RES / "step14_score_dose_response"
fig, ax = plt.subplots(2, 2, figsize=(7.4, 5.6)); ok = False
for j, k in enumerate(LAB):
    fq, fc = D / f"{k}_score_quartiles.csv", D / f"{k}_score_cumincidence.csv"
    if not (fq.exists() and fc.exists()): continue
    ok = True; q = pd.read_csv(fq); ci = pd.read_csv(fc)
    for qq, col in zip(["Q1", "Q2", "Q3", "Q4"], ["#4C78A8", "#9DB7D5", "#F2A57A", UP]):
        d_ = ci[ci.quartile == qq]; ax[0, j].step(d_.years, 100 * d_.cum_incidence, where="post", color=col, lw=1.2, label=qq)
    ax[0, j].set_xlabel("Years since baseline"); ax[0, j].set_ylabel("Cumulative incidence (%)"); ax[0, j].legend(title="Score quartile", frameon=False, fontsize=6.5); ax[0, j].set_title(LAB[k], fontsize=8)
    rows = q[q.term.isin(["Q2 vs Q1", "Q3 vs Q1", "Q4 vs Q1"])]; y = np.arange(len(rows))[::-1]
    ax[1, j].errorbar(rows.hr, y, xerr=[rows.hr - rows.lower95, rows.upper95 - rows.hr], fmt="s", c="k", ms=4, lw=1); ax[1, j].axvline(1, c=GR, lw=.7)
    ax[1, j].set_yticks(y, [t.replace(" vs Q1", "") for t in rows.term]); ax[1, j].set_xlabel("HR vs Q1 (95% CI)")
    sd_ = q[q.term == "per 1 SD"].iloc[0]; pt = q[q.term.str.startswith("P-trend")].iloc[0]
    ax[1, j].set_title(f"Per SD: HR {sd_.hr:.2f} ({sd_.lower95:.2f}-{sd_.upper95:.2f}); P-trend {pt.p:.1e}", fontsize=7.5)
if ok: fig.tight_layout(); save(fig, str(FIG / "Fig10_metabolite_score"))
for k in LAB:
    f = D / f"{k}_dose_response.csv"
    if not f.exists(): continue
    dr = pd.read_csv(f); nl = pd.read_csv(D / f"{k}_nonlinearity.csv").set_index("metabolite"); ms = list(dict.fromkeys(dr.metabolite)); n = len(ms)
    fig, ax = plt.subplots(3, 3, figsize=(7.4, 6.2), sharex=True)
    for a, m in zip(ax.ravel(), ms):
        d_ = dr[dr.metabolite == m]; a.fill_between(d_.z, d_.lower95, d_.upper95, color="#9DB7D5", alpha=.7, lw=0); a.plot(d_.z, d_.hr, c="k", lw=1); a.axhline(1, c=GR, lw=.6)
        a.set_title(f"{m}, P-nonlinear {nl.loc[m, 'p_nonlinear']:.2g}", fontsize=7); a.set_ylabel("HR vs median")
    for a in ax[-1]: a.set_xlabel("Metabolite (SD, normalised)")
    for a in ax.ravel()[n:]: a.axis("off")
    fig.suptitle(LAB[k], fontsize=9); fig.tight_layout(); save(fig, str(SF / f"FigS8_dose_response_{k}"))
    print(k, "nonlinearity:", nl.p_nonlinear.round(4).to_dict())
