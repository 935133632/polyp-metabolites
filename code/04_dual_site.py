"""Dual-site overlap: participants eligible for both K31 and K63 incident analyses; metabolites separating dual from single-site cases."""
import sys; sys.path.insert(0, r"D:\research\Xingang Li\polyp_metabolites\code")
import statsmodels.api as sm, warnings; warnings.filterwarnings("ignore")
from common import *
d = RES / "step07_dual_site"; d.mkdir(exist_ok=True)
a = build_cohort("K31"); b = build_cohort("K63")
m = a.merge(b[["eid", "event", "event_date", "follow_up_days"]], on="eid", suffixes=("_k31", "_k63"))
m["group"] = np.select([(m.event_k31 == 1) & (m.event_k63 == 1), (m.event_k31 == 1), (m.event_k63 == 1)], ["both", "K31_only", "K63_only"], "none")
def timing(r):
    if r.group != "both": return ""
    return "same_day" if r.event_date_k31 == r.event_date_k63 else ("K31_first" if r.event_date_k31 < r.event_date_k63 else "K63_first")
m["timing"] = m.apply(timing, axis=1)
m.group.value_counts().rename_axis("group").reset_index(name="n").to_csv(d / "group_counts.csv", index=False)
m[m.group == "both"].timing.value_counts().rename_axis("timing").reset_index(name="n").to_csv(d / "timing_counts.csv", index=False)
p1 = pd.read_csv(RES / "step01_primary/model1_K31_cox_all_metabolites.csv"); p2 = pd.read_csv(RES / "step01_primary/model1_K63_cox_all_metabolites.csv")
j = p1.merge(p2, on="metabolite", suffixes=("_k31", "_k63"))
j.to_csv(d / "cross_outcome_cox_comparison.csv", index=False)
shared = j[(j.p_bonf_k31 < .05) & (j.p_bonf_k63 < .05) & (np.sign(j.coef_k31) == np.sign(j.coef_k63))]
shared_set = set(shared.metabolite)
sub = m[m.group.isin(["both", "K31_only", "K63_only"])].copy(); sub["y"] = (sub.group == "both").astype(int)
M = load_metabolome().set_index("eid"); X = design_covars(sub); Mc = M.reindex(sub.eid.values); rows = []
for t in M.columns:   # all 251 metabolites
    df = X.copy(); df["x"] = Mc[t].values; df["y"] = sub.y.values; df = df.dropna(); df["x"] = (df.x - df.x.mean()) / df.x.std()
    f = sm.Logit(df.y, sm.add_constant(df.drop(columns="y"))).fit(disp=0); ci = f.conf_int().loc["x"]
    rows.append({"metabolite": t, "shared_bonf_both_outcomes": t in shared_set, "n": len(df), "n_both": int(df.y.sum()), "or": np.exp(f.params["x"]), "or_lower95": np.exp(ci[0]), "or_upper95": np.exp(ci[1]), "p": f.pvalues["x"]})
r = pd.DataFrame(rows); r["p_bonf"] = bonferroni(r.p.values); r = r.sort_values("p"); r.to_csv(d / "dual_vs_single_logistic.csv", index=False)
print(m.group.value_counts().to_dict(), "shared sig metabolites:", len(shared))
