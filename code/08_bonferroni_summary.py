import sys; sys.path.insert(0, r"D:\research\Xingang Li\polyp_metabolites\code")
from common import *
from metab_annot import metab_class
thr = .05 / 251; print("threshold", thr)
for k in ["K63", "K31"]:
    p = pd.read_csv(RES / f"step01_primary/model1_{k}_cox_all_metabolites.csv").query("status=='ok'"); bs = p[p.p < thr]
    print(f"\n== {k}: primary Bonferroni {len(bs)}")
    print(bs.assign(c=bs.metabolite.map(metab_class)).groupby("c").size().sort_values(ascending=False).to_string())
    g = pd.read_csv(RES / f"step03_subgroup/{k}_subgroups.csv").query("status=='ok'")
    print("subgroups Bonf per stratum:", {s: int((d.p < thr).sum()) for s, d in g.groupby("subgroup")}, "| main-Bonf sig in ALL 4 strata:",
          int((g[g.metabolite.isin(bs.metabolite)].assign(s=lambda d: d.p < thr).groupby("metabolite").s.all()).sum()), "| same direction in all 4:",
          int(g[g.metabolite.isin(bs.metabolite)].assign(m=lambda d: d.metabolite.map(bs.set_index("metabolite").coef)).assign(ok=lambda d: d.coef * d.m > 0).groupby("metabolite").ok.all().sum()))
    s = pd.read_csv(RES / f"step06_internal_validation/{k}_split_sample.csv"); d = s[s.p_disc < thr]
    print("split: discovery Bonf", len(d), "| same dir", int(d.same_direction.sum()), "| validation P<.05", int(d.replicated_p05.sum()), "| validation Bonf", int(((d.p_valid < thr) & d.same_direction).sum()))
    b = pd.read_csv(RES / f"step06_internal_validation/{k}_bootstrap_all_metabolites.csv"); bb = b[b.metabolite.isin(bs.metabolite)]
    print("bootstrap among main-Bonf: P<.05 in >=95% reps", int((bb.freq_p05 >= .95).sum()), "/", len(bb))
