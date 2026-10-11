import sys, time; sys.path.insert(0, r"D:\research\Xingang Li\polyp_metabolites\code")
from common import *
out = RES / "step01_primary"; out.mkdir(exist_ok=True)
for k in OUTCOMES:
    t = time.time(); c = build_cohort(k)
    r = cox_scan(c); r.to_csv(out / f"model1_{k}_cox_all_metabolites.csv", index=False)
    ok = r[r.status == "ok"]
    print(k, "N", len(c), "events", c.event.sum(), "fit", len(ok), "Bonferroni<0.05", (ok.p_bonf < .05).sum(), f"{time.time()-t:.0f}s")
    print(ok.head(5)[["metabolite", "hr", "hr_lower95", "hr_upper95", "p", "p_bonf"]].to_string(index=False))
