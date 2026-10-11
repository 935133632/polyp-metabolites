"""Specificity analysis: Model 1 Cox scan of all 251 metabolites for other digestive outcomes.
usage: python 13_specificity.py <code>      (code = 3-character ICD-10 K code, or CRC = C18-C20 combined)
Outcomes with >= 2000 incident events are analysed; first-diagnosis dates come from processed/digestive_first_dates.pkl (built by 13b)."""
import sys; sys.path.insert(0, r"D:\research\Xingang Li\polyp_metabolites\code")
from common import *

code = sys.argv[1]
df = pd.read_pickle(PROC / "digestive_first_dates.pkl")
sel = df[df.code.isin(["C18", "C19", "C20"])] if code == "CRC" else df[df.code == code]
first = sel.groupby("eid").date.min()
c = cohort_from_first_dates(first)
d = RES / "step13_specificity"; d.mkdir(exist_ok=True)
r = cox_scan(c); r["outcome"] = code; r.to_csv(d / f"{code}.csv", index=False)
print(code, "N", len(c), "events", int(c.event.sum()), "Bonferroni-significant", int((r.p_bonf < .05).sum()))
