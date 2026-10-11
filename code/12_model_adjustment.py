"""Extended adjustment of the primary Cox scan (all 251 metabolites), always in the same complete-case sample.
usage: python 12_model_adjustment.py <K63|K31> <M1|M2a|M2b|M2c>
  M1  = Model 1 (age, sex, TDI, BMI, SBP, DBP, smoking, drinking)
  M2a = M1 + waist circumference + body fat %
  M2b = M1 + HbA1c + log CRP
  M2c = M1 + waist + body fat % + HbA1c + log CRP"""
import sys; sys.path.insert(0, r"D:\research\Xingang Li\polyp_metabolites\code")
from common import *

k, model = sys.argv[1], sys.argv[2]
EXTRA = {"M1": [], "M2a": ["Waist", "BodyFatPct"], "M2b": ["HbA1c", "CRP"], "M2c": ["Waist", "BodyFatPct", "HbA1c", "CRP"]}
c = build_cohort(k).merge(pd.read_pickle(PROC / "extra_covars.pkl"), on="eid", how="left")
c = c.dropna(subset=["Waist", "BodyFatPct", "HbA1c", "CRP"]).reset_index(drop=True)      # identical sample for every model
d = RES / "step12_model_adjustment"; d.mkdir(exist_ok=True)
r = cox_scan(c, extra=EXTRA[model]); r["model"] = model; r.to_csv(d / f"{k}_{model}.csv", index=False)
print(k, model, "N", len(c), "events", int(c.event.sum()), "Bonferroni-significant", int((r.p_bonf < .05).sum()))
