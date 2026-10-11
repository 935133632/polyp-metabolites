"""Extra covariates for the extended adjustment models: waist circumference, body fat %, CRP (log), HbA1c (baseline instance)."""
import sys; sys.path.insert(0, r"D:\research\Xingang Li\polyp_metabolites\code")
from common import *

base = load_base()
an = pd.read_csv(DATA / "WH-PhyMeasure-Anthropometry-20251128.csv", usecols=["eid", "p48_i0", "p23099_i0"]).rename(columns={"p48_i0": "Waist", "p23099_i0": "BodyFatPct"})
bio = pd.read_csv(DATA / "FRS-Biochemistry-Inst0-20251127.csv", usecols=["eid", "p30710_i0", "p30750_i0"]).rename(columns={"p30710_i0": "CRP", "p30750_i0": "HbA1c"})
x = base[["eid"]].merge(an, on="eid", how="left").merge(bio, on="eid", how="left")
x["CRP"] = np.log(x.CRP.where(x.CRP > 0))        # skewed; natural log (mg/L)
x.to_pickle(PROC / "extra_covars.pkl")
print(x.drop(columns="eid").describe().round(2).T[["count", "mean", "std", "min", "max"]]); print("missing share:", x.drop(columns="eid").isna().mean().round(3).to_dict())
c = build_cohort("K63"); print("check K63 cohort:", len(c), int(c.event.sum()))
