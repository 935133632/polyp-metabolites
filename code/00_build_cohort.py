"""Build analysis base table, outcome dates and metabolite matrices from WH-eid UKB files (read-only inputs)."""
import sys; sys.path.insert(0, r"D:\research\Xingang Li\polyp_metabolites\code")
import numpy as np, pandas as pd, pyreadr
from common import *

# ---- metabolome
m = pyreadr.read_r(DATA / "09.UKB Metabolome Clean Data.Rdata")["df1"]
m = m.loc[m.drop(columns="eid").notna().any(axis=1)].reset_index(drop=True)
m["eid"] = m.eid.astype(int)
m.to_pickle(PROC / "metabolome_raw.pkl")
mt = m.copy()
for c in mt.columns[1:]:
    mt[c] = rint(mt[c])
mt.to_pickle(PROC / "metabolome_rint.pkl")
eids = set(m.eid)
print("metabolome participants", len(m), "metabolites", m.shape[1] - 1)

# ---- baseline / covariates
b = pd.read_csv(DATA / "WH-BasicChara-20251127.csv", usecols=["eid", "p21022", "p31", "p22189", "p53_i0"])
b = b[b.eid.isin(eids)].rename(columns={"p21022": "AgeRecr", "p22189": "TDI"})
b["Sex"] = b.p31.map({"Female": 0.0, "Male": 1.0}); b["baseline_date"] = pd.to_datetime(b.p53_i0)
b = b.drop(columns=["p31", "p53_i0"])
an = pd.read_csv(DATA / "WH-PhyMeasure-Anthropometry-20251128.csv", usecols=["eid", "p21001_i0"]).rename(columns={"p21001_i0": "BMI"})
bp = pd.read_csv(DATA / "WH-PhyMeasure-BP-20251128.csv", usecols=["eid", "p4080_i0_a0", "p4080_i0_a1", "p93_i0_a0", "p93_i0_a1", "p4079_i0_a0", "p4079_i0_a1", "p94_i0_a0", "p94_i0_a1"])
bp["SBPmean"] = bp[["p4080_i0_a0", "p4080_i0_a1"]].mean(axis=1).fillna(bp[["p93_i0_a0", "p93_i0_a1"]].mean(axis=1))
# automated (p4079) with manual (p94) fallback, mirroring SBP
bp["DBPmean"] = bp[["p4079_i0_a0", "p4079_i0_a1"]].mean(axis=1).fillna(bp[["p94_i0_a0", "p94_i0_a1"]].mean(axis=1))
base = b.merge(an, on="eid", how="left").merge(bp[["eid", "SBPmean", "DBPmean"]], on="eid", how="left")

# ---- smoking / drinking (baseline instance; "Prefer not to answer" -> missing)
sm = pd.read_csv(DATA / "WH-TouchScreen-LifeStyle-Smoke-20251128.csv", usecols=["eid", "p20116_i0"]).rename(columns={"p20116_i0": "SmokeStatus"})
al = pd.read_csv(DATA / "WH-TouchScreen-LifeStyle-Alcohol-20251128.csv", usecols=["eid", "p20117_i0"]).rename(columns={"p20117_i0": "DrinkStatus"})
for df_, c in ((sm, "SmokeStatus"), (al, "DrinkStatus")):
    df_[c] = df_[c].replace({"Prefer not to answer": np.nan})
base = base.merge(sm, on="eid", how="left").merge(al, on="eid", how="left")
print(base.SmokeStatus.value_counts(dropna=False).to_dict(), base.DrinkStatus.value_counts(dropna=False).to_dict())

# ---- death
d = pd.read_csv(DATA / "WH-Outcome-Death-20251128.csv", usecols=["eid", "p40000_i0", "p40000_i1"])
d["death_date"] = pd.to_datetime(d.p40000_i0.fillna(d.p40000_i1))
base = base.merge(d[["eid", "death_date"]], on="eid", how="left")

# ---- ICD-10 first dates
dcols_a = [f"p41280_a{i}" for i in range(198)]; dcols_b = [f"p41280_a{i}" for i in range(198, 259)]
PB = pd.read_csv(DATA / "WH-Outcome-Diagnosis-ICD10-PartB-20251128.csv", usecols=["eid", *dcols_b], low_memory=False).set_index("eid")
ra = pd.read_csv(DATA / "WH-Outcome-Diagnosis-ICD10-PartA-20251128.csv", usecols=["eid", "p41270", *dcols_a], chunksize=20000, low_memory=False)
rows, end_date = [], pd.Timestamp.min
for ca in ra:
    ca = ca[ca.eid.isin(eids)].reset_index(drop=True)
    cb = PB.reindex(ca.eid.values).reset_index(drop=True)
    dates = pd.concat([ca[dcols_a], cb[dcols_b]], axis=1).apply(pd.to_datetime, errors="coerce")
    end_date = max(end_date, dates.max().max())
    D = dates.to_numpy()
    for i, (eid, s) in enumerate(zip(ca.eid.values, ca.p41270.values)):
        if not isinstance(s, str): continue
        first = {}
        for j, code in enumerate(s.split("|")):
            code = code.strip(); dt = D[i, j] if j < D.shape[1] else np.datetime64("NaT")
            if pd.isna(dt): continue
            for k, spec in OUTCOMES.items():
                if code.startswith(spec["prefix"]):
                    if k not in first or dt < first[k]: first[k] = dt
        if first: rows.append({"eid": int(eid), **first})
od = pd.DataFrame(rows)
for k in OUTCOMES:
    if k not in od: od[k] = pd.NaT
od.to_pickle(PROC / "outcome_dates.pkl")
print("cohort end date", end_date)
base["censor_date"] = base.death_date.where(base.death_date.notna() & (base.death_date < end_date), end_date)
base = base.sort_values("eid").reset_index(drop=True)
base.to_pickle(PROC / "base.pkl")
print(base.isna().mean().round(4).to_dict())
for k in OUTCOMES:
    c = build_cohort(k); p = prevalent_cohort(k).prevalent.sum()
    print(k, "N", len(c), "incident events", c.event.sum(), "prevalent", p, "median FU", round(c.follow_up_years.median(), 2))
