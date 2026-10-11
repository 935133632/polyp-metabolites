"""Shared paths, loaders and statistics for the polyp x NMR-metabolite project."""
from __future__ import annotations
import numpy as np, pandas as pd
from pathlib import Path
from scipy.stats import norm, rankdata

WORK = Path(r"D:\research\Xingang Li\polyp_metabolites")
DATA = Path(r"D:\research\UKB_olink\Data")
PROC, RES, FIG, TAB = (WORK / d for d in ("processed", "results", "figures", "tables"))

# Outcome definitions. "polyp" = ICD-10 polyp-specific code; "any3" = original 3-character code of the Olink project.
OUTCOMES = {
    "K63": {"label": "Colon polyp (K63.5)", "prefix": "K63.5"},
    "K31": {"label": "Stomach/duodenal polyp (K31.7)", "prefix": "K31.7"},
    "K63any": {"label": "K63 any (original 3-char definition)", "prefix": "K63"},
    "K31any": {"label": "K31 any (original 3-char definition)", "prefix": "K31"},
}
NUM_COVARS = ["AgeRecr", "Sex", "TDI", "BMI", "SBPmean", "DBPmean"]
CAT_COVARS = ["SmokeStatus", "DrinkStatus"]


def design_covars(df: pd.DataFrame, drop=()) -> pd.DataFrame:
    """Numeric covariates + dummy-coded smoking/drinking (reference = first level alphabetically)."""
    cols = [c for c in NUM_COVARS if c not in drop]
    X = df[cols].copy()
    for c in CAT_COVARS:
        d = pd.get_dummies(df[c], prefix=c, dtype=float)
        d = d.iloc[:, 1:]                      # drop reference
        d[df[c].isna()] = np.nan
        X = pd.concat([X, d], axis=1)
    return X
MIN_COMPLETE, MIN_EVENTS = 100, 5


def bonferroni(p):
    p = np.asarray(p, float); return np.minimum(p * np.isfinite(p).sum(), 1.0)


def rint(x: pd.Series) -> pd.Series:
    """Rank-based inverse-normal transform (ties averaged); NaN preserved."""
    ok = x.notna()
    r = rankdata(x[ok], method="average")
    out = pd.Series(np.nan, index=x.index)
    out[ok] = norm.ppf((r - 0.5) / ok.sum())
    return out


def load_metabolome(transformed: bool = True) -> pd.DataFrame:
    f = PROC / ("metabolome_rint.pkl" if transformed else "metabolome_raw.pkl")
    return pd.read_pickle(f)


def load_base() -> pd.DataFrame:
    return pd.read_pickle(PROC / "base.pkl")


def load_outcome_dates() -> pd.DataFrame:
    return pd.read_pickle(PROC / "outcome_dates.pkl")


def build_cohort(outcome: str, exclude_first_years: float = 0.0) -> pd.DataFrame:
    """Incident cohort: prevalent cases removed; time = days from baseline to event/censoring."""
    base = load_base(); od = load_outcome_dates()
    c = base.merge(od[["eid", outcome]].rename(columns={outcome: "disease_first_date"}), on="eid", how="left")
    c["prevalent"] = c.disease_first_date.notna() & (c.disease_first_date <= c.baseline_date)
    c["incident"] = c.disease_first_date.notna() & (c.disease_first_date > c.baseline_date)
    c = c.loc[~c.prevalent].copy()
    # incident events after censoring date (death) are not counted
    c["incident"] &= c.disease_first_date <= c.censor_date
    c["event_date"] = np.where(c.incident, c.disease_first_date, c.censor_date)
    c["event_date"] = pd.to_datetime(c.event_date)
    c["follow_up_days"] = (c.event_date - c.baseline_date).dt.days.astype(float)
    c = c[c.follow_up_days > 0].copy()
    c["event"] = c.incident.astype(int)
    c["follow_up_years"] = c.follow_up_days / 365.25
    if exclude_first_years > 0:
        c = c[~((c.event == 1) & (c.follow_up_years < exclude_first_years))].copy()
    return c.reset_index(drop=True)


def prevalent_cohort(outcome: str) -> pd.DataFrame:
    base = load_base(); od = load_outcome_dates()
    c = base.merge(od[["eid", outcome]].rename(columns={outcome: "d"}), on="eid", how="left")
    c["prevalent"] = (c.d.notna() & (c.d <= c.baseline_date)).astype(int)
    return c


def _fit_one(name, x, T, E, X):
    from lifelines import CoxPHFitter
    df = X.copy(); df["x"] = x; df["T"] = T; df["E"] = E
    df = df.replace([np.inf, -np.inf], np.nan).dropna()
    n, ne = len(df), int(df.E.sum())
    r = {"metabolite": name, "n_complete": n, "n_events": ne}
    if n < MIN_COMPLETE or ne < MIN_EVENTS or df.x.std() == 0:
        return {**r, "status": "skipped"}
    df["x"] = (df.x - df.x.mean()) / df.x.std()
    df = df.loc[:, (df.nunique() > 1)]
    try:
        s = CoxPHFitter().fit(df, "T", "E").summary.loc["x"]
    except Exception as e:
        return {**r, "status": f"fail:{type(e).__name__}"}
    return {**r, "coef": s["coef"], "hr": s["exp(coef)"], "hr_lower95": s["exp(coef) lower 95%"],
            "hr_upper95": s["exp(coef) upper 95%"], "se": s["se(coef)"], "z": s["z"], "p": s["p"], "status": "ok"}


def cox_scan(cohort, metabolites=None, drop_covars=(), n_jobs=18):
    """Per-metabolite adjusted Cox (metabolite already RINT; re-standardised within complete cases)."""
    from joblib import Parallel, delayed
    M = load_metabolome().set_index("eid")
    names = list(metabolites) if metabolites is not None else list(M.columns)
    c = cohort.reset_index(drop=True)
    X = design_covars(c, drop_covars)
    Mc = M.reindex(c.eid.values)
    T, E = c.follow_up_days.values, c.event.values
    res = [_fit_one(m, Mc[m].values, T, E, X) for m in names]
    out = pd.DataFrame(res)
    out["p_bonf"] = np.minimum(out["p"] * out["p"].notna().sum(), 1.0)   # Bonferroni across fitted metabolites
    return out.sort_values("p").reset_index(drop=True)



def roc_points(y, score, n=200):
    """ROC curve thinned to at most n points (aggregate, not individual-level)."""
    from sklearn.metrics import roc_curve
    f, t, _ = roc_curve(y, score); idx = np.unique(np.linspace(0, len(f) - 1, n).astype(int))
    return pd.DataFrame({"fpr": f[idx], "tpr": t[idx]})
