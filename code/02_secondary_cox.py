"""Secondary Cox analyses on all 251 metabolites.
usage:  python 02_secondary_cox.py main K63          -> sensitivity (exclude events <2y), subgroups, split-sample validation
        python 02_secondary_cox.py age K63           -> age strata only (60 main, 65 supplementary)
        python 02_secondary_cox.py boot K63 0 5      -> bootstrap replicates 0..4 (run several chunks in parallel)
        python 02_secondary_cox.py merge K63         -> merge bootstrap chunks into summary
"""
import sys; sys.path.insert(0, r"D:\research\Xingang Li\polyp_metabolites\code")
from common import *

mode, k = sys.argv[1], sys.argv[2]
NBOOT = 30


def seed(*a):
    import hashlib
    return np.random.default_rng(int(hashlib.md5("|".join(map(str, a)).encode()).hexdigest()[:8], 16))


def run_subgroups(c, k, ages=True, sexes=True):
    """Sex strata and age strata at 60 (main) and 65 (supplementary), all 251 metabolites per stratum."""
    d = RES / "step03_subgroup"; d.mkdir(exist_ok=True)
    f_main, f_65 = d / f"{k}_subgroups.csv", d / f"{k}_subgroups_age65.csv"
    spec = []
    if sexes: spec += [("Female", c.Sex == 0, ("Sex",), "main"), ("Male", c.Sex == 1, ("Sex",), "main")]
    if ages: spec += [("Age<60", c.AgeRecr < 60, ("AgeRecr",), "main"), ("Age>=60", c.AgeRecr >= 60, ("AgeRecr",), "main"),
                      ("Age<65", c.AgeRecr < 65, ("AgeRecr",), "65"), ("Age>=65", c.AgeRecr >= 65, ("AgeRecr",), "65")]
    old = pd.read_csv(f_main) if (f_main.exists() and not (ages and sexes)) else pd.DataFrame()
    if len(old): old = old[~old.subgroup.str.startswith("Age")] if ages else old[old.subgroup.str.startswith("Age")]   # keep the strata not being re-run
    out_main, out_65 = [old] if len(old) else [], []
    for name, mask, drop, where in spec:
        r = cox_scan(c[mask], None, drop); r["subgroup"] = name
        (out_main if where == "main" else out_65).append(r)
        if out_main: pd.concat(out_main).to_csv(f_main, index=False)
        if out_65: pd.concat(out_65).to_csv(f_65, index=False)


if mode == "main":
    d = RES / "step02_sensitivity_exclude_first2y"; d.mkdir(exist_ok=True)
    if not (d / f"{k}_sensitivity.csv").exists():
        cox_scan(build_cohort(k, exclude_first_years=2)).to_csv(d / f"{k}_sensitivity.csv", index=False)
    c = build_cohort(k)
    # subgroups: all metabolites; sex or age removed from covariates inside its own stratum
    run_subgroups(c, k, ages=True, sexes=True)
    # split-sample: 70/30 stratified by event, all metabolites in both halves
    d = RES / "step06_internal_validation"; d.mkdir(exist_ok=True)
    rng = seed("split", k)
    ev = c.index[c.event == 1].to_numpy(); ne = c.index[c.event == 0].to_numpy()
    disc = np.concatenate([rng.choice(ev, int(.7 * len(ev)), replace=False), rng.choice(ne, int(.7 * len(ne)), replace=False)])
    mask = c.index.isin(disc)
    a = cox_scan(c[mask]); b = cox_scan(c[~mask])
    m = a.merge(b, on="metabolite", suffixes=("_disc", "_valid"))
    m["same_direction"] = np.sign(m.coef_disc) == np.sign(m.coef_valid)
    m["discovery_bonf05"] = m.p_bonf_disc < .05
    m["replicated_p05"] = m.same_direction & (m.p_valid < .05)
    m.to_csv(d / f"{k}_split_sample.csv", index=False)

elif mode == "age":
    run_subgroups(build_cohort(k), k, ages=True, sexes=False)

elif mode == "boot":
    d = RES / "step06_internal_validation" / "bootstrap_chunks"; d.mkdir(parents=True, exist_ok=True)
    c = build_cohort(k)
    for b in range(int(sys.argv[3]), int(sys.argv[4])):
        f = d / f"{k}_rep{b:02d}.csv"
        if f.exists(): continue
        rng = seed("boot", k, b)
        cb = c.iloc[rng.integers(0, len(c), len(c))].reset_index(drop=True)
        r = cox_scan(cb); r["rep"] = b; r.to_csv(f, index=False)

elif mode == "merge":
    d = RES / "step06_internal_validation"
    bs = pd.concat([pd.read_csv(f) for f in sorted((d / "bootstrap_chunks").glob(f"{k}_rep*.csv"))])
    bs = bs[bs.status == "ok"]
    g = bs.groupby("metabolite").agg(n_reps=("rep", "nunique"), freq_p05=("p", lambda x: (x < .05).mean()), freq_bonf05=("p_bonf", lambda x: (x < .05).mean()),
                                     prop_positive=("coef", lambda x: (x > 0).mean()), median_coef=("coef", "median")).reset_index()
    g.to_csv(d / f"{k}_bootstrap_all_metabolites.csv", index=False)
    print(k, "reps merged:", bs.rep.nunique(), "metabolites:", len(g))
