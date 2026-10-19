"""Step 3: does GPC3 predict survival? Done the way a clinical statistician would.

Phase 1 split LIHC patients at the median and found no difference. Here:
1. Kaplan-Meier with a median split, for overall, disease-specific and progression-free survival
   (TCGA Clinical Data Resource endpoints, recommended for TCGA by Liu et al. 2018).
2. Cox regression with GPC3 as a continuous variable (per standard deviation), alone and
   adjusted for age, sex and stage, with a check of the proportional-hazards assumption.
3. The "optimal cutpoint" trap: many web tools pick the cutpoint with the smallest p-value.
   A permutation test shows how small that minimum p gets by chance alone.
4. The same continuous Cox model in every TCGA cancer, as a forest plot.
"""
import numpy as np
import pandas as pd

from stats import bh, cutpoint_permutation, cutpoint_scan, km, logrank
from util import DERIVED, RESULTS, config, log, save_json

ENDPOINTS = {"OS": "overall survival", "DSS": "disease-specific survival", "PFI": "progression-free interval"}


def stage_group(s):
    s = str(s).upper()
    if "IV" in s or "III" in s:
        return 1.0
    if "STAGE I" in s:
        return 0.0
    return np.nan


def data():
    pan = pd.read_csv(DERIVED / "pancancer_expression.csv.gz")
    surv = pd.read_csv(DERIVED / "survival.csv.gz")
    tum = pan[(pan.study == "TCGA") & (pan.sample_type.isin(["Primary Tumor", "Primary Blood Derived Cancer - Peripheral Blood"]))]
    tum = tum.assign(patient=tum["sample"].str[:12]).drop_duplicates("patient")
    df = tum[["patient", config.GENE]].merge(surv.drop_duplicates("_PATIENT"), left_on="patient", right_on="_PATIENT")
    df["age"] = pd.to_numeric(df["age_at_initial_pathologic_diagnosis"], errors="coerce")
    df["male"] = (df["gender"].astype(str).str.upper() == "MALE").astype(float)
    df["late_stage"] = df["ajcc_pathologic_tumor_stage"].map(stage_group)
    return df


def cox(d, cols, time, event):
    from lifelines import CoxPHFitter
    d = d[cols + [time, event]].dropna()
    d = d[d[time] > 0]
    if len(d) < 30 or d[event].sum() < 10:
        return None
    cph = CoxPHFitter().fit(d, duration_col=time, event_col=event)
    s = cph.summary.loc["z"]
    out = {"n": len(d), "events": int(d[event].sum()), "hr": float(s["exp(coef)"]), "lo": float(s["exp(coef) lower 95%"]),
           "hi": float(s["exp(coef) upper 95%"]), "p": float(s["p"])}
    try:
        from lifelines.statistics import proportional_hazard_test
        out["ph_p"] = float(proportional_hazard_test(cph, d, time_transform="rank").summary.loc["z", "p"])
    except Exception:
        out["ph_p"] = None
    return out


def main():
    rng = np.random.default_rng(config.SEED)
    df = data()
    g = config.GENE
    f = df[df["cancer type abbreviation"] == config.FOCUS_CANCER].copy()
    f["z"] = (f[g] - f[g].mean()) / f[g].std()
    res = {"n_patients": len(f), "endpoints": {}}
    for ep, label in ENDPOINTS.items():
        t, e = f"{ep}.time", ep
        d = f.dropna(subset=[t, e])
        d = d[d[t] > 0]
        high = d[g] > d[g].median()
        chi, p = logrank(d[t], d[e], high)
        uni = cox(d.assign(z=d.z), ["z"], t, e)
        adj = cox(d, ["z", "age", "male", "late_stage"], t, e)
        qs, ps = cutpoint_scan(d[g].values, d[t].values, d[e].values.astype(bool), config.CUTPOINT_RANGE)
        obs, p_adj, null = cutpoint_permutation(d[g].values, d[t].values, d[e].values.astype(bool), config.N_PERM, rng, config.CUTPOINT_RANGE)
        res["endpoints"][ep] = {
            "label": label, "n": int(len(d)), "events": int(d[e].sum()), "median_split_p": p,
            "km_high": km(d.loc[high, t] / 365.25, d.loc[high, e]), "km_low": km(d.loc[~high, t] / 365.25, d.loc[~high, e]),
            "n_high": int(high.sum()), "n_low": int((~high).sum()),
            "cox_univariable": uni, "cox_adjusted": adj,
            "cutpoint": {"quantiles": qs.tolist(), "p": ps.tolist(), "min_p": obs, "best_quantile": float(qs[ps.argmin()]),
                         "permutation_p": p_adj, "null_fraction_below_005": float((null < 0.05).mean())}}
        log(f"{ep}: median split p={p:.3f}; Cox per SD HR={uni['hr']:.2f} p={uni['p']:.3f}; adjusted HR={adj['hr'] if adj else float('nan'):.2f}; "
            f"best cutpoint p={obs:.4f} -> permutation p={p_adj:.3f}")

    # every cancer, overall survival, continuous GPC3
    rows = []
    for c, d in df.groupby("cancer type abbreviation"):
        if d[g].std() == 0 or d[g].isna().all():
            continue
        d = d.assign(z=(d[g] - d[g].mean()) / d[g].std())
        r = cox(d, ["z"], "OS.time", "OS")
        if r:
            rows.append({"cancer": c, **r})
    pc = pd.DataFrame(rows)
    pc["q"] = bh(pc.p.values)
    pc.sort_values("hr").to_csv(RESULTS / "survival_pancancer.csv", index=False)
    res["pancancer_sig"] = pc[pc.q < config.FDR].cancer.tolist()
    save_json(RESULTS / "survival.json", res)
    log(f"pan-cancer OS: significant (q<0.05) in {res['pancancer_sig']}")


if __name__ == "__main__":
    main()
