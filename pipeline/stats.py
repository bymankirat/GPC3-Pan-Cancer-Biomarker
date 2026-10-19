"""Statistics written out by hand so they can be tested: AUC, bootstrap, BH, a fast log-rank test, cutpoint scan."""
import numpy as np


def auc(pos, neg):
    """P(random positive > random negative), ties count half (equals Mann-Whitney U / n1 n2)."""
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    allv = np.concatenate([pos, neg])
    order = allv.argsort(kind="mergesort")
    ranks = np.empty(len(allv))
    sv = allv[order]
    i = 0
    while i < len(sv):  # average ranks for ties
        j = i
        while j + 1 < len(sv) and sv[j + 1] == sv[i]:
            j += 1
        ranks[order[i:j + 1]] = (i + j) / 2 + 1
        i = j + 1
    r_pos = ranks[:len(pos)].sum()
    u = r_pos - len(pos) * (len(pos) + 1) / 2
    return u / (len(pos) * len(neg))


def bootstrap_auc(pos, neg, n, rng):
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    vals = [auc(rng.choice(pos, len(pos)), rng.choice(neg, len(neg))) for _ in range(n)]
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def bh(p):
    p = np.asarray(p, float)
    out = np.full(len(p), np.nan)
    ok = ~np.isnan(p)
    pv = p[ok]
    order = np.argsort(pv)
    q = pv[order] * len(pv) / (np.arange(len(pv)) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    res = np.empty(len(pv))
    res[order] = np.minimum(q, 1)
    out[ok] = res
    return out


def logrank(time, event, group):
    """Two-group log-rank test. Returns (chi-square, p). group is boolean."""
    from scipy.stats import chi2
    time, event, group = np.asarray(time, float), np.asarray(event, bool), np.asarray(group, bool)
    ut = np.unique(time[event])
    if len(ut) == 0 or group.all() or (~group).all():
        return 0.0, 1.0
    # at-risk counts at each event time, via sorting
    ts = np.sort(time)
    t1 = np.sort(time[group])
    n = len(ts) - np.searchsorted(ts, ut, side="left")
    n1 = len(t1) - np.searchsorted(t1, ut, side="left")
    d = _counts(time[event], ut)
    d1 = _counts(time[event & group], ut)
    e1 = d * n1 / n
    v = np.where(n > 1, d * (n1 / n) * (1 - n1 / n) * (n - d) / np.maximum(n - 1, 1), 0)
    var = v.sum()
    if var <= 0:
        return 0.0, 1.0
    stat = (d1.sum() - e1.sum()) ** 2 / var
    return float(stat), float(chi2.sf(stat, 1))


def _counts(values, grid):
    idx = np.searchsorted(grid, values)
    return np.bincount(idx, minlength=len(grid))[:len(grid)].astype(float)


def cutpoint_scan(x, time, event, qrange=(0.2, 0.8), step=0.01):
    """Log-rank p at every cutpoint between two quantiles. Returns (quantiles, p-values)."""
    qs = np.round(np.arange(qrange[0], qrange[1] + 1e-9, step), 3)
    cuts = np.quantile(x, qs)
    ps = np.array([logrank(time, event, x > c)[1] for c in cuts])
    return qs, ps


def cutpoint_permutation(x, time, event, n_perm, rng, qrange=(0.2, 0.8)):
    """How often does a random relabelling give a minimum p as small as the real one?"""
    _, ps = cutpoint_scan(x, time, event, qrange)
    obs = ps.min()
    null = np.empty(n_perm)
    for i in range(n_perm):
        _, pp = cutpoint_scan(rng.permutation(x), time, event, qrange)
        null[i] = pp.min()
    return float(obs), float((1 + (null <= obs).sum()) / (n_perm + 1)), null


def km(time, event):
    """Kaplan-Meier curve as step points (t, S)."""
    time, event = np.asarray(time, float), np.asarray(event, bool)
    ts = np.unique(time[event])
    s, out = 1.0, [(0.0, 1.0)]
    for t in ts:
        n = (time >= t).sum()
        d = ((time == t) & event).sum()
        s *= 1 - d / n
        out.append((float(t), float(s)))
    return out
