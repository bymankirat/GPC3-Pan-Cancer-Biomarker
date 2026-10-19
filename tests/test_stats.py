"""Unit tests of the statistics on inputs with known answers."""
import numpy as np
import pytest

from stats import auc, bh, bootstrap_auc, cutpoint_permutation, cutpoint_scan, km, logrank


def test_auc_matches_pair_counting():
    rng = np.random.default_rng(0)
    a, b = rng.integers(0, 5, 40).astype(float), rng.integers(0, 5, 30).astype(float)
    pairs = np.mean([(x > y) + 0.5 * (x == y) for x in a for y in b])
    assert auc(a, b) == pytest.approx(pairs)
    assert auc([3, 4], [1, 2]) == 1.0 and auc([1, 2], [3, 4]) == 0.0


def test_bootstrap_ci_contains_estimate_and_is_reproducible():
    rng = np.random.default_rng(1)
    a, b = rng.normal(1, 1, 80), rng.normal(0, 1, 60)
    lo, hi = bootstrap_auc(a, b, 300, np.random.default_rng(5))
    assert lo < auc(a, b) < hi
    assert (lo, hi) == bootstrap_auc(a, b, 300, np.random.default_rng(5))


def test_bh_known_values():
    q = bh([0.01, 0.04, 0.03, np.nan])
    assert q[:3].tolist() == pytest.approx([0.03, 0.04, 0.04]) and np.isnan(q[3])


def test_logrank_matches_lifelines():
    lifelines = pytest.importorskip("lifelines")
    from lifelines.statistics import logrank_test
    rng = np.random.default_rng(2)
    t = np.round(rng.exponential(500, 200))
    e = rng.random(200) < 0.6
    g = rng.random(200) < 0.5
    t[g] *= 0.7
    ref = logrank_test(t[g], t[~g], e[g], e[~g])
    stat, p = logrank(t, e, g)
    assert stat == pytest.approx(ref.test_statistic, rel=1e-6) and p == pytest.approx(ref.p_value, rel=1e-6)


def test_minimum_p_cutpoint_is_optimistic_and_permutation_corrects_it():
    rng = np.random.default_rng(3)
    n = 200
    x, t, e = rng.normal(size=n), rng.exponential(1000, n), rng.random(n) < 0.5   # no real effect
    qs, ps = cutpoint_scan(x, t, e)
    assert ps.min() <= ps[np.argmin(np.abs(qs - 0.5))]
    obs, p_adj, null = cutpoint_permutation(x, t, e, 100, rng)
    assert p_adj >= obs and p_adj > 0.05
    assert (null < 0.05).mean() > 0.1          # chance alone often gives a "significant" best cutpoint


def test_km_simple_case():
    curve = km([1, 2, 2, 3], [1, 1, 0, 1])
    assert curve == [(0.0, 1.0), (1.0, 0.75), (2.0, 0.5), (3.0, 0.0)]
