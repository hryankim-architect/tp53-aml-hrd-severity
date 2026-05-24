"""Tests for the v0.3 nested Cox model comparison module."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("lifelines")

from tp53_hrd.model_compare import compare_nested, run_nested_model_suite


def _synthetic_two_axis_cohort(n: int = 80, seed: int = 42) -> pd.DataFrame:
    """Cohort where BOTH axes truly drive hazard (so M2 should beat M1)."""
    rng = np.random.default_rng(seed)
    sev = rng.uniform(0, 1, n)
    hrd = rng.integers(0, 60, n).astype(float)
    log_h = 1.2 * sev + 0.04 * hrd
    t = rng.exponential(scale=np.exp(-log_h) * 365.0)
    e = (t < 600).astype(int)
    t = np.clip(t, 1, 1500)
    return pd.DataFrame({
        "severity_score": sev,
        "hrd_score": hrd,
        "os_days": t,
        "os_event": e,
    })


def _synthetic_single_axis_cohort(n: int = 80, seed: int = 42) -> pd.DataFrame:
    """Cohort where ONLY TP53 drives hazard; HRD is noise (so M2 should NOT beat M1)."""
    rng = np.random.default_rng(seed)
    sev = rng.uniform(0, 1, n)
    hrd = rng.integers(0, 60, n).astype(float)  # uncorrelated noise
    log_h = 1.5 * sev  # HRD coefficient = 0
    t = rng.exponential(scale=np.exp(-log_h) * 365.0)
    e = (t < 600).astype(int)
    t = np.clip(t, 1, 1500)
    return pd.DataFrame({
        "severity_score": sev,
        "hrd_score": hrd,
        "os_days": t,
        "os_event": e,
    })


def test_compare_nested_requires_complex_has_more_params():
    from lifelines import CoxPHFitter
    df = _synthetic_two_axis_cohort()
    cph_a = CoxPHFitter()
    cph_a.fit(df[["severity_score", "os_days", "os_event"]], "os_days", "os_event")
    cph_b = CoxPHFitter()
    cph_b.fit(df[["severity_score", "os_days", "os_event"]], "os_days", "os_event")
    with pytest.raises(ValueError):
        compare_nested(cph_a, cph_b, "same", "same")


def test_run_nested_model_suite_returns_three_pairs():
    df = _synthetic_two_axis_cohort()
    result = run_nested_model_suite(df)
    assert set(result.keys()) == {"m2_vs_m1", "m3_vs_m2", "m3_vs_m1"}
    for _k, v in result.items():
        assert "lrt_p" in v
        assert "aic_delta" in v
        assert "cindex_delta" in v
        assert "justifies_complex" in v


def test_run_nested_model_suite_prefers_m2_when_both_axes_drive():
    """When both axes truly contribute, M2 should beat M1 (LRT p < 0.05, AIC delta < 0)."""
    df = _synthetic_two_axis_cohort(n=200, seed=42)
    result = run_nested_model_suite(df)
    m2m1 = result["m2_vs_m1"]
    assert m2m1["lrt_p"] < 0.05, f"expected significant LRT, got p={m2m1['lrt_p']}"
    assert m2m1["aic_delta"] < 0, f"expected AIC improvement, got delta={m2m1['aic_delta']}"
    assert m2m1["justifies_complex"] is True


def test_run_nested_model_suite_does_not_prefer_m2_when_hrd_is_noise():
    """When HRD is pure noise, M2 should NOT meaningfully beat M1."""
    df = _synthetic_single_axis_cohort(n=200, seed=42)
    result = run_nested_model_suite(df)
    m2m1 = result["m2_vs_m1"]
    # With pure-noise HRD, LRT p should be high or AIC delta should be small/positive
    assert not m2m1["justifies_complex"], (
        f"complex model wrongly preferred: lrt_p={m2m1['lrt_p']}, aic_delta={m2m1['aic_delta']}"
    )


def test_run_nested_model_suite_missing_columns():
    df = pd.DataFrame({"severity_score": [0.1], "os_days": [10], "os_event": [1]})
    with pytest.raises(KeyError):
        run_nested_model_suite(df)
