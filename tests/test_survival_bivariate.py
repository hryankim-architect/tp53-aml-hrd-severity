"""Smoke tests for `cox_bivariate` (v0.2 multi-covariate Cox)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("lifelines")

from tp53_hrd.survival import cox_bivariate


def _synthetic_cohort(n: int = 60, seed: int = 42) -> pd.DataFrame:
    """Two-axis cohort: severity_score + hrd_score independently predict event."""
    rng = np.random.default_rng(seed)
    sev = rng.uniform(0, 1, n)
    hrd = rng.integers(0, 60, n)
    # Linear log-hazard from both axes -> exponential time
    log_h = 0.8 * sev + 0.02 * hrd
    t = rng.exponential(scale=np.exp(-log_h) * 365.0)
    e = (t < 600).astype(int)
    t = np.clip(t, 1, 1500)
    return pd.DataFrame({
        "severity_score": sev,
        "hrd_score": hrd.astype(float),
        "os_days": t,
        "os_event": e,
    })


def test_cox_bivariate_returns_per_covariate_hr():
    df = _synthetic_cohort()
    result = cox_bivariate(df, ["severity_score", "hrd_score"])
    assert set(result["model_summary"].keys()) == {"severity_score", "hrd_score"}
    for cov in ("severity_score", "hrd_score"):
        m = result["model_summary"][cov]
        assert m["hr"] > 0
        assert m["ci_low"] <= m["hr"] <= m["ci_high"]
        assert 0 <= m["p_value"] <= 1
    assert result["n"] == 60
    assert result["interaction_term"] is None


def test_cox_bivariate_with_interaction_adds_interaction_term():
    df = _synthetic_cohort()
    result = cox_bivariate(
        df, ["severity_score", "hrd_score"], interaction=True
    )
    assert result["interaction_term"] == "severity_score:hrd_score"
    assert "severity_score:hrd_score" in result["model_summary"]


def test_cox_bivariate_raises_on_missing_columns():
    df = pd.DataFrame({"severity_score": [0.5], "os_days": [100], "os_event": [1]})
    with pytest.raises(KeyError):
        cox_bivariate(df, ["severity_score", "hrd_score"])
