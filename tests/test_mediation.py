"""Tests for the v0.3 Baron-Kenny mediation analysis."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("lifelines")
pytest.importorskip("sklearn")

from tp53_hrd.mediation import run_mediation


def _build_mediated_cohort(n: int = 200, seed: int = 42) -> pd.DataFrame:
    """Cohort where TP53 -> HRD -> survival with FULL mediation by design.

    Generative model:
      severity_score = U[0, 1]
      hrd_score      = 50 * severity_score + N(0, 5)   # path a (strong)
      log_hazard     = 0.05 * hrd_score                 # path b only
      => TP53's effect on survival is ENTIRELY via HRD
    """
    rng = np.random.default_rng(seed)
    sev = rng.uniform(0, 1, n)
    hrd = 50.0 * sev + rng.normal(0, 5, n)
    log_h = 0.05 * hrd
    t = rng.exponential(scale=np.exp(-log_h) * 365.0)
    e = (t < 600).astype(int)
    t = np.clip(t, 1, 1500)
    return pd.DataFrame({
        "severity_score": sev,
        "hrd_score": hrd,
        "os_days": t,
        "os_event": e,
    })


def _build_unmediated_cohort(n: int = 200, seed: int = 42) -> pd.DataFrame:
    """Cohort where TP53 -> survival DIRECTLY, with HRD independent (no mediation)."""
    rng = np.random.default_rng(seed)
    sev = rng.uniform(0, 1, n)
    hrd = rng.normal(25, 5, n)  # independent of sev
    log_h = 1.5 * sev            # TP53 effect, no HRD contribution
    t = rng.exponential(scale=np.exp(-log_h) * 365.0)
    e = (t < 600).astype(int)
    t = np.clip(t, 1, 1500)
    return pd.DataFrame({
        "severity_score": sev,
        "hrd_score": hrd,
        "os_days": t,
        "os_event": e,
    })


def test_mediation_returns_all_paths():
    df = _build_mediated_cohort(n=100)
    result = run_mediation(df, n_bootstrap=200, seed=42)
    assert result.n == 100
    # Path a (TP53 -> HRD) should be strongly positive (we wired 50x)
    assert result.a_coef > 30, f"expected strong path-a (~50), got {result.a_coef}"
    # Path b (HRD -> survival | TP53) should be positive in log-HR space
    assert result.b_log_hr > 0


def test_mediation_indirect_effect_present_in_mediated_cohort():
    df = _build_mediated_cohort(n=300, seed=42)
    result = run_mediation(df, n_bootstrap=500, seed=42)
    # Indirect effect = a * b should be substantially positive
    assert result.indirect_log_hr > 0.1, (
        f"expected positive indirect log-HR, got {result.indirect_log_hr}"
    )
    # Bootstrap CI should NOT cross zero
    assert result.indirect_ci_low is not None
    assert result.indirect_ci_low > 0, (
        f"95% CI on indirect effect crosses zero: [{result.indirect_ci_low}, {result.indirect_ci_high}]"
    )


def test_mediation_proportion_mediated_close_to_one_when_fully_mediated():
    df = _build_mediated_cohort(n=400, seed=42)
    result = run_mediation(df, n_bootstrap=300, seed=42)
    assert result.proportion_mediated is not None
    # With full mediation, indirect should be >= 50% of total
    assert result.proportion_mediated > 0.5, (
        f"expected substantial mediation (>0.5), got {result.proportion_mediated}"
    )


def test_mediation_proportion_small_when_no_mediation():
    df = _build_unmediated_cohort(n=400, seed=42)
    result = run_mediation(df, n_bootstrap=300, seed=42)
    # When HRD is independent of TP53, path a should be ~0, so indirect ~0
    assert abs(result.a_coef) < 5, f"expected near-zero path-a, got {result.a_coef}"
    assert abs(result.indirect_log_hr) < 0.2


def test_mediation_raises_on_missing_columns():
    df = pd.DataFrame({"severity_score": [0.5], "os_days": [100], "os_event": [1]})
    with pytest.raises(KeyError):
        run_mediation(df, n_bootstrap=0)
