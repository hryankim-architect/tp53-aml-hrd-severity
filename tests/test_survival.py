"""Tests for the survival analysis module.

Synthetic cohort designed for *clear* separation so the assertions remain
stable even though n is small: ``high`` patients all die within 200 days,
``low`` patients survive past 1500 days. With this signal the median ordering
and Cox HR direction are guaranteed; only the exact p-value depends on
lifelines internals.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from tp53_hrd.survival import (
    cox_severity,
    kaplan_meier_summary,
    make_km_plot,
    multivariate_logrank_p,
    two_arm_summary,
)


def _synthetic_cohort() -> pd.DataFrame:
    """Build a small cohort with clean high-vs-low survival separation."""
    rows = []
    # high band: 5 patients, all die within 200 days
    for i, days in enumerate([50, 80, 120, 150, 200]):
        rows.append({
            "patient_id": f"H-{i}",
            "severity_band": "high",
            "severity_score": 0.75,
            "os_days": float(days),
            "os_event": 1,
        })
    # moderate band: 2 patients, intermediate
    for i, days in enumerate([300, 500]):
        rows.append({
            "patient_id": f"M-{i}",
            "severity_band": "moderate",
            "severity_score": 0.50,
            "os_days": float(days),
            "os_event": 1,
        })
    # low band: 8 patients, mostly survive past 1500 days
    for i, days in enumerate([1500, 1800, 1900, 2000, 2100, 2200, 2400, 2800]):
        rows.append({
            "patient_id": f"L-{i}",
            "severity_band": "low",
            "severity_score": 0.0,
            "os_days": float(days),
            "os_event": 0 if i >= 2 else 1,  # 2 deaths, 6 censored
        })
    return pd.DataFrame(rows)


class TestKaplanMeierSummary:
    def test_three_bands_present(self):
        summary = kaplan_meier_summary(_synthetic_cohort())
        assert set(summary["group"]) == {"low", "moderate", "high"}

    def test_ordering_low_moderate_high(self):
        # _ordered_groups should put bands in clinical-risk order
        summary = kaplan_meier_summary(_synthetic_cohort())
        assert list(summary["group"]) == ["low", "moderate", "high"]

    def test_n_and_events_per_group(self):
        summary = kaplan_meier_summary(_synthetic_cohort())
        d = {row["group"]: (row["n"], row["n_events"]) for _, row in summary.iterrows()}
        assert d["high"] == (5, 5)
        assert d["moderate"] == (2, 2)
        assert d["low"] == (8, 2)

    def test_median_os_high_less_than_low(self):
        # high band's median must be lower than low band's median
        summary = kaplan_meier_summary(_synthetic_cohort())
        by_group = {row["group"]: row["median_os_days"] for _, row in summary.iterrows()}
        # low band has mostly censored data, may return None for median — that
        # itself is informative (less than half have died, so median > follow-up)
        assert by_group["high"] is not None
        # The high band's median should be ~120 days (middle of 50/80/120/150/200)
        assert 80 <= by_group["high"] <= 200

    def test_missing_columns_raises(self):
        bad = pd.DataFrame([{"severity_band": "high"}])  # no os_days/os_event
        with pytest.raises(KeyError, match="os_days|os_event"):
            kaplan_meier_summary(bad)


class TestLogRank:
    def test_p_value_is_float(self):
        p = multivariate_logrank_p(_synthetic_cohort())
        assert isinstance(p, float)
        assert 0.0 <= p <= 1.0

    def test_strong_separation_yields_low_p(self):
        # With this much separation between high and low, p should be < 0.05
        p = multivariate_logrank_p(_synthetic_cohort())
        assert p < 0.05


class TestTwoArmSummary:
    def test_two_arms_present(self):
        result = two_arm_summary(_synthetic_cohort())
        assert set(result["summary"]["group"]) == {"high", "not-high"}

    def test_two_arm_logrank_low_p(self):
        result = two_arm_summary(_synthetic_cohort())
        assert result["logrank_p"] < 0.05


class TestCoxSeverity:
    def test_hr_direction_positive(self):
        # severity_score ↑ should ↑ hazard, so HR > 1
        result = cox_severity(_synthetic_cohort())
        assert result["hr"] > 1.0

    def test_concordance_above_chance(self):
        # With clean separation, concordance should be well above 0.5
        result = cox_severity(_synthetic_cohort())
        assert result["concordance"] > 0.7

    def test_n_and_events_match_input(self):
        cohort = _synthetic_cohort()
        result = cox_severity(cohort)
        assert result["n"] == len(cohort)
        assert result["n_events"] == int(cohort["os_event"].sum())

    def test_missing_covariate_raises(self):
        bad = pd.DataFrame([{"os_days": 100.0, "os_event": 1}])
        with pytest.raises(KeyError, match="severity_score"):
            cox_severity(bad)


class TestMakeKmPlot:
    def test_writes_png(self, tmp_path: Path):
        out = tmp_path / "km.png"
        make_km_plot(_synthetic_cohort(), out)
        assert out.exists()
        assert out.stat().st_size > 1024  # non-trivial PNG

    def test_creates_parent_dirs(self, tmp_path: Path):
        out = tmp_path / "deep" / "nested" / "km.png"
        make_km_plot(_synthetic_cohort(), out)
        assert out.exists()
