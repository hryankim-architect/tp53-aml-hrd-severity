"""Tests for the v0.2 composite TP53 + HRD-scar severity score."""

from __future__ import annotations

import pandas as pd
import pytest

from tp53_hrd import severity


def test_hrd_norm_clipped_to_unit_interval():
    assert severity.hrd_norm(0) == 0.0
    assert severity.hrd_norm(21) == pytest.approx(0.5, abs=1e-6)
    assert severity.hrd_norm(42) == 1.0
    assert severity.hrd_norm(100) == 1.0
    assert severity.hrd_norm(None) == 0.0


def test_composite_severity_equal_weighting():
    # tp53=1.0 + hrd_norm=0.0 -> 0.5
    assert severity.composite_severity(1.0, 0) == pytest.approx(0.5)
    # tp53=0.0 + hrd_score=42 (norm=1.0) -> 0.5
    assert severity.composite_severity(0.0, 42) == pytest.approx(0.5)
    # both max -> 1.0
    assert severity.composite_severity(1.0, 42) == pytest.approx(1.0)
    # both zero -> 0.0
    assert severity.composite_severity(0.0, 0) == pytest.approx(0.0)


def test_compute_severity_layers_in_hrd_when_provided():
    cohort = pd.DataFrame({
        "patient_id": ["TCGA-AA-0001", "TCGA-BB-0002"],
        "tier": ["A", "WT"],
        "group": ["TP53-mut", "TP53-WT"],
    })
    maf = pd.DataFrame({
        "patient_id": ["TCGA-AA-0001"],
        "Hugo_Symbol": ["TP53"],
        "Variant_Classification": ["Missense_Mutation"],
        "HGVSp_Short": ["p.R175H"],
        "t_alt_count": [40],
        "t_depth": [80],
    })
    hrd = pd.DataFrame({
        "loh": [10, 0],
        "tai": [8, 0],
        "lst": [6, 0],
        "hrd_score": [24, 0],
        "hrd_positive": [False, False],
        "n_segments_input": [80, 75],
    }, index=["TCGA-AA-0001", "TCGA-BB-0002"])
    scored = severity.compute_severity(cohort, maf, hrd_scar=hrd)
    assert "composite_severity" in scored.columns
    assert "hrd_score" in scored.columns
    row0 = scored.iloc[0]
    # TP53 tier A + biallelic VAF 0.5 -> tp53_score = (3+1)/4 = 1.0
    # HRD = 24 -> hrd_norm = 24/42 ≈ 0.571
    # composite = 0.5*1.0 + 0.5*0.571 ≈ 0.786
    assert float(row0["composite_severity"]) == pytest.approx(0.786, abs=0.01)


def test_compute_severity_v01_only_when_no_hrd():
    cohort = pd.DataFrame({
        "patient_id": ["P1"],
        "tier": ["A"],
        "group": ["TP53-mut"],
    })
    maf = pd.DataFrame({
        "patient_id": ["P1"],
        "Hugo_Symbol": ["TP53"],
        "Variant_Classification": ["Missense_Mutation"],
        "HGVSp_Short": ["p.R175H"],
        "t_alt_count": [10],
        "t_depth": [40],  # VAF = 0.25, NOT biallelic
    })
    scored = severity.compute_severity(cohort, maf, hrd_scar=None)
    assert "composite_severity" not in scored.columns
    assert "hrd_score" not in scored.columns
    assert "severity_score" in scored.columns
