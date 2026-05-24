"""Tests for the v0.3 VAF-rank clonal hierarchy approximation."""

from __future__ import annotations

import pandas as pd

from tp53_hrd.clonal import (
    clonal_call_for_patient,
    cohort_clonal_calls,
    cohort_clonal_summary,
)


def _maf_row(symbol: str, alt: int, depth: int) -> dict:
    return {
        "Hugo_Symbol": symbol,
        "t_alt_count": alt,
        "t_depth": depth,
        "Variant_Classification": "Missense_Mutation",
    }


def test_clonal_call_founder_when_tp53_is_highest_vaf():
    maf = pd.DataFrame([
        _maf_row("TP53", 45, 100),    # VAF 0.45 — highest
        _maf_row("DNMT3A", 30, 100),  # VAF 0.30
        _maf_row("NPM1", 15, 100),    # VAF 0.15
    ])
    call = clonal_call_for_patient("P1", maf)
    assert call.call == "founder_consistent"
    assert call.tp53_max_vaf == 0.45
    assert call.tp53_rank_among_ranked == 1


def test_clonal_call_subclonal_when_tp53_is_lowest_vaf():
    maf = pd.DataFrame([
        _maf_row("DNMT3A", 45, 100),
        _maf_row("FLT3", 40, 100),
        _maf_row("NPM1", 35, 100),
        _maf_row("IDH1", 30, 100),
        _maf_row("RUNX1", 25, 100),
        _maf_row("TP53", 10, 100),   # VAF 0.10 — last, below median of others
    ])
    call = clonal_call_for_patient("P2", maf)
    assert call.call == "subclonal_consistent"
    assert call.tp53_max_vaf == 0.10


def test_clonal_call_no_tp53_when_patient_has_no_tp53_mutation():
    maf = pd.DataFrame([
        _maf_row("DNMT3A", 45, 100),
        _maf_row("NPM1", 30, 100),
    ])
    call = clonal_call_for_patient("P3", maf)
    assert call.call == "no_tp53"
    assert call.tp53_max_vaf is None


def test_clonal_call_ambiguous_between_extremes():
    # TP53 ranks 4 of 5; not top-K (K=3), not below median-of-others
    maf = pd.DataFrame([
        _maf_row("DNMT3A", 50, 100),  # 0.50
        _maf_row("NPM1",   45, 100),  # 0.45
        _maf_row("FLT3",   40, 100),  # 0.40
        _maf_row("TP53",   35, 100),  # 0.35 -> rank 4
        _maf_row("RUNX1",  30, 100),  # 0.30
    ])
    call = clonal_call_for_patient("P4", maf)
    # median of others = median(0.50, 0.45, 0.40, 0.30) = 0.425; TP53 (0.35) < 0.425
    # So subclonal_consistent, not ambiguous
    assert call.call == "subclonal_consistent"


def test_clonal_call_drops_below_min_vaf():
    maf = pd.DataFrame([
        _maf_row("TP53", 50, 100),    # 0.50 — kept
        _maf_row("DNMT3A", 2, 100),   # 0.02 — dropped (below MIN_VAF_FOR_RANKING)
        _maf_row("NPM1", 3, 100),     # 0.03 — dropped
    ])
    call = clonal_call_for_patient("P5", maf)
    # Only TP53 survives ranking -> founder_consistent
    assert call.n_mutations_ranked == 1
    assert call.tp53_rank_among_ranked == 1
    assert call.call == "founder_consistent"


def test_cohort_clonal_calls_round_trip():
    cohort = pd.DataFrame({
        "patient_id": ["P1", "P2"],
        "tier": ["A", "WT"],
    })
    maf = pd.DataFrame([
        {"patient_id": "P1", "Hugo_Symbol": "TP53", "t_alt_count": 45, "t_depth": 100},
        {"patient_id": "P1", "Hugo_Symbol": "DNMT3A", "t_alt_count": 30, "t_depth": 100},
        {"patient_id": "P2", "Hugo_Symbol": "NPM1", "t_alt_count": 40, "t_depth": 100},
    ])
    df = cohort_clonal_calls(cohort, maf)
    assert "P1" in df.index
    assert "P2" in df.index
    assert df.loc["P1", "call"] == "founder_consistent"
    assert df.loc["P2", "call"] == "no_tp53"


def test_cohort_clonal_summary_computes_founder_rate():
    cohort = pd.DataFrame({"patient_id": ["P1", "P2", "P3", "P4"], "tier": ["A"] * 4})
    maf = pd.DataFrame([
        # P1, P2: TP53 high-VAF -> founder
        {"patient_id": "P1", "Hugo_Symbol": "TP53", "t_alt_count": 50, "t_depth": 100},
        {"patient_id": "P1", "Hugo_Symbol": "DNMT3A", "t_alt_count": 20, "t_depth": 100},
        {"patient_id": "P2", "Hugo_Symbol": "TP53", "t_alt_count": 48, "t_depth": 100},
        {"patient_id": "P2", "Hugo_Symbol": "NPM1", "t_alt_count": 25, "t_depth": 100},
        # P3: TP53 subclonal (low VAF)
        {"patient_id": "P3", "Hugo_Symbol": "DNMT3A", "t_alt_count": 50, "t_depth": 100},
        {"patient_id": "P3", "Hugo_Symbol": "NPM1", "t_alt_count": 45, "t_depth": 100},
        {"patient_id": "P3", "Hugo_Symbol": "FLT3", "t_alt_count": 40, "t_depth": 100},
        {"patient_id": "P3", "Hugo_Symbol": "IDH1", "t_alt_count": 35, "t_depth": 100},
        {"patient_id": "P3", "Hugo_Symbol": "TP53", "t_alt_count": 5, "t_depth": 100},
        # P4: no TP53
        {"patient_id": "P4", "Hugo_Symbol": "DNMT3A", "t_alt_count": 30, "t_depth": 100},
    ])
    df = cohort_clonal_calls(cohort, maf)
    summary = cohort_clonal_summary(df)
    assert summary["n_tp53_mut"] == 3
    assert summary["founder_consistent"] == 2
    assert summary["subclonal_consistent"] == 1
    assert summary["no_tp53"] == 1
    assert summary["founder_rate_among_tp53_mut"] == 2 / 3
