"""Tests for TP53 variant tiering.

Fixture variants are the *real* TP53 calls discovered in the TCGA-LAML
open-tier MAFs on 2026-05-23 (n=11 variants across 8 unique patients). This
keeps the test suite anchored to actual data instead of synthetic edge cases.
"""

from __future__ import annotations

import pandas as pd

from tp53_hrd.annotate import (
    TIER_A_HOTSPOTS,
    TIER_RANK,
    classify_tp53_variant,
    patient_tier,
    tier_per_patient,
)

# (patient_id, HGVSp_Short, Variant_Classification, expected_tier)
SATURDAY_VARIANTS = [
    ("TCGA-AB-2813", "p.C176Y",       "Missense_Mutation", "B"),
    ("TCGA-AB-2829", "p.R280G",       "Missense_Mutation", "B"),
    ("TCGA-AB-2829", "p.X225_splice", "Splice_Site",       "C"),
    ("TCGA-AB-2868", "p.X126_splice", "Splice_Site",       "C"),
    ("TCGA-AB-2885", "p.H193Y",       "Missense_Mutation", "B"),
    ("TCGA-AB-2904", "p.R337C",       "Missense_Mutation", "B"),
    ("TCGA-AB-2935", "p.R248Q",       "Missense_Mutation", "A"),  # hotspot
    ("TCGA-AB-2935", "p.R248Q",       "Missense_Mutation", "A"),  # duplicate aliquot
    ("TCGA-AB-2938", "p.R342Efs*3",   "Frame_Shift_Del",   "C"),
    ("TCGA-AB-2938", "p.H179R",       "Missense_Mutation", "B"),
    ("TCGA-AB-2943", "p.R273C",       "Missense_Mutation", "A"),  # hotspot
]


def _row(hgvsp: str, vc: str, *, hugo: str = "TP53") -> dict:
    return {"Hugo_Symbol": hugo, "HGVSp_Short": hgvsp, "Variant_Classification": vc}


class TestClassifyVariant:
    def test_hotspot_r248q_is_tier_a(self):
        assert classify_tp53_variant(_row("p.R248Q", "Missense_Mutation")) == "A"

    def test_hotspot_r273c_is_tier_a(self):
        assert classify_tp53_variant(_row("p.R273C", "Missense_Mutation")) == "A"

    def test_non_hotspot_missense_is_tier_b(self):
        assert classify_tp53_variant(_row("p.C176Y", "Missense_Mutation")) == "B"
        assert classify_tp53_variant(_row("p.H193Y", "Missense_Mutation")) == "B"

    def test_splice_site_is_tier_c(self):
        assert classify_tp53_variant(_row("p.X126_splice", "Splice_Site")) == "C"

    def test_frameshift_is_tier_c(self):
        assert classify_tp53_variant(_row("p.R342Efs*3", "Frame_Shift_Del")) == "C"

    def test_non_tp53_returns_none(self):
        assert classify_tp53_variant(_row("p.V600E", "Missense_Mutation", hugo="BRAF")) is None

    def test_handles_missing_hgvsp(self):
        # Splice variants may have null HGVSp_Short; should fall back to vc
        assert classify_tp53_variant(_row(None, "Splice_Site")) == "C"

    def test_handles_missing_p_prefix(self):
        # Some annotators omit the "p." prefix
        assert classify_tp53_variant(_row("R248Q", "Missense_Mutation")) == "A"

    def test_all_ten_canonical_hotspots_recognized(self):
        for hp in TIER_A_HOTSPOTS:
            assert classify_tp53_variant(_row(f"p.{hp}", "Missense_Mutation")) == "A"


class TestPatientTier:
    def test_empty_is_wt(self):
        assert patient_tier(pd.DataFrame()) == "WT"

    def test_single_hotspot_is_a(self):
        df = pd.DataFrame([_row("p.R248Q", "Missense_Mutation")])
        assert patient_tier(df) == "A"

    def test_missense_plus_splice_is_c_over_b(self):
        # TCGA-AB-2829 case from Saturday data
        df = pd.DataFrame([
            _row("p.R280G", "Missense_Mutation"),
            _row("p.X225_splice", "Splice_Site"),
        ])
        assert patient_tier(df) == "C"

    def test_frameshift_plus_missense_is_c_over_b(self):
        # TCGA-AB-2938 case from Saturday data
        df = pd.DataFrame([
            _row("p.R342Efs*3", "Frame_Shift_Del"),
            _row("p.H179R", "Missense_Mutation"),
        ])
        assert patient_tier(df) == "C"

    def test_hotspot_beats_everything(self):
        df = pd.DataFrame([
            _row("p.R248Q", "Missense_Mutation"),   # A
            _row("p.X1_splice", "Splice_Site"),     # C
            _row("p.S100T", "Missense_Mutation"),   # B
        ])
        assert patient_tier(df) == "A"

    def test_tier_rank_ordering(self):
        # A is most severe, then C, then B, then WT
        assert TIER_RANK["A"] > TIER_RANK["C"] > TIER_RANK["B"] > TIER_RANK["WT"]


class TestTierPerPatient:
    def _build_maf(self) -> pd.DataFrame:
        """Build a small MAF DataFrame from the Saturday variant list."""
        rows = []
        for pid, hgvsp, vc, _expected in SATURDAY_VARIANTS:
            rows.append({
                "Hugo_Symbol": "TP53",
                "HGVSp_Short": hgvsp,
                "Variant_Classification": vc,
                "patient_id": pid,
            })
        # Add some non-TP53 noise — should be filtered out by tier_per_patient
        rows.append({
            "Hugo_Symbol": "DNMT3A",
            "HGVSp_Short": "p.R882H",
            "Variant_Classification": "Missense_Mutation",
            "patient_id": "TCGA-AB-2813",
        })
        return pd.DataFrame(rows)

    def test_eight_patients_with_known_tiers(self):
        maf = self._build_maf()
        result = tier_per_patient(maf)
        assert len(result) == 8, f"expected 8 unique TP53-mut patients, got {len(result)}"

        expected = {
            "TCGA-AB-2813": "B",
            "TCGA-AB-2829": "C",   # missense + splice → C
            "TCGA-AB-2868": "C",
            "TCGA-AB-2885": "B",
            "TCGA-AB-2904": "B",
            "TCGA-AB-2935": "A",   # R248Q hotspot
            "TCGA-AB-2938": "C",   # frameshift + missense → C
            "TCGA-AB-2943": "A",   # R273C hotspot
        }
        actual = dict(zip(result["patient_id"], result["tier"], strict=True))
        assert actual == expected

    def test_distribution_2A_3B_3C(self):
        # The Saturday-discovered distribution that anchors the README
        maf = self._build_maf()
        result = tier_per_patient(maf)
        counts = result["tier"].value_counts().to_dict()
        assert counts == {"A": 2, "B": 3, "C": 3}

    def test_missing_patient_col_raises(self):
        maf = pd.DataFrame([{"Hugo_Symbol": "TP53", "HGVSp_Short": "p.R248Q",
                             "Variant_Classification": "Missense_Mutation"}])
        try:
            tier_per_patient(maf)
        except KeyError as e:
            assert "patient_id" in str(e)
        else:
            raise AssertionError("expected KeyError for missing patient_id column")
