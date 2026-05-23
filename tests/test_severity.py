"""Tests for the TP53-HRD severity score.

Fixture values are the actual t_alt_count / t_depth pairs from the Saturday
TCGA-LAML MAF inspection — not synthetic numbers. This keeps the regression
test honest to real-data behavior.
"""

from __future__ import annotations

import pandas as pd

from tp53_hrd.severity import (
    MAX_RAW_SCORE,
    VAF_BIALLELIC_THRESHOLD,
    _safe_vaf,
    compute_severity,
    patient_max_vaf,
    severity_band,
    severity_score,
)


def _v(hgvsp: str, vc: str, alt: int, depth: int, *, hugo: str = "TP53") -> dict:
    return {
        "Hugo_Symbol": hugo,
        "HGVSp_Short": hgvsp,
        "Variant_Classification": vc,
        "t_alt_count": alt,
        "t_depth": depth,
    }


class TestSafeVaf:
    def test_normal_case(self):
        assert _safe_vaf(45, 103) == 45 / 103

    def test_zero_depth_returns_none(self):
        assert _safe_vaf(10, 0) is None

    def test_negative_depth_returns_none(self):
        assert _safe_vaf(10, -5) is None

    def test_string_inputs_coerced(self):
        assert _safe_vaf("45", "103") == 45 / 103

    def test_none_inputs(self):
        assert _safe_vaf(None, 100) is None
        assert _safe_vaf(45, None) is None


class TestPatientMaxVaf:
    def test_single_variant(self):
        df = pd.DataFrame([_v("p.R248Q", "Missense_Mutation", 45, 103)])
        vaf = patient_max_vaf(df)
        assert abs(vaf - 45 / 103) < 1e-9

    def test_two_variants_takes_max(self):
        # TCGA-AB-2938: H179R (193/316 = 0.611) vs R342Efs (50/296 = 0.169)
        df = pd.DataFrame([
            _v("p.H179R", "Missense_Mutation", 193, 316),
            _v("p.R342Efs*3", "Frame_Shift_Del", 50, 296),
        ])
        vaf = patient_max_vaf(df)
        assert abs(vaf - 193 / 316) < 1e-9

    def test_empty(self):
        assert patient_max_vaf(pd.DataFrame()) is None


class TestSeverityScore:
    def test_wt_is_zero(self):
        assert severity_score("WT", None) == 0.0
        # Note: severity_score("WT", x) for x != None is never called in practice
        # because WT patients have no TP53 variants in the MAF, so the cohort's
        # max_vaf for WT rows is always None. The function is pure math, though,
        # so a VAF bonus would still apply if someone passed it directly:
        assert severity_score("WT", 0.8) == 1 / MAX_RAW_SCORE  # = 0.25

    def test_tier_a_no_vaf_bonus(self):
        # TCGA-AB-2935 R248Q: VAF = 45/103 ≈ 0.437 (subclonal)
        score = severity_score("A", 45 / 103)
        assert score == 3 / MAX_RAW_SCORE  # = 0.75

    def test_tier_a_with_vaf_bonus(self):
        score = severity_score("A", 0.85)
        assert score == 4 / MAX_RAW_SCORE  # = 1.0

    def test_tier_b_no_vaf_bonus(self):
        # TCGA-AB-2829 R280G: VAF = 59/158 ≈ 0.373 (subclonal)
        score = severity_score("B", 59 / 158)
        assert score == 2 / MAX_RAW_SCORE  # = 0.5

    def test_tier_b_with_vaf_bonus(self):
        # TCGA-AB-2813 C176Y: VAF = 89/116 ≈ 0.767 (clonal)
        score = severity_score("B", 89 / 116)
        assert score == 3 / MAX_RAW_SCORE  # = 0.75

    def test_tier_c_at_threshold(self):
        # Boundary check: exactly 0.5 earns the bonus
        score = severity_score("C", VAF_BIALLELIC_THRESHOLD)
        assert score == 3 / MAX_RAW_SCORE


class TestSeverityBand:
    def test_low(self):
        assert severity_band(0.0) == "low"
        assert severity_band(0.24) == "low"

    def test_moderate(self):
        assert severity_band(0.25) == "moderate"
        assert severity_band(0.5) == "moderate"
        assert severity_band(0.59) == "moderate"

    def test_high(self):
        assert severity_band(0.60) == "high"
        assert severity_band(0.75) == "high"
        assert severity_band(1.0) == "high"


class TestComputeSeverity:
    def _build_cohort_and_maf(self):
        # Saturday's 7 eligible TP53-mut patients + 2 WT controls
        cohort = pd.DataFrame([
            {"patient_id": "TCGA-AB-2813", "group": "TP53-mut", "tier": "B"},
            {"patient_id": "TCGA-AB-2829", "group": "TP53-mut", "tier": "C"},
            {"patient_id": "TCGA-AB-2868", "group": "TP53-mut", "tier": "C"},
            {"patient_id": "TCGA-AB-2885", "group": "TP53-mut", "tier": "B"},
            {"patient_id": "TCGA-AB-2904", "group": "TP53-mut", "tier": "B"},
            {"patient_id": "TCGA-AB-2935", "group": "TP53-mut", "tier": "A"},
            {"patient_id": "TCGA-AB-2938", "group": "TP53-mut", "tier": "C"},
            {"patient_id": "TCGA-AB-9001", "group": "WT", "tier": "WT"},
            {"patient_id": "TCGA-AB-9002", "group": "WT", "tier": "WT"},
        ])

        maf_rows = [
            # patient_id is added via patient_from_barcode in real pipeline; here we set directly
            (_v("p.C176Y",       "Missense_Mutation", 89,  116), "TCGA-AB-2813"),
            (_v("p.R280G",       "Missense_Mutation", 59,  158), "TCGA-AB-2829"),
            (_v("p.X225_splice", "Splice_Site",       130, 404), "TCGA-AB-2829"),
            (_v("p.X126_splice", "Splice_Site",       43,  106), "TCGA-AB-2868"),
            (_v("p.H193Y",       "Missense_Mutation", 77,  93),  "TCGA-AB-2885"),
            (_v("p.R337C",       "Missense_Mutation", 87,  109), "TCGA-AB-2904"),
            (_v("p.R248Q",       "Missense_Mutation", 45,  103), "TCGA-AB-2935"),
            (_v("p.R342Efs*3",   "Frame_Shift_Del",   50,  296), "TCGA-AB-2938"),
            (_v("p.H179R",       "Missense_Mutation", 193, 316), "TCGA-AB-2938"),
        ]
        maf = pd.DataFrame([{**v, "patient_id": pid} for v, pid in maf_rows])
        return cohort, maf

    def test_seven_mut_patients_scored(self):
        cohort, maf = self._build_cohort_and_maf()
        result = compute_severity(cohort, maf)
        assert len(result) == 9
        for col in ("max_vaf", "vaf_biallelic", "severity_score", "severity_band"):
            assert col in result.columns

    def test_wt_patients_have_no_vaf(self):
        cohort, maf = self._build_cohort_and_maf()
        result = compute_severity(cohort, maf)
        wt_rows = result[result["group"] == "WT"]
        assert all(pd.isna(wt_rows["max_vaf"]))
        assert all(wt_rows["severity_score"] == 0.0)
        assert all(wt_rows["severity_band"] == "low")

    def test_tier_a_patient_2935_is_high(self):
        # R248Q VAF 0.437 → no bonus, but Tier A alone = 0.75 → high
        # (DataFrame returns numpy.bool_, so use truthiness rather than `is False`)
        cohort, maf = self._build_cohort_and_maf()
        result = compute_severity(cohort, maf)
        row = result[result["patient_id"] == "TCGA-AB-2935"].iloc[0]
        assert abs(row["max_vaf"] - 45 / 103) < 1e-9
        assert not row["vaf_biallelic"]
        assert row["severity_score"] == 0.75
        assert row["severity_band"] == "high"

    def test_tier_c_with_max_vaf_2938_is_high(self):
        # 2938 max VAF = 193/316 = 0.611 (above 0.5) → Tier C + bonus = 0.75 → high
        cohort, maf = self._build_cohort_and_maf()
        result = compute_severity(cohort, maf)
        row = result[result["patient_id"] == "TCGA-AB-2938"].iloc[0]
        assert row["vaf_biallelic"]
        assert row["severity_score"] == 0.75
        assert row["severity_band"] == "high"

    def test_tier_c_subclonal_2829_is_moderate(self):
        # 2829 max VAF = 130/404 ≈ 0.322 → no bonus → Tier C = 0.5 → moderate
        cohort, maf = self._build_cohort_and_maf()
        result = compute_severity(cohort, maf)
        row = result[result["patient_id"] == "TCGA-AB-2829"].iloc[0]
        assert not row["vaf_biallelic"]
        assert row["severity_score"] == 0.5
        assert row["severity_band"] == "moderate"

    def test_expected_band_distribution_real_data(self):
        # Saturday-discovered distribution on the real TCGA-LAML cohort:
        # high: 2813 (B+bonus), 2885 (B+bonus), 2904 (B+bonus), 2935 (A no-bonus), 2938 (C+bonus) = 5
        # moderate: 2829 (C subclonal), 2868 (C subclonal) = 2
        # low: 2 WT
        cohort, maf = self._build_cohort_and_maf()
        result = compute_severity(cohort, maf)
        counts = result["severity_band"].value_counts().to_dict()
        assert counts == {"high": 5, "low": 2, "moderate": 2}
