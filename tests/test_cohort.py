"""Tests for cohort selection."""

from __future__ import annotations

import pandas as pd
import pytest

from tp53_hrd.cohort import CohortSummary, select_cohort


def _build_maf() -> pd.DataFrame:
    """Build a small MAF with the Saturday-discovered TP53 variants.

    Mixed in are a few non-TP53 variants (so non-TP53-mutant patients exist
    in the MAF and can serve as the WT pool) and 20 'background' patients
    with no TP53 (genuine WT candidates).
    """
    rows = []

    # TP53-mutant patients (8 from Saturday)
    tp53_variants = [
        ("TCGA-AB-2813", "p.C176Y", "Missense_Mutation"),
        ("TCGA-AB-2829", "p.X225_splice", "Splice_Site"),
        ("TCGA-AB-2868", "p.X126_splice", "Splice_Site"),
        ("TCGA-AB-2885", "p.H193Y", "Missense_Mutation"),
        ("TCGA-AB-2904", "p.R337C", "Missense_Mutation"),
        ("TCGA-AB-2935", "p.R248Q", "Missense_Mutation"),  # hotspot
        ("TCGA-AB-2938", "p.R342Efs*3", "Frame_Shift_Del"),
        ("TCGA-AB-2943", "p.R273C", "Missense_Mutation"),  # hotspot
    ]
    for pid, hgvsp, vc in tp53_variants:
        rows.append({
            "Hugo_Symbol": "TP53",
            "HGVSp_Short": hgvsp,
            "Variant_Classification": vc,
            "patient_id": pid,
            "Tumor_Sample_Barcode": f"{pid}-03A",
        })

    # 20 WT patients each with one non-TP53 variant (so they appear in the MAF)
    for i in range(20):
        pid = f"TCGA-AB-{6000 + i:04d}"
        rows.append({
            "Hugo_Symbol": "DNMT3A",
            "HGVSp_Short": "p.R882H",
            "Variant_Classification": "Missense_Mutation",
            "patient_id": pid,
            "Tumor_Sample_Barcode": f"{pid}-03A",
        })

    return pd.DataFrame(rows)


def _build_clinical(n_extra_wt: int = 20) -> pd.DataFrame:
    """Build clinical for all 8 TP53-mut + 20 WT patients."""
    rows = []
    tp53_patients = [
        "TCGA-AB-2813", "TCGA-AB-2829", "TCGA-AB-2868", "TCGA-AB-2885",
        "TCGA-AB-2904", "TCGA-AB-2935", "TCGA-AB-2938", "TCGA-AB-2943",
    ]
    for i, pid in enumerate(tp53_patients):
        rows.append({
            "patient_id": pid,
            "vital_status": "Dead" if i % 2 == 0 else "Alive",
            "os_days": float(300 + i * 50),
            "os_event": 1 if i % 2 == 0 else 0,
            "age_at_diagnosis": 20000 + i * 100,
        })
    for i in range(n_extra_wt):
        pid = f"TCGA-AB-{6000 + i:04d}"
        rows.append({
            "patient_id": pid,
            "vital_status": "Alive" if i % 3 != 0 else "Dead",
            "os_days": float(500 + i * 30),
            "os_event": 1 if i % 3 == 0 else 0,
            "age_at_diagnosis": 18000 + i * 200,
        })
    return pd.DataFrame(rows)


class TestSelectCohort:
    def test_eight_tp53_mut_plus_eight_wt(self):
        maf = _build_maf()
        clinical = _build_clinical()
        cohort, summary = select_cohort(maf, clinical)

        assert summary.n_total == 16
        assert summary.n_mut == 8
        assert summary.n_wt == 8

    def test_returns_cohort_summary_dataclass(self):
        maf = _build_maf()
        clinical = _build_clinical()
        _, summary = select_cohort(maf, clinical)
        assert isinstance(summary, CohortSummary)

    def test_deterministic_with_same_seed(self):
        maf = _build_maf()
        clinical = _build_clinical()
        cohort_a, _ = select_cohort(maf, clinical, seed=42)
        cohort_b, _ = select_cohort(maf, clinical, seed=42)
        # WT membership must be identical across runs with the same seed
        wt_a = sorted(cohort_a[cohort_a["group"] == "WT"]["patient_id"])
        wt_b = sorted(cohort_b[cohort_b["group"] == "WT"]["patient_id"])
        assert wt_a == wt_b

    def test_different_seed_different_wt(self):
        maf = _build_maf()
        clinical = _build_clinical()
        cohort_a, _ = select_cohort(maf, clinical, seed=42)
        cohort_b, _ = select_cohort(maf, clinical, seed=7)
        # With enough WT pool size (20) and only 8 picks, ~always differs
        wt_a = sorted(cohort_a[cohort_a["group"] == "WT"]["patient_id"])
        wt_b = sorted(cohort_b[cohort_b["group"] == "WT"]["patient_id"])
        assert wt_a != wt_b

    def test_tier_distribution_matches_saturday(self):
        maf = _build_maf()
        clinical = _build_clinical()
        cohort, summary = select_cohort(maf, clinical)
        # Saturday's distribution: 2 A + 3 B + 3 C TP53-mut, plus 8 WT
        assert summary.tier_counts == {"A": 2, "B": 3, "C": 3, "WT": 8}

    def test_os_days_joined_from_clinical(self):
        maf = _build_maf()
        clinical = _build_clinical()
        cohort, _ = select_cohort(maf, clinical)
        # Every selected patient must have a non-null OS days
        assert cohort["os_days"].notna().all()
        assert cohort["os_event"].notna().all()

    def test_missing_patient_id_in_maf_raises(self):
        maf = pd.DataFrame([{"Hugo_Symbol": "TP53", "Tumor_Sample_Barcode": "X"}])
        clinical = _build_clinical()
        with pytest.raises(KeyError, match="patient_id"):
            select_cohort(maf, clinical)

    def test_missing_os_col_in_clinical_raises(self):
        maf = _build_maf()
        bad_clinical = pd.DataFrame([{"patient_id": "X"}])
        with pytest.raises(KeyError, match="os_days|os_event"):
            select_cohort(maf, bad_clinical)

    def test_drops_tp53_patients_without_survival(self):
        maf = _build_maf()
        clinical = _build_clinical()
        # Remove TCGA-AB-2935 (the hotspot patient) from clinical
        clinical_missing = clinical[clinical["patient_id"] != "TCGA-AB-2935"]
        _, summary = select_cohort(maf, clinical_missing)
        assert summary.n_mut == 7  # one TP53-mut dropped for missing survival

    def test_caps_at_pool_size(self):
        maf = _build_maf()
        # Only 3 WT patients in clinical
        clinical_few_wt = _build_clinical(n_extra_wt=3)
        _, summary = select_cohort(maf, clinical_few_wt)
        assert summary.n_wt == 3  # capped at pool size
        assert summary.n_mut == 8
