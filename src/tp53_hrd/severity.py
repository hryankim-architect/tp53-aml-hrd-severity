"""TP53-driven HRD severity score (v0.1).

The score is a deliberately small *composite* of two signals available in the
open-tier MAF data:

1. **Tier weight** — Tier A = 3 (canonical hotspot), B = 2 (non-hotspot
   missense), C = 2 (truncating LOF), WT = 0. A and {B, C} differ but B and
   C are weighted equally because truncating LOF and non-hotspot missense
   are both moderately damaging in the absence of further evidence.

2. **VAF bi-allelic bonus** — TP53 acts more like a tumor suppressor when
   both alleles are disrupted. A VAF ≥ 0.5 in a diploid tumor is a proxy
   for either clonal heterozygous *plus* LOH (loss of the WT allele) or a
   bi-allelic hit; this earns a +1 bonus. The check is intentionally crude
   (we have no copy-number data in the open MAF tier) but matches how a
   working clinical bioinformatics scientist would scope a v0.1.

Max raw score is therefore 4 (Tier A + bonus). The final ``severity_score``
is the raw score normalized to ``[0, 1]``.

Severity bands:

* ``"low"`` — score < 0.25 (essentially WT)
* ``"moderate"`` — 0.25 ≤ score < 0.60 (TP53-mut without VAF support)
* ``"high"`` — score ≥ 0.60 (Tier A regardless of VAF, OR Tier B/C with
  bi-allelic VAF support)

What's **out of scope for v0.1** (documented in ``docs/what-is-out-of-scope.md``):

* Co-occurring chromosomal events (chr5/7 loss, complex karyotype) — needs
  cytogenetics data not in the basic clinical fetch.
* Copy-number-based LOH verification — needs SCNV calls.
* scarHRD R-package integration for genome-wide HRD signatures.
"""

from __future__ import annotations

from typing import Any

import pandas as pd

from tp53_hrd.annotate import patient_tier

TIER_WEIGHTS: dict[str, int] = {"A": 3, "C": 2, "B": 2, "WT": 0}
VAF_BIALLELIC_THRESHOLD = 0.5
MAX_RAW_SCORE = max(TIER_WEIGHTS.values()) + 1  # max tier + max VAF bonus = 4


def _safe_vaf(t_alt_count: Any, t_depth: Any) -> float | None:
    """VAF = t_alt_count / t_depth, with sane handling of zeros / NaNs."""
    try:
        depth = float(t_depth)
        alt = float(t_alt_count)
    except (TypeError, ValueError):
        return None
    if depth <= 0:
        return None
    return alt / depth


def patient_max_vaf(tp53_variants: pd.DataFrame) -> float | None:
    """Return the highest VAF across all TP53 variants in a patient.

    Returns ``None`` if there are no usable VAF values (e.g. all rows have
    zero depth or missing alt counts).
    """
    if tp53_variants.empty:
        return None
    vafs = [
        v
        for v in (
            _safe_vaf(row.get("t_alt_count"), row.get("t_depth"))
            for _, row in tp53_variants.iterrows()
        )
        if v is not None
    ]
    return max(vafs) if vafs else None


def severity_score(tier: str, max_vaf: float | None) -> float:
    """Compute the normalized severity score ``[0, 1]`` from tier + VAF."""
    tier_w = TIER_WEIGHTS.get(tier, 0)
    vaf_bonus = (
        1
        if (max_vaf is not None and max_vaf >= VAF_BIALLELIC_THRESHOLD)
        else 0
    )
    return (tier_w + vaf_bonus) / MAX_RAW_SCORE


def severity_band(score: float) -> str:
    """Map a numeric severity score to a clinical band label."""
    if score < 0.25:
        return "low"
    if score < 0.60:
        return "moderate"
    return "high"


def compute_severity(
    cohort: pd.DataFrame, maf: pd.DataFrame
) -> pd.DataFrame:
    """Add severity-score columns to a cohort DataFrame.

    Returns a copy of ``cohort`` with new columns:

    * ``max_vaf`` — highest TP53 VAF across all of the patient's variants
      (``None`` for WT patients).
    * ``vaf_biallelic`` — boolean, ``True`` if ``max_vaf >= 0.5``.
    * ``severity_score`` — float in ``[0, 1]``.
    * ``severity_band`` — ``"low"``, ``"moderate"``, ``"high"``.
    """
    if "patient_id" not in cohort.columns:
        raise KeyError("cohort missing 'patient_id'")
    if "tier" not in cohort.columns:
        raise KeyError("cohort missing 'tier' — run select_cohort first")
    if "patient_id" not in maf.columns:
        raise KeyError(
            "maf missing 'patient_id' — apply patient_from_barcode first"
        )

    tp53 = maf[maf["Hugo_Symbol"] == "TP53"]

    records = []
    for _, row in cohort.iterrows():
        pid = row["patient_id"]
        tier = row["tier"]
        patient_variants = tp53[tp53["patient_id"] == pid]
        max_vaf = patient_max_vaf(patient_variants)
        score = severity_score(tier, max_vaf)
        records.append(
            {
                "max_vaf": max_vaf,
                "vaf_biallelic": bool(
                    max_vaf is not None and max_vaf >= VAF_BIALLELIC_THRESHOLD
                ),
                "severity_score": score,
                "severity_band": severity_band(score),
            }
        )

    extras = pd.DataFrame(records, index=cohort.index)
    return pd.concat([cohort.reset_index(drop=True),
                      extras.reset_index(drop=True)], axis=1)


def per_patient_records(
    cohort_with_severity: pd.DataFrame, maf: pd.DataFrame
) -> list[dict[str, Any]]:
    """Build the per-patient JSON records the README references.

    Each record is fully self-contained for downstream survival / report code
    and is the canonical output of the P3 pipeline body.
    """
    tp53 = maf[maf["Hugo_Symbol"] == "TP53"]
    records: list[dict[str, Any]] = []
    for _, row in cohort_with_severity.iterrows():
        pid = row["patient_id"]
        patient_variants = tp53[tp53["patient_id"] == pid]
        variant_list = []
        for _, v in patient_variants.iterrows():
            variant_list.append(
                {
                    "hgvsp": v.get("HGVSp_Short"),
                    "variant_classification": v.get("Variant_Classification"),
                    "t_alt_count": _to_int_or_none(v.get("t_alt_count")),
                    "t_depth": _to_int_or_none(v.get("t_depth")),
                    "vaf": _safe_vaf(v.get("t_alt_count"), v.get("t_depth")),
                }
            )

        # Re-derive tier from MAF to keep cohort and MAF strictly aligned
        # (this catches drift between the two stages of the pipeline).
        derived_tier = (
            patient_tier(patient_variants) if row["group"] == "TP53-mut" else "WT"
        )

        records.append(
            {
                "patient_id": pid,
                "group": row["group"],
                "tier": row["tier"],
                "tier_derived": derived_tier,
                "variants": variant_list,
                "max_vaf": (
                    float(row["max_vaf"])
                    if pd.notna(row.get("max_vaf"))
                    else None
                ),
                "vaf_biallelic": bool(row.get("vaf_biallelic", False)),
                "severity_score": float(row["severity_score"]),
                "severity_band": row["severity_band"],
                "os_days": (
                    float(row["os_days"])
                    if pd.notna(row.get("os_days"))
                    else None
                ),
                "os_event": (
                    int(row["os_event"])
                    if pd.notna(row.get("os_event"))
                    else None
                ),
            }
        )
    return records


def _to_int_or_none(value: Any) -> int | None:
    try:
        return int(value) if pd.notna(value) else None
    except (TypeError, ValueError):
        return None
