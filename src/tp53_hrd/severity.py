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


# ---------------------------------------------------------------------------
# v0.2: composite TP53 + HRD-scar severity (Telli 2016 axis added)
# ---------------------------------------------------------------------------
# The HRD-scar score (`scar.compute_hrd_score`) is a count in the
# range ~0-100+ following the LOH + TAI + LST sum convention; we
# rescale to [0, 1] using Telli's HRD-positive threshold of 42 as
# the anchor: hrd_norm = min(1.0, hrd_score / 42).
# Composite weighting (v0.2): 50% TP53, 50% HRD-scar. The choice is
# deliberately equal because (a) AML-specific weight learning would
# need n much greater than 15, and (b) equal-weight is the most
# transparent default for a v0.1 proof-of-concept scope.

HRD_NORM_ANCHOR: int = 42  # Telli 2016 HRD-positive threshold
TP53_WEIGHT: float = 0.5
HRD_WEIGHT: float = 0.5


def hrd_norm(hrd_score: int | None) -> float:
    """Normalise a Telli HRD count to [0, 1] using 42 as the unit anchor."""
    if hrd_score is None:
        return 0.0
    return min(1.0, max(0.0, hrd_score / HRD_NORM_ANCHOR))


def composite_severity(
    tp53_score: float, hrd_score_value: int | None
) -> float:
    """Composite severity = TP53_WEIGHT * tp53 + HRD_WEIGHT * hrd_norm(scar)."""
    return TP53_WEIGHT * tp53_score + HRD_WEIGHT * hrd_norm(hrd_score_value)


def compute_severity(
    cohort: pd.DataFrame,
    maf: pd.DataFrame,
    hrd_scar: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Add severity-score columns to a cohort DataFrame.

    Args:
        cohort: must have ``patient_id`` + ``tier``.
        maf: must have ``patient_id`` (call ``patient_from_barcode`` first).
        hrd_scar: optional, DataFrame indexed by patient_id with columns
            from ``scar.cohort_hrd_scores`` (loh, tai, lst, hrd_score,
            hrd_positive, n_segments_input). When provided, the v0.2
            composite TP53+HRD axis is added; otherwise the cohort
            keeps only the v0.1 TP53-only score.

    Returns a copy of ``cohort`` with new columns:

    * ``max_vaf`` / ``vaf_biallelic`` / ``severity_score`` / ``severity_band``
      — v0.1 TP53-only columns.
    * ``loh`` / ``tai`` / ``lst`` / ``hrd_score`` / ``hrd_positive`` /
      ``n_segments_input`` — v0.2 HRD-scar columns (only when ``hrd_scar``
      is provided; NaN/False for patients with no ASCAT data).
    * ``composite_severity`` — v0.2 0.5*TP53 + 0.5*HRD_norm composite
      (only when ``hrd_scar`` is provided; falls back to ``severity_score``
      when not).
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
        rec: dict[str, Any] = {
            "max_vaf": max_vaf,
            "vaf_biallelic": bool(
                max_vaf is not None and max_vaf >= VAF_BIALLELIC_THRESHOLD
            ),
            "severity_score": score,
            "severity_band": severity_band(score),
        }
        # v0.2: layer in HRD-scar if available for this patient
        if hrd_scar is not None and pid in hrd_scar.index:
            sr = hrd_scar.loc[pid]
            rec.update({
                "loh": int(sr["loh"]),
                "tai": int(sr["tai"]),
                "lst": int(sr["lst"]),
                "hrd_score": int(sr["hrd_score"]),
                "hrd_positive": bool(sr["hrd_positive"]),
                "n_segments_input": int(sr["n_segments_input"]),
                "composite_severity": composite_severity(score, int(sr["hrd_score"])),
            })
        elif hrd_scar is not None:
            # cohort patient but no ASCAT data — composite falls back to TP53-only
            rec.update({
                "loh": None,
                "tai": None,
                "lst": None,
                "hrd_score": None,
                "hrd_positive": False,
                "n_segments_input": 0,
                "composite_severity": composite_severity(score, None),
            })
        records.append(rec)

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
