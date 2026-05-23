"""TP53 variant tiering for the AML HRD severity score.

Three tiers, applied per variant. When a patient has multiple TP53 variants the
*maximum* severity tier wins (see :func:`patient_tier`).

* **Tier A** — canonical hotspots in the DNA-binding domain. These are the ten
  most-recurrent TP53 variants across solid tumors and AML, listed in COSMIC
  and OncoKB as high-confidence pathogenic. See
  ``docs/architecture.md`` for the source list.
* **Tier B** — non-hotspot missense in TP53. Functional impact uncertain
  without further evidence (in v0.1 we treat them as moderately damaging).
* **Tier C** — truncating variants (nonsense, frameshift, splice-site). These
  are confidently loss-of-function and biologically equivalent to or worse
  than a single hotspot allele.

Severity order for max-tier collapse is **A > C > B** — truncating LOF is
more confidently damaging than non-hotspot missense, even though Tier A is
the most-studied class.
"""

from __future__ import annotations

from typing import Any

import pandas as pd

# Canonical AML-relevant TP53 hotspots (HGVSp short form, without "p." prefix)
TIER_A_HOTSPOTS = frozenset(
    {
        "R175H",  # DNA-binding (β-sandwich), most common gain-of-function
        "Y220C",  # structural; rescuable by PhiKan083-class drugs
        "G245S",  # DNA-binding; complex karyotype association
        "R248Q",  # DNA contact (major groove)
        "R248W",  # DNA contact (major groove)
        "R249S",  # DNA contact
        "R273H",  # DNA contact
        "R273C",  # DNA contact
        "R282W",  # DNA-binding; often with chr5/7 loss
        "V157F",  # DNA-binding
    }
)

# MAF Variant_Classification values that indicate loss-of-function. Splice_Site
# implicitly covers ±2bp at canonical donor/acceptor; deeper-intron splice
# disruptions are categorized as 'Splice_Region' and are not included here
# because their functional impact is less certain.
TRUNCATING_CLASSES = frozenset(
    {
        "Nonsense_Mutation",
        "Frame_Shift_Del",
        "Frame_Shift_Ins",
        "Splice_Site",
        "Translation_Start_Site",
    }
)

# Tier severity ranking (higher = more severe). 'C' beats 'B' because
# truncating is more confidently damaging than non-hotspot missense.
TIER_RANK = {"A": 3, "C": 2, "B": 1, "WT": 0}


def _strip_p_prefix(hgvsp: str | None) -> str:
    if not hgvsp:
        return ""
    return hgvsp[2:] if hgvsp.startswith("p.") else hgvsp


def classify_tp53_variant(row: dict[str, Any] | pd.Series) -> str | None:
    """Classify one TP53 variant row as Tier A, B, or C.

    Returns ``None`` if the row is not a TP53 variant. Non-TP53 callers should
    filter the DataFrame before calling this.
    """
    if row.get("Hugo_Symbol") != "TP53":
        return None

    hgvsp = _strip_p_prefix(row.get("HGVSp_Short"))
    vc = row.get("Variant_Classification") or ""

    if hgvsp in TIER_A_HOTSPOTS:
        return "A"
    if vc in TRUNCATING_CLASSES:
        return "C"
    # Default everything else (missense, in-frame indel, silent, etc.) to B.
    # Silent variants are unusual in somatic calls; if they appear they likely
    # reflect a low-impact call.
    return "B"


def patient_tier(tp53_variants: pd.DataFrame) -> str:
    """Collapse multiple TP53 variants per patient to a single max-severity tier.

    Empty input returns ``"WT"``.
    """
    if tp53_variants.empty:
        return "WT"

    tiers = [
        classify_tp53_variant(row)
        for _, row in tp53_variants.iterrows()
    ]
    tiers = [t for t in tiers if t is not None]
    if not tiers:
        return "WT"

    return max(tiers, key=lambda t: TIER_RANK.get(t, -1))


def tier_per_patient(maf: pd.DataFrame, patient_col: str = "patient_id") -> pd.DataFrame:
    """Return a DataFrame indexed by patient with their TP53 max tier.

    Patients without TP53 calls are not in the output; merge with the full
    clinical cohort to add ``"WT"`` rows for controls.
    """
    if patient_col not in maf.columns:
        raise KeyError(
            f"column {patient_col!r} not in MAF — call patient_from_barcode "
            f"to derive it from Tumor_Sample_Barcode before tier_per_patient"
        )

    tp53 = maf[maf["Hugo_Symbol"] == "TP53"]
    if tp53.empty:
        return pd.DataFrame(columns=[patient_col, "tier", "n_variants"])

    records = []
    for patient, group in tp53.groupby(patient_col):
        records.append(
            {
                patient_col: patient,
                "tier": patient_tier(group),
                "n_variants": len(group),
            }
        )
    return pd.DataFrame(records).sort_values(patient_col).reset_index(drop=True)
