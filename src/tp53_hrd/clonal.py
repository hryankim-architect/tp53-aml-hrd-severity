"""Clonal hierarchy approximation from per-mutation VAF.

The "TP53 is upstream of HRD" hypothesis predicts that in TP53-mutant
patients, the TP53 mutation should be a *founding clonal event* (present
in essentially all tumour cells, high VAF) rather than a late subclonal
acquisition. The classical signal for this is **VAF ordering within a
patient**: in a roughly diploid, near-pure tumour, the founding mutation
typically has the highest VAF (close to 0.5 for heterozygous; higher
with LOH); later subclonal mutations have systematically lower VAF.

This module computes, for each TP53-mutant patient in the cohort, a
crude founder-vs-subclone call based on relative VAF ranking against
the patient's other somatic mutations in the open-tier MAF:

    "founder_consistent"  — TP53 max-VAF is among the top-K mutations
    "subclonal_consistent" — TP53 max-VAF is below the cohort median
    "ambiguous"            — neither extreme

Caveats (v0.3):

    Open-tier MAF VAF is a *proxy* for clonal architecture, not a
    measurement of it. Proper clonal-evolution analysis needs either:

    * Allele-specific copy-number correction (CCF = cancer cell
      fraction, e.g. PyClone-VI / SciClone)
    * Single-cell DNA-seq or scRNA-seq
    * Longitudinal samples (diagnosis vs relapse)

    All three are out of scope for this open-tier proof of concept.
    The v0.3 output is labelled *"VAF-rank consistency"* in the README
    rather than "founder/subclone call" to keep the proxy nature
    visible to the reader.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
import pandas as pd

# Mutations with VAF below this are dropped as likely sequencing noise /
# very-low-frequency subclones we cannot reliably rank.
MIN_VAF_FOR_RANKING: float = 0.05

# Founder-consistent: TP53 max-VAF is in the top-K of the patient's
# rank-ordered mutations. K=3 is intentionally generous (most TP53-AML
# patients have <20 high-quality somatic calls in the open MAF).
TOP_K_FOR_FOUNDER: int = 3


@dataclass
class ClonalCall:
    """Per-patient VAF-rank-consistency call."""

    patient_id: str
    n_mutations_total: int
    n_mutations_ranked: int  # passed MIN_VAF_FOR_RANKING
    tp53_max_vaf: float | None
    tp53_rank_among_ranked: int | None
    median_vaf_other: float | None
    call: str  # founder_consistent | subclonal_consistent | ambiguous | no_tp53


def _safe_vaf(t_alt: object, t_depth: object) -> float | None:
    try:
        alt = float(t_alt)
        depth = float(t_depth)
    except (TypeError, ValueError):
        return None
    if depth <= 0:
        return None
    return alt / depth


def clonal_call_for_patient(patient_id: str, patient_maf: pd.DataFrame) -> ClonalCall:
    """Compute the VAF-rank-consistency call for one patient."""
    n_total = int(len(patient_maf))
    vafs = []
    tp53_vafs = []
    for _, row in patient_maf.iterrows():
        vaf = _safe_vaf(row.get("t_alt_count"), row.get("t_depth"))
        if vaf is None or vaf < MIN_VAF_FOR_RANKING:
            continue
        vafs.append(vaf)
        if str(row.get("Hugo_Symbol", "")).upper() == "TP53":
            tp53_vafs.append(vaf)

    if not tp53_vafs:
        return ClonalCall(
            patient_id=patient_id,
            n_mutations_total=n_total,
            n_mutations_ranked=len(vafs),
            tp53_max_vaf=None,
            tp53_rank_among_ranked=None,
            median_vaf_other=float(np.median(vafs)) if vafs else None,
            call="no_tp53",
        )

    tp53_max = max(tp53_vafs)
    ranked_desc = sorted(vafs, reverse=True)
    # 1-indexed rank of TP53 max within the patient's ranked mutations
    tp53_rank = next(
        (i + 1 for i, v in enumerate(ranked_desc) if abs(v - tp53_max) < 1e-9), None
    )

    non_tp53_vafs = [v for v in vafs if v not in tp53_vafs]
    median_other = float(np.median(non_tp53_vafs)) if non_tp53_vafs else None

    if tp53_rank is not None and tp53_rank <= TOP_K_FOR_FOUNDER:
        call = "founder_consistent"
    elif median_other is not None and tp53_max < median_other:
        call = "subclonal_consistent"
    else:
        call = "ambiguous"

    return ClonalCall(
        patient_id=patient_id,
        n_mutations_total=n_total,
        n_mutations_ranked=len(vafs),
        tp53_max_vaf=float(tp53_max),
        tp53_rank_among_ranked=int(tp53_rank) if tp53_rank else None,
        median_vaf_other=median_other,
        call=call,
    )


def cohort_clonal_calls(
    cohort: pd.DataFrame, maf: pd.DataFrame
) -> pd.DataFrame:
    """Apply clonal_call_for_patient across every cohort patient.

    Args:
        cohort: must have ``patient_id`` (and ideally ``tier`` for
                downstream filtering).
        maf: must have ``patient_id``, ``Hugo_Symbol``, ``t_alt_count``,
             ``t_depth``.

    Returns:
        DataFrame indexed by ``patient_id`` with the ClonalCall fields
        as columns.
    """
    required_cohort = {"patient_id"}
    required_maf = {"patient_id", "Hugo_Symbol", "t_alt_count", "t_depth"}
    missing_c = required_cohort - set(cohort.columns)
    missing_m = required_maf - set(maf.columns)
    if missing_c or missing_m:
        raise KeyError(
            f"missing: cohort needs {sorted(missing_c)}, maf needs {sorted(missing_m)}"
        )

    rows = []
    for _, row in cohort.iterrows():
        pid = row["patient_id"]
        patient_maf = maf[maf["patient_id"] == pid]
        rows.append(asdict(clonal_call_for_patient(pid, patient_maf)))
    return pd.DataFrame(rows).set_index("patient_id")


def cohort_clonal_summary(calls_df: pd.DataFrame) -> dict:
    """Aggregate the per-patient calls into a single summary dict."""
    counts = calls_df["call"].value_counts().to_dict()
    total = int(len(calls_df))
    tp53_mut_only = calls_df[calls_df["call"] != "no_tp53"]
    n_tp53_mut = int(len(tp53_mut_only))
    founder_n = int(counts.get("founder_consistent", 0))
    subclone_n = int(counts.get("subclonal_consistent", 0))
    ambig_n = int(counts.get("ambiguous", 0))

    # The headline "TP53 upstream consistent" rate: among TP53-mut patients,
    # fraction whose TP53 mutation ranks as a likely founding clonal event.
    founder_rate: float | None = founder_n / n_tp53_mut if n_tp53_mut > 0 else None

    return {
        "n_total": total,
        "n_tp53_mut": n_tp53_mut,
        "founder_consistent": founder_n,
        "subclonal_consistent": subclone_n,
        "ambiguous": ambig_n,
        "no_tp53": int(counts.get("no_tp53", 0)),
        "founder_rate_among_tp53_mut": founder_rate,
    }
