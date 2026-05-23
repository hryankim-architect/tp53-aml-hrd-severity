"""Patient-level cohort selection for the TP53-HRD severity demo.

The goal is to assemble a *small, balanced, deterministic* cohort that
exercises the full pipeline (tiering → severity score → survival) without
inflating ``n`` past the open-tier ceiling discovered on 2026-05-23.

* **n_mut**: up to ``n_max_mut`` TP53-mutant patients with usable survival data.
  Saturday's analysis found 8 such patients across all 153 aliquot MAFs.
* **n_wt**: a deterministic random sample of WT controls drawn from patients
  with no TP53 calls and with usable survival data.
* **Seed**: ``random.seed(seed)`` before WT sampling. The same seed produces
  the same cohort on any host, which is the reproducibility contract the
  README promises.

The cohort manifest (one row per patient) is the input both for the severity
score and the survival analysis.
"""

from __future__ import annotations

import random
from dataclasses import dataclass

import pandas as pd

from tp53_hrd.annotate import patient_tier
from tp53_hrd.maf import patient_from_barcode

DEFAULT_SEED = 42
DEFAULT_N_MAX_MUT = 8
DEFAULT_N_WT = 8


@dataclass(frozen=True)
class CohortSummary:
    """Counts surfaced for logging / audit emit."""

    n_total: int
    n_mut: int
    n_wt: int
    tier_counts: dict[str, int]


def select_cohort(
    maf: pd.DataFrame,
    clinical: pd.DataFrame,
    *,
    seed: int = DEFAULT_SEED,
    n_max_mut: int = DEFAULT_N_MAX_MUT,
    n_wt: int = DEFAULT_N_WT,
) -> tuple[pd.DataFrame, CohortSummary]:
    """Assemble the cohort manifest plus a summary count.

    The MAF must already have a ``patient_id`` column (use
    :func:`patient_from_barcode` to derive it). The clinical DataFrame must
    have ``patient_id``, ``os_days``, ``os_event``.

    Returns a (cohort, summary) tuple. ``cohort`` columns:

    * ``patient_id``
    * ``group`` — ``"TP53-mut"`` or ``"WT"``
    * ``tier`` — ``"A"``, ``"B"``, ``"C"``, or ``"WT"``
    * ``os_days``, ``os_event`` (joined from clinical)
    * ``age_at_diagnosis`` (joined from clinical if present)
    """
    if "patient_id" not in maf.columns:
        raise KeyError("MAF must have 'patient_id' column — call patient_from_barcode first")
    for col in ("patient_id", "os_days", "os_event"):
        if col not in clinical.columns:
            raise KeyError(f"clinical missing required column: {col}")

    # --- 1. Patient-level TP53 tier
    tp53 = maf[maf["Hugo_Symbol"] == "TP53"]
    tp53_patient_to_tier: dict[str, str] = {}
    for patient, group in tp53.groupby("patient_id"):
        tp53_patient_to_tier[patient] = patient_tier(group)

    # --- 2. TP53-mut patients with usable survival
    clinical_with_survival = clinical[clinical["os_days"].notna()].copy()
    eligible_mut = sorted(
        p for p in tp53_patient_to_tier if p in set(clinical_with_survival["patient_id"])
    )
    selected_mut = eligible_mut[:n_max_mut]

    # --- 3. WT pool — patients with NO TP53 call in the MAF, with survival
    all_patients_in_maf = set(maf["patient_id"].dropna())
    tp53_mut_patients = set(tp53_patient_to_tier)
    wt_pool = sorted(
        p
        for p in clinical_with_survival["patient_id"]
        if p in all_patients_in_maf and p not in tp53_mut_patients
    )

    rng = random.Random(seed)
    n_wt_actual = min(n_wt, len(wt_pool))
    selected_wt = rng.sample(wt_pool, n_wt_actual)

    # --- 4. Build cohort DataFrame
    records: list[dict] = []
    for p in selected_mut:
        records.append(
            {"patient_id": p, "group": "TP53-mut", "tier": tp53_patient_to_tier[p]}
        )
    for p in selected_wt:
        records.append({"patient_id": p, "group": "WT", "tier": "WT"})

    cohort = pd.DataFrame(records)
    cohort = cohort.merge(
        clinical_with_survival[["patient_id", "os_days", "os_event", "age_at_diagnosis"]]
        if "age_at_diagnosis" in clinical_with_survival.columns
        else clinical_with_survival[["patient_id", "os_days", "os_event"]],
        on="patient_id",
        how="left",
    )

    # --- 5. Summary
    tier_counts: dict[str, int] = {}
    for t in cohort["tier"]:
        tier_counts[t] = tier_counts.get(t, 0) + 1

    summary = CohortSummary(
        n_total=len(cohort),
        n_mut=len(selected_mut),
        n_wt=len(selected_wt),
        tier_counts=tier_counts,
    )

    return cohort, summary


def cohort_to_tsv(cohort: pd.DataFrame, out_path) -> None:
    """Write the cohort manifest as TSV for committing to ``tests/fixtures``."""
    out_path = str(out_path)
    cohort.to_csv(out_path, sep="\t", index=False)


# Re-export for caller convenience
__all__ = [
    "CohortSummary",
    "DEFAULT_N_MAX_MUT",
    "DEFAULT_N_WT",
    "DEFAULT_SEED",
    "cohort_to_tsv",
    "patient_from_barcode",
    "select_cohort",
]
