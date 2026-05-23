"""tp53_hrd: TP53-driven HRD severity score for AML (capability portrait).

A 30-patient TCGA-LAML demo of:
  - TP53 hotspot annotation (Tier A/B/C, see docs/architecture.md)
  - HRD-influenced severity score (composite of tier + VAF + co-occurring events)
  - Kaplan-Meier + Cox proportional hazards survival analysis

Substrate hooks live in :mod:`tp53_hrd.audit` (NDJSON hash-chained ledger),
:mod:`tp53_hrd.tracking` (MLflow wrapper, no-op fallback), and
:mod:`tp53_hrd.canary` (deterministic smoke test for lab_semantic_check.py).
"""

__version__ = "0.1.0"
