# Architecture

This repo uses a small, deliberate architecture. The goal is that a reviewer
can trace the full pipeline from raw GDC data to survival output in a single
sitting.

## Pipeline stages

| Stage | Module | Audit event | Primary artifact |
|---|---|---|---|
| Entry point | `tp53_hrd.pipeline` | `pipeline_start` | — |
| MAF load + combine | `tp53_hrd.maf` | `maf_loaded` | in-memory DataFrame |
| TP53 tier annotation | `tp53_hrd.annotate` | `annotate_complete` | tier labels per variant |
| Clinical fetch / parse | `tp53_hrd.clinical` | `clinical_loaded` | patient OS fields |
| Cohort selection (seed=42) | `tp53_hrd.cohort` | `cohort_selected` | `tests/fixtures/cohort-15.tsv` |
| HRD scar scoring (LOH+TAI+LST) | `tp53_hrd.scar` | `hrd_scored` | per-patient HRD score |
| Composite severity scoring | `tp53_hrd.severity` | `severity_scored` | `severity_score`, `severity_band` |
| Survival analysis (KM + Cox) | `tp53_hrd.survival` | `survival_complete` | `artifacts/survival_summary.json` |
| Artifact write | `tp53_hrd.pipeline` | `pipeline_complete` | `artifacts/cohort-15-results.json`, `artifacts/km-severity-bands.png` |

Each row in the table maps to one or more `audit.emit()` calls. The NDJSON
ledger at `audit/local-demo.ndjson` is the persistent record of every run.

## Substrate integration

Three loosely-coupled channels connect this repo to the Polish-Phase5 substrate:

| Channel | Module | Env var | Substrate endpoint |
|---|---|---|---|
| Audit (immutable record) | `tp53_hrd.audit` | `AUDIT_HOST` | `http://${AUDIT_HOST}/events` |
| Experiment tracking | `tp53_hrd.tracking` | `MLFLOW_TRACKING_URI` | configurable |
| Canary smoke test | `tp53_hrd.canary` | `TP53_HRD_CANARY_FIXTURE` | invoked by `lab_semantic_check.py` |

All three degrade to no-ops when the substrate is absent. The local NDJSON
ledger remains the source of truth for audit even if the remote POST fails.

## Hash-chained audit ledger

Each NDJSON entry carries a `prev_hash` field: the SHA-256 of the canonical
(sorted keys, controlled separators) JSON of the previous entry. Any
modification to a past entry invalidates the hash of every entry that follows.
The `audit.verify()` function walks the chain and returns `(ok, n_entries,
first_bad_ts)`.

On the Polish-Phase5 substrate this verification runs at roughly 6.19 µs per
entry up to 10 k entries, with full-chain tamper detection at about 6 ms. This
repo's audit volume is far smaller; it uses the same format so the substrate's
`gatk_audit.py` verifier works against it without modification.

## MLflow wrapper

The `tp53_hrd.tracking` module wraps an MLflow active-run context. Three
reasons it exists:

1. Parameters and metrics are version-controlled alongside the run, so the
   demo output is reproducible.
2. When the substrate MLflow server is reachable, a reviewer can compare runs
   across projects from a single UI.
3. When the server is absent, the wrapper is a no-op and `make run` completes
   cleanly on a plain laptop.

## Canary

`tp53_hrd.canary` is a fixture-driven smoke test. It loads a committed
JSON fixture, runs the core scoring path, and exits 0 if output matches
the expected values. Requirements: deterministic input, completes in under
30 seconds, no external services.

The Polish-Phase5 `lab_semantic_check.py` probe calls this daily. A green
canary signals that the substrate-level monitoring can detect regressions
in this repo without custom alerting code here.

## What this architecture avoids

- No microservices.
- No async runtime.
- No process supervisor.
- No container per pipeline stage (everything runs in a single Python process).
- No schema-validation layer; Pydantic is used only where a check is cheap.
- No DAG engine; if a future version needs one, it belongs inside the body
  of the pipeline, not in the scaffold shared with other repos.

The scaffold defines the contract. Implementation details belong to the
individual pipeline modules.
