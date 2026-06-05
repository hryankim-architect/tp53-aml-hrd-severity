# `tp53-aml-hrd-severity`

![ci](https://github.com/hryankim-architect/tp53-aml-hrd-severity/actions/workflows/ci.yml/badge.svg)

> **One principle, applied here.** Pick the smallest, most interpretable representation that could carry the signal; measure it against an honest baseline; report the verdict faithfully — whether the compact choice wins, ties, or loses. *That last step is why AI safety is needed: knowing a capability is real rather than a flattering benchmark.*
>
> In this repo: **representation** a bounded composite (tier+VAF, HRD-norm), nothing learned → **baseline** a many-parameter / learned Cox → **verdict** *by design* there is nothing to over-fit at n=15; the Cox HR 8.39 (p=0.024) is emitted but labeled descriptive.

*This repo uses the open-tier TCGA-LAML subset (n=15, 7 TP53-mutant + 8 WT controls); it demonstrates the method end-to-end, not production-scale statistical power.*

**What this shows (v0.2)**: TP53-mutation-driven **plus** HRD-genomic-scar
severity scoring for AML, end-to-end from open-tier TCGA-LAML mutation
calls + ASCAT2 allele-specific copy-number segments through a per-patient
composite severity score (TP53 axis × HRD-scar axis) into a Cox /
Kaplan-Meier survival readout including univariate, bivariate, and
interaction-term models. The HRD-scar component follows the Telli et al.
2016 (Clin Cancer Res, PMC6773427) definition:

> **HRD score = LOH + TAI + LST**
>
> - LOH = segments with Minor_Copy_Number == 0, length > 15 Mb,
>   NOT spanning the whole chromosome (Abkevich 2012)
> - TAI = subchromosomal segments with allelic imbalance extending to
>   a telomere, not crossing the centromere (Birkbak 2012)
> - LST = chromosomal transitions between adjacent segments ≥ 10 Mb
>   long, after smoothing across gaps < 3 Mb (Popova 2012)
> - HRD-positive call at score ≥ 42 (Telli 2016 TNBC validation;
>   AML-specific cutoff TBD)

Implemented in pure Python (no R / scarHRD dependency) over the GDC
TCGA-LAML ASCAT2 open-tier allele-specific segment files.

**v0.3 adds three interpretive analyses** on top of v0.2's Cox table,
so the README's bivariate / interaction results are now framed by
formal model selection + a mechanism-decomposition test:

- **Nested-model comparison** (`src/tp53_hrd/model_compare.py`),
  LRT + AIC delta + C-index delta for the three pairwise comparisons
  among M1 (univariate TP53), M2 (bivariate), M3 (with interaction).
  Per-run boolean `justifies_complex` for each pair.
- **Causal mediation analysis** (`src/tp53_hrd/mediation.py`),
  Baron-Kenny 1986 path decomposition with 1000-bootstrap CI on the
  indirect effect; reports `proportion_mediated` directly.
- **Clonal hierarchy approximation** (`src/tp53_hrd/clonal.py`),
  VAF-rank-based "founder vs subclonal" call per TP53-mutant
  patient; cohort-level `founder_rate_among_tp53_mut` summary.

These three together let the README's TP53→HRD upstream-downstream
narrative be reported as *converging-evidence triangulation* rather
than hand-waving from a single Cox table. Scope: at n=15 each
analysis is descriptive, not confirmatory; see `docs/release-notes/v0.3.md`.

**Reproducibility**: `make data && make run` produces the demo output in under
two minutes on a single Mac/Linux box. No GPU, no cloud credentials.

**Substrate**: emits audit (NDJSON hash-chained ledger), tracks MLflow runs,
and exposes a canary smoke test that the `lab_semantic_check.py` probe can call.

**Prior work**: A version of this method ran at full cohort scale on TCGA +
BeatAML + internal cohorts during my time directing clinical bioinformatics at
Gilead, calibrated against scarHRD signatures. This repo covers the open-shareable
TCGA-LAML slice only; the numbers reflect that subset, not the full production
run. See [`docs/what-is-out-of-scope.md`](docs/what-is-out-of-scope.md).

---

## The clinical question

Acute myeloid leukemia patients with TP53 mutations have notoriously short
overall survival, and a subset behaves like a homologous-recombination-deficient
(HRD) phenotype, even though AML is not the disease that defined the HRD
concept. The hypothesis this method codifies:

> A composite of TP53 variant *tier* (canonical hotspot vs non-hotspot missense
> vs truncating LOF) and the *clonality* of the variant (VAF as a bi-allelic
> proxy) should rank patients by HRD-like severity, and that ranking should
> separate overall survival.

The pipeline computes the score, bins it into low / moderate / high bands,
and runs both a 3-band Kaplan-Meier with multivariate log-rank and a Cox
proportional-hazards regression on the continuous score.

---

## End-to-end pipeline

```
data/tcga-laml/mafs/*.maf.gz           data/tcga-laml/clinical.json
        │                                       │
        ▼                                       ▼
load_combined_maf                       load_or_fetch_clinical
        │                                       │
        ▼                                       │
patient_from_barcode (collapse aliquots)        │
        │                                       │
        └───────────────┬───────────────────────┘
                        ▼
                  select_cohort (seed=42)
                        │
                        ▼
                 compute_severity (tier + VAF bonus)
                        │
       ┌────────────────┼──────────────────────────┐
       ▼                ▼                          ▼
  per_patient_records   make_km_plot          cox_severity +
   (JSON output)        (PNG, 3-band KM)      logrank tests
       │                │                          │
       └────────────────┴──────────────────────────┘
                        ▼
                  artifacts/
                  ├── cohort-15-results.json
                  ├── km-severity-bands.png
                  └── survival_summary.json
```

Every stage emits a NDJSON audit entry to `audit/local-demo.ndjson`. If
`AUDIT_HOST` is set, entries also POST to the substrate audit-API. MLflow
metrics flow to `MLFLOW_TRACKING_URI` if configured. Both default to no-ops,
so the demo runs cleanly on a fresh checkout.

---

## Quickstart

```bash
# 1. Install pinned dependencies
make install                  # or: uv sync --extra dev

# 2. Fetch TCGA-LAML aliquot MAFs (153 files, ~1.5 MB tarball) and clinical
make data                     # populates data/tcga-laml/

# 3. Run the end-to-end pipeline
make run                      # writes 3 files to artifacts/, < 2 seconds

# 4. Run the test suite (~90 tests, includes fixture-based integration)
make test

# 5. Run the canary smoke test (used by lab_semantic_check.py)
make canary
```

---

## Real-data results (the climax)

Pipeline output on the n=15 TCGA-LAML cohort (7 TP53-mutant + 8 wild-type
controls, seed=42 deterministic selection):

| Cohort metric | Value |
|---|---|
| TP53-mutant patients (eligible) | 7 |
| WT controls | 8 |
| Tier distribution | 1A · 3B · 3C · 8 WT |
| Severity band distribution | 5 high · 2 moderate · 8 low |

Kaplan-Meier per severity band:

| Band | n | Events | Median OS (days) | 1-yr surv | 3-yr surv |
|---|---|---|---|---|---|
| low  | 8 | 7 | 822 | 0.625 | 0.375 |
| moderate | 2 | 2 | 151 | — | — |
| high | 5 | 5 | 214 | 0.20 | — |

Statistical tests:

| Test | Value |
|---|---|
| Multivariate log-rank p (3-band) | **0.031** |
| Log-rank p (high vs not-high) | 0.065 |
| Cox HR (severity_score) | **8.39** (95% CI 1.33 – 52.94) |
| Cox p-value | **0.024** |
| Concordance index | 0.67 |

The high-band median OS is **3.8× shorter** than the low-band median (214
vs 822 days), with a statistically significant 3-band log-rank separation
even at n=15. The Cox HR direction confirms the score is monotonic with
hazard: a unit increase in `severity_score` (0 → 1) corresponds to an 8.4×
hazard increase.

These numbers are **demonstrative**, not publication-grade, n=15 is below
the cohort size most clinical studies use, but they show the method
*works on the subset of TCGA-LAML that can be fully shared in a public repo*.

---

## Per-patient JSON output (sample)

`artifacts/cohort-15-results.json` contains one record per patient. Example
for the Tier A hotspot patient (R248Q, subclonal VAF):

```json
{
  "patient_id": "TCGA-AB-2935",
  "group": "TP53-mut",
  "tier": "A",
  "variants": [
    {
      "hgvsp": "p.R248Q",
      "variant_classification": "Missense_Mutation",
      "t_alt_count": 45,
      "t_depth": 103,
      "vaf": 0.437
    }
  ],
  "max_vaf": 0.437,
  "vaf_biallelic": false,
  "severity_score": 0.75,
  "severity_band": "high",
  "os_days": 61.0,
  "os_event": 1
}
```

The full 15-record JSON is the canonical input for any downstream comparison
(e.g. plotting against external HRD scores).

---

## Sample selection, scope

The plan started at n=30 (17 TP53-mutant + 13 WT). Two open-tier ceilings
shrank the eligible cohort:

1. **MAF-tier ceiling**: TCGA-LAML's open-access MAF tarball contains
   variant calls for ~131 aliquots, of which **only 9** carry a TP53
   variant. The controlled-access tier likely has more (~20 expected from
   AML literature at ~10–15% TP53 mutation rate) but requires dbGaP
   credentials.

2. **Clinical-tier ceiling**: of those 9 TP53-mutant aliquots, **8** collapse
   to unique patients (one patient was sequenced twice as separate aliquots).
   Of those 8 patients, **1 has `vital_status=Dead` but no `days_to_death`**
   in the open clinical record. Without a death date the patient is
   unusable for survival analysis and is dropped by `clinical.parse_clinical`.

Net eligible cohort: **7 TP53-mutant patients with usable survival** +
**8 randomly selected WT controls** (seed=42) = n=15.

The pipeline still demonstrates every method step end-to-end. The README
just doesn't pretend the n is bigger than it is.

---

## Substrate environment variables

The substrate hooks read these at runtime; the defaults are no-ops, so the
demo runs cleanly without the substrate present:

| Var | Default | What it does |
|---|---|---|
| `AUDIT_HOST` | unset | If set, audit entries POST to `http://${AUDIT_HOST}/events`. |
| `MLFLOW_TRACKING_URI` | unset | If set, MLflow runs are tracked at this URI. |
| `TP53_HRD_CANARY_FIXTURE` | `tests/fixtures/canary.json` | Path used by `canary.py` for the deterministic smoke test. |

On a lab node, `scripts/run_lab.sh` exports these to the lab
defaults (`chi-mac-m:8081`, `chi-mac-m:5050`) before invoking `make run`.

---

## Repo layout

```
.
├── README.md                       # This file
├── LICENSE                         # MIT
├── Makefile                        # install | data | run | test | report | clean
├── pyproject.toml                  # uv-managed; pinned versions
├── .github/workflows/
│   └── ci.yml                      # ruff + pytest + canary
├── data/
│   ├── manifest.yaml               # (unused for P3, data is fetched dynamically)
│   └── tcga-laml/                  # populated by `make data`, git-ignored
├── src/tp53_hrd/
│   ├── audit.py                    # NDJSON hash-chained ledger emit
│   ├── tracking.py                 # MLflow run wrapper (no-op fallback)
│   ├── canary.py                   # deterministic smoke test
│   ├── maf.py                      # load + combine aliquot MAFs
│   ├── annotate.py                 # TP53 tier classification
│   ├── clinical.py                 # GDC clinical fetch + parse
│   ├── cohort.py                   # patient-level selection (seed=42)
│   ├── severity.py                 # composite tier + VAF score
│   ├── survival.py                 # KM + Cox + log-rank + plot
│   └── pipeline.py                 # end-to-end CLI entry
├── tests/                          # 87 tests across all modules
│   ├── fixtures/
│   │   ├── canary.json
│   │   └── cohort-15.tsv           # committed cohort manifest (seed=42)
│   └── test_*.py
├── docs/
│   ├── architecture.md             # substrate integration diagram
│   └── what-is-out-of-scope.md     # scope boundary ledger
└── scripts/
    └── run_lab.sh                  # one-liner for lab nodes
```

---

## What this repo does not do

See [`docs/what-is-out-of-scope.md`](docs/what-is-out-of-scope.md) for the
full ledger. Short version: no production-scale claims, no copy-number-based
LOH verification, no scarHRD R-package integration, no multi-cohort
meta-analysis, no therapy-response prediction. Those belong to the
production version of this method, not this demo.

---

## Lineage

This repo was created from
[`bioinformatics-repo-scaffold-template`](https://github.com/hryankim-architect/bioinformatics-repo-scaffold-template),
the shared scaffold that other repos in this portfolio also inherit.

---

## License

MIT. See [`LICENSE`](LICENSE).
