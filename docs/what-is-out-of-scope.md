# What is out of scope (P3 — `tp53-aml-hrd-severity`)

This file is the anti-scope-creep ledger for the P3 capability portrait.
The repo's value comes from being *small and complete* — every item below
is something a reviewer might reasonably ask for that the v0.1 demo
deliberately does not attempt.

If a future PR proposes any of these, the contributor must answer one
question: **why is this still out of scope?** If the answer is good, edit
this file in the same PR. If not, the PR doesn't land.

---

## Statistical-power claims

n=15 (7 TP53-mutant + 8 WT) is below standard clinical-study power. The
Cox HR (8.39, 95% CI 1.33–52.94, p=0.024) and the 3-band log-rank p
(0.031) are *demonstrative*, not conclusive. The wide confidence intervals
reflect this — they are the right output for the input sample size, not
an artifact to be hidden.

**Why out of scope**: Inflating the cohort would require either the
controlled-access tier (dbGaP credentials) or pooling external cohorts.
Either move expands the repo's data plane beyond the open subset, which
breaks the "small and reproducible on a single workstation" contract.

---

## Controlled-access tier expansion

The TCGA-LAML controlled-access tier likely contains ~20 TP53-mutant
patients (vs the 8 in the open tier) and richer clinical fields
(cytogenetics, blast counts, FAB subtypes). Pulling them would close the
"open-tier ceiling" gap described in the README.

**Why out of scope**: Controlled-access requires a dbGaP DAR application,
PI sign-off, and a private data plane. The capability portrait's purpose
is to be runnable by anyone who clones the repo.

---

## BeatAML extension

The BeatAML cohort (~600+ patients with WES + clinical + ex-vivo drug
response) would give the score a second cohort to calibrate against.

**Why out of scope**: BeatAML adds a second data plane, doubles the
download budget, and would force the severity score to negotiate
cross-platform calling differences. The production version of this method
already ran on BeatAML internally; the lab demo proves the method works
on one cohort.

---

## Multi-cohort survival meta-analysis

A real production HRD severity claim would fit the score on TCGA, validate
on BeatAML, and report a pooled HR with random-effects meta-analysis.

**Why out of scope**: Meta-analysis multiplies the cohort handling, the
statistical assumptions, and the "did you really beat existing scores?"
defense surface — all of which belong to a paper, not a capability
portrait.

---

## scarHRD signature calibration

scarHRD is the standard R package for genome-wide HRD signature scoring
from SNP-array or sequencing copy-number data. The production version of
this severity score was calibrated against scarHRD output as the ground
truth.

**Why out of scope**: scarHRD requires SNP-array or matched-normal WGS
inputs that are not in TCGA-LAML's open tier, plus a heavy R dependency
(R + Bioconductor) that breaks the "uv sync, single Python venv"
reproducibility promise.

---

## Copy-number-based LOH verification

The current severity score uses VAF ≥ 0.5 as a bi-allelic proxy. A real
LOH call would compare the VAF to local copy number from segmented SCNV
data, accounting for tumor purity.

**Why out of scope**: SCNV calls and tumor-purity estimates would require
a second data tier (allele-specific copy number) and a purity estimator
(e.g. ASCAT or FACETS), both of which break the "VAF as a clean proxy
that the open MAF gives us for free" contract.

---

## Co-occurring chromosomal events (chr5/7 loss, complex karyotype)

TP53-mutant AML is enriched for chr5/7 deletions and complex karyotypes,
both of which independently worsen prognosis. A real severity score would
add a +1 bonus for these.

**Why out of scope**: Cytogenetics is in TCGA-LAML controlled access; the
open clinical fetch returns vital status + age + follow-up but not
karyotype annotations. Adding this would require either the controlled-
access tier or a separate cytogenetics manifest.

---

## Therapy-response prediction

The TP53-HRD severity score *describes* prognosis under standard-of-care
induction. It does not predict response to specific drug classes (HMAs,
venetoclax combinations, PARP inhibitors).

**Why out of scope**: Response prediction needs per-regimen treatment
records and post-induction MRD endpoints, neither of which are in the
open clinical record at the granularity required.

---

## Production hardening

The pipeline runs in a single Python process. There is no HA, no RBAC,
no multi-tenant isolation, no input streaming, no retry/backoff, no
distributed orchestration.

**Why out of scope**: The substrate (`audit.py`, `tracking.py`) provides
the building blocks; the capability portrait does not re-implement them.
Production hardening belongs to the orchestration project (P1
`healthomics-lab-orchestrator`), not the analytical method demo.

---

## Adding an item

Open a PR that:

1. Adds the item to the appropriate section above (or creates a new
   section if none fits).
2. Adds a one-sentence reason in italics for why it remains out of scope.
3. Links to the upstream feature request or issue if there is one.

That's it. The friction is intentional.
