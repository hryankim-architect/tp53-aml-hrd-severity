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

## R-resident scarHRD package itself

scarHRD is the canonical R/Bioconductor implementation of LOH + TAI + LST
counting (the same Telli 2016 metrics this repo computes).

**Why out of scope**: v0.2 reimplements the Telli 2016 definition directly
in Python over the GDC ASCAT2 open-tier segment files (see `src/tp53_hrd/scar.py`),
so the R dependency is no longer needed for the demo. A side-by-side
agreement check between this Python implementation and the R scarHRD
output on a shared reference cohort (e.g. TCGA-OVCA, where scarHRD has
published reference scores) is a defensible v0.3 add — but it requires
R + Bioconductor install, breaking the "uv sync, single Python venv"
reproducibility promise. Deferred until a reviewer asks for it.

---

## AML-specific HRD threshold calibration

The HRD-positive call uses Telli 2016's threshold of ≥ 42, which was
validated on a Triple-Negative Breast Cancer cohort. The threshold has
not been formally validated in AML.

**Why out of scope**: per-cohort threshold calibration needs a labelled
ground-truth set (e.g. functional HRD assay, RAD51 foci, BRCA1/2 methylation
status) that is not available in TCGA-LAML open tier. The README states
the threshold inheritance honestly and reports per-cohort distribution
statistics so a reader can see where this cohort actually falls.

---

## Allele-specific copy-number validation against ASCAT3 / FACETS

ASCAT2 and ASCAT3 sometimes call the same chromosome differently;
FACETS is a popular alternative caller. A full validation would compute
the HRD score under all three callers and report the agreement.

**Why out of scope**: ASCAT2 is the GDC default and has the broadest
TCGA coverage; agreement-statistics analysis is a methodological side-
project. v0.2 records which workflow_type each patient's segment file
came from in the audit ledger so the comparison is straightforward to
add if requested.

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
