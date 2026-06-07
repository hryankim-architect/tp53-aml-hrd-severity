# `data/`

This directory is **not** for committed data. Public inputs are downloaded
on demand via `make data`, which reads `manifest.yaml` and verifies each
file's SHA-256.

`manifest.yaml` is the single source of truth and pins two NIH GDC open-tier
sources (no controlled-access auth):

- **`clinical`** — one POST to the GDC `/cases` endpoint for the TCGA-LAML
  survival fields. The raw JSON is byte-stable for a fixed query, so its
  sha256 is pinned directly. Cached at `tcga-laml/clinical.json`.
- **`inputs`** — 153 per-aliquot WXS "aliquot ensemble masked" somatic MAFs,
  each fetched by GDC file UUID. GDC stores these as immutable artifacts, so
  every file's sha256 is stable across machines. Cached under `tcga-laml/mafs/`.

Everything except `manifest.yaml` and this README is git-ignored. To regenerate
the manifest after GDC reprocesses the project:

```bash
python -m tp53_hrd.pipeline refresh-manifest        # re-enumerate the file set
python -m tp53_hrd.pipeline fetch --write-checksums  # download + fill sha256
```

then review the diff and commit. Keep inputs small and necessary: **every entry
needs a reason**. Adding one is a PR-sized decision, not a casual edit.

The manifest also pins the **15 Arm-3 ASCAT2 allele-specific CN segments** (one
per seed=42 cohort patient) under `tcga-laml/ascat/`. `make data` fetches +
checksum-verifies them like the MAFs; `scar_data` reads that cache and only
falls back to a run-time GDC fetch for anything missing.

If you need a tiny fixture inside the repo (e.g. for tests), put it under
`tests/fixtures/`, not here.
