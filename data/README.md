# `data/`

This directory is **not** for committed data. Public inputs are downloaded
on demand via `make data`, which reads `manifest.yaml` and verifies each
file's SHA-256.

- `manifest.yaml`, required, schema documented inline.
- `*.fastq.gz`, `*.vcf.gz`, `*.tsv` etc. — git-ignored.

If you need to ship a tiny fixture inside the repo (e.g. for tests), put it
under `tests/fixtures/`, not here.

Keep inputs small and necessary: **every entry in `manifest.yaml` needs a
reason**. Adding a new input is a PR-sized decision, not a casual edit.
