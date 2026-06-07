"""End-to-end TP53-HRD severity pipeline for TCGA-LAML.

The shape mirrors the scaffold template so the audit / tracking / canary
substrate hooks fire in the same shape across every
repo::

    audit_start  →  tracking_start  →  body  →  tracking_end  →  audit_end

The body wires together every P3 module:

    load_combined_maf  →  patient_from_barcode  →  load_or_fetch (clinical)
    →  select_cohort (seed=42)  →  compute_severity  →  per_patient_records
    →  kaplan_meier_summary / multivariate_logrank_p / two_arm_summary
    →  cox_severity  →  make_km_plot

Outputs written to ``--out`` (default ``artifacts/``):

* ``cohort-15-results.json`` — one record per patient with tier, VAF,
  severity score, severity band, OS days, and OS event.
* ``survival_summary.json`` — KM summary, log-rank p-values, Cox HR / 95% CI.
* ``km-severity-bands.png`` — Kaplan-Meier curves per severity band.

The pipeline is **idempotent on cached data**: ``make data`` populates
``data/tcga-laml/mafs/`` and ``data/tcga-laml/clinical.json``; ``make run``
re-derives the cohort and downstream artifacts deterministically.
"""

from __future__ import annotations

import hashlib
import json
import re
import time
import urllib.error
import urllib.request
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import click
import yaml

from tp53_hrd import audit, gdc, tracking
from tp53_hrd.clinical import load_or_fetch as load_or_fetch_clinical
from tp53_hrd.cohort import DEFAULT_SEED, select_cohort
from tp53_hrd.maf import load_combined_maf, patient_from_barcode
from tp53_hrd.severity import compute_severity, per_patient_records
from tp53_hrd.survival import (
    cox_severity,
    kaplan_meier_summary,
    make_km_plot,
    multivariate_logrank_p,
    two_arm_summary,
)

# Default data layout — see data/manifest.yaml and clinical fetch
MAF_DIR = Path("data/tcga-laml/mafs")
CLINICAL_JSON = Path("data/tcga-laml/clinical.json")
# Arm-3 ASCAT segments live under data/ too, so `make data` (manifest) can
# pre-fetch + checksum them; scar_data caches here and skips refetch on run.
ASCAT_DIR = Path("data/tcga-laml/ascat")


def _run_id(name: str) -> str:
    stamp = datetime.now(UTC).strftime("%Y%m%d-%H%M%S")
    return f"{name}-{stamp}"


def _checksum(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def _download(url: str, dest: Path, *, timeout: float = 60.0, retries: int = 4) -> None:
    """Download ``url`` to ``dest`` with a per-attempt timeout and retries.

    GDC occasionally stalls a single connection; without a timeout a bulk fetch
    of 150+ small files can hang indefinitely. We read the whole response (these
    inputs are KB-scale) and write atomically only on success.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    last: Exception | None = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as resp:
                data = resp.read()
            dest.write_bytes(data)
            return
        except (urllib.error.URLError, TimeoutError, OSError) as exc:  # noqa: PERF203
            last = exc
            time.sleep(min(2**attempt, 8))
    raise RuntimeError(f"failed to download {url} after {retries} attempts: {last}")


def _verify_or_fetch(
    *, dest: Path, expected: str | None, rel: str, fetch: Any
) -> dict[str, Any]:
    """Shared cache-check / download / verify logic for one input.

    ``fetch`` is a zero-arg callable that writes ``dest`` (URL download for MAFs,
    POST for clinical). Returns a result dict with a ``status`` of ``cached``,
    ``downloaded`` or ``checksum_mismatch``.
    """
    if dest.exists() and expected and _checksum(dest) == expected:
        return {"rel": rel, "path": str(dest), "status": "cached"}
    fetch()
    actual = _checksum(dest)
    if expected and actual != expected:
        return {
            "rel": rel,
            "path": str(dest),
            "status": "checksum_mismatch",
            "expected": expected,
            "actual": actual,
        }
    return {"rel": rel, "path": str(dest), "status": "downloaded", "sha256": actual}


def fetch_manifest(manifest_path: Path, out_dir: Path) -> dict[str, Any]:
    """Download every input declared in the manifest; verify SHA-256.

    Handles both the ``clinical`` block (one GDC ``/cases`` POST, byte-stable for
    a fixed query) and the ``inputs`` list (per-aliquot MAFs fetched by GDC file
    UUID). Re-running on populated data is a no-op for any input whose on-disk
    sha256 already matches the manifest.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    with manifest_path.open("r", encoding="utf-8") as fh:
        manifest = yaml.safe_load(fh) or {}

    result: dict[str, Any] = {}

    clinical = manifest.get("clinical")
    if clinical:
        rel = clinical["path"]
        dest = out_dir / rel

        def _fetch_clinical(dest: Path = dest) -> None:
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(gdc.fetch_clinical_raw())

        result["clinical"] = _verify_or_fetch(
            dest=dest, expected=clinical.get("sha256"), rel=rel, fetch=_fetch_clinical
        )

    inputs: list[dict[str, Any]] = []
    for entry in manifest.get("inputs", []):
        url = entry["url"]
        rel = entry["path"]
        dest = out_dir / rel
        res = _verify_or_fetch(
            dest=dest,
            expected=entry.get("sha256"),
            rel=rel,
            fetch=lambda url=url, dest=dest: _download(url, dest),
        )
        if res["status"] == "downloaded":
            res["size_bytes"] = entry.get("size_bytes")
        inputs.append(res)
    result["inputs"] = inputs

    return result


def write_manifest_checksums(manifest_path: Path, result: dict[str, Any]) -> int:
    """Write freshly-computed sha256 values back into the manifest.

    Comment- and order-preserving (line-based, like the multiqc-gate sibling).
    Fills the ``clinical`` block's sha256 and every ``inputs`` entry matched by
    ``path``. Returns the number of sha256 fields written.
    """
    by_rel = {
        r["rel"]: r["sha256"]
        for r in result.get("inputs", [])
        if r.get("status") == "downloaded" and r.get("rel") and r.get("sha256")
    }
    clin = result.get("clinical") or {}
    clin_sha = clin.get("sha256") if clin.get("status") == "downloaded" else None
    if not by_rel and not clin_sha:
        return 0

    lines = manifest_path.read_text(encoding="utf-8").splitlines(keepends=True)
    url_re = re.compile(r"^\s*-\s*url:")
    path_re = re.compile(r"^\s*path:\s*(.+?)\s*$")
    # ``.*?`` (not ``.+?``) so a blank ``sha256:`` line — the post-refresh shape —
    # is matched and filled, not just an already-populated one.
    sha_re = re.compile(r"^(\s*)sha256:\s*(.*?)\s*$")

    # Boundary: `inputs:` key starts the list; anything before it is the
    # clinical block. (Both are top-level keys in the manifest.)
    inputs_key = next(
        (i for i, ln in enumerate(lines) if re.match(r"^inputs:\s*$", ln)), len(lines)
    )
    filled = 0

    # --- clinical block: fill the sha256 that precedes `inputs:` -------------
    if clin_sha is not None:
        for j in range(inputs_key):
            if lines[j].lstrip().startswith("#"):
                continue
            ms = sha_re.match(lines[j])
            if ms:
                lines[j] = f"{ms.group(1)}sha256: {clin_sha}\n"
                filled += 1
                break

    # --- inputs: fill each `- url:` block by its path -----------------------
    starts = [
        i
        for i, ln in enumerate(lines)
        if i >= inputs_key and url_re.match(ln) and not ln.lstrip().startswith("#")
    ]
    for k, si in enumerate(starts):
        ei = starts[k + 1] if k + 1 < len(starts) else len(lines)
        rel = sha_i = None
        indent = ""
        for j in range(si, ei):
            if lines[j].lstrip().startswith("#"):
                continue
            mp, ms = path_re.match(lines[j]), sha_re.match(lines[j])
            if mp and rel is None:
                rel = mp.group(1).strip().strip('"').strip("'")
            if ms and sha_i is None:
                sha_i, indent = j, ms.group(1)
        if rel in by_rel and sha_i is not None:
            lines[sha_i] = f"{indent}sha256: {by_rel[rel]}\n"
            filled += 1

    if filled:
        manifest_path.write_text("".join(lines), encoding="utf-8")
    return filled


def refresh_manifest(manifest_path: Path) -> int:
    """Re-enumerate the GDC TCGA-LAML MAF set and rewrite the ``inputs`` block.

    Preserves the header comments and the ``clinical`` block; replaces only the
    ``inputs:`` list with the current GDC file set (sha256 left blank, to be
    filled by ``fetch --write-checksums``). Returns the input count written.
    """
    files = gdc.list_laml_maf_files()
    entries = gdc.build_maf_inputs(files)

    lines = manifest_path.read_text(encoding="utf-8").splitlines(keepends=True)
    inputs_key = next(
        (i for i, ln in enumerate(lines) if re.match(r"^inputs:\s*$", ln)), None
    )
    head = lines[: inputs_key + 1] if inputs_key is not None else [*lines, "inputs:\n"]

    body: list[str] = []
    for e in entries:
        body.append(f"  - url: {e['url']}\n")
        body.append(f"    path: {e['path']}\n")
        body.append(f"    sha256: {e['sha256']}\n")
        body.append(f"    size_bytes: {e['size_bytes']}\n")
        body.append(f"    license: {e['license']}\n")
        body.append(f'    source: "{e["source"]}"\n')

    manifest_path.write_text("".join(head) + "".join(body), encoding="utf-8")
    return len(entries)


def run_pipeline(
    run_name: str,
    out_dir: Path,
    *,
    maf_dir: Path = MAF_DIR,
    clinical_json: Path = CLINICAL_JSON,
    seed: int = DEFAULT_SEED,
) -> dict[str, Any]:
    """Run the P3 TP53-HRD pipeline end-to-end.

    Returns a summary dict with cohort counts, severity distribution, and
    survival statistics. Writes three artifacts under ``out_dir``.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    job_id = _run_id(run_name)

    if not maf_dir.exists():
        raise FileNotFoundError(
            f"MAF directory {maf_dir} not found. Run `make data` first."
        )

    audit.emit(
        action="pipeline_start",
        job_id=job_id,
        fields={
            "out_dir": str(out_dir),
            "maf_dir": str(maf_dir),
            "clinical_json": str(clinical_json),
            "seed": seed,
        },
    )

    metrics: dict[str, float] = {}
    summary: dict[str, Any] = {}

    with tracking.run(name=job_id, experiment="tp53_hrd"):
        tracking.log_params({
            "run_name": run_name,
            "seed": seed,
            "maf_dir": str(maf_dir),
        })

        # --- begin body --------------------------------------------------

        # 1. Load + combine MAFs, derive patient IDs
        t0 = time.time()
        maf = load_combined_maf(maf_dir)
        maf["patient_id"] = maf["Tumor_Sample_Barcode"].apply(patient_from_barcode)
        metrics["n_variant_rows"] = float(len(maf))
        metrics["n_aliquots"] = float(maf["Tumor_Sample_Barcode"].nunique())

        # 2. Load clinical (cached on disk by load_or_fetch)
        clinical = load_or_fetch_clinical(clinical_json)
        metrics["n_clinical_patients"] = float(len(clinical))

        # 3. Cohort selection (deterministic by seed)
        cohort, cohort_summary = select_cohort(maf, clinical, seed=seed)
        metrics["n_cohort_total"] = float(cohort_summary.n_total)
        metrics["n_cohort_mut"] = float(cohort_summary.n_mut)
        metrics["n_cohort_wt"] = float(cohort_summary.n_wt)
        for tier, count in cohort_summary.tier_counts.items():
            metrics[f"tier_{tier}_n"] = float(count)

        audit.emit(
            action="cohort_built",
            job_id=job_id,
            fields={
                "n_total": cohort_summary.n_total,
                "n_mut": cohort_summary.n_mut,
                "n_wt": cohort_summary.n_wt,
                "tier_counts": cohort_summary.tier_counts,
            },
        )

        # 4a. (v0.2) HRD genomic-scar score per Telli 2016
        # Fetch ASCAT allele-specific segments for each cohort patient
        # from GDC (open tier), compute LOH + TAI + LST per patient.
        # Graceful skip if scar module / network is unavailable so the
        # v0.1 TP53-only path still completes.
        hrd_scar_df = None
        try:
            from tp53_hrd import scar, scar_data
            # Prefer the data/ cache populated by `make data` (manifest-pinned,
            # checksum-verified); scar_data skips refetch when a patient's
            # segment file is already present, and falls back to a runtime GDC
            # fetch for any that are missing.
            ascat_dir = ASCAT_DIR
            patient_paths = scar_data.fetch_cohort_ascat(
                list(cohort["patient_id"]), ascat_dir
            )
            audit.emit(
                action="ascat_segments_fetched",
                job_id=job_id,
                fields={
                    "n_patients_with_ascat": len(patient_paths),
                    "n_patients_in_cohort": int(len(cohort)),
                },
            )
            hrd_scar_df = scar.cohort_hrd_scores(patient_paths)
            if not hrd_scar_df.empty:
                metrics["hrd_n_patients_scored"] = float(len(hrd_scar_df))
                metrics["hrd_n_positive"] = float(hrd_scar_df["hrd_positive"].sum())
                metrics["hrd_score_mean"] = float(hrd_scar_df["hrd_score"].mean())
                metrics["hrd_loh_mean"] = float(hrd_scar_df["loh"].mean())
                metrics["hrd_tai_mean"] = float(hrd_scar_df["tai"].mean())
                metrics["hrd_lst_mean"] = float(hrd_scar_df["lst"].mean())
                audit.emit(
                    action="hrd.scar_scores.computed",
                    job_id=job_id,
                    fields={
                        "n_patients_scored": int(len(hrd_scar_df)),
                        "n_hrd_positive": int(hrd_scar_df["hrd_positive"].sum()),
                        "hrd_score_mean": float(hrd_scar_df["hrd_score"].mean()),
                    },
                )
        except Exception as exc:  # noqa: BLE001 — Arm 3 isolation
            audit.emit(
                action="hrd_scar_skipped",
                job_id=job_id,
                fields={"reason": f"{type(exc).__name__}: {exc}"},
            )

        # 4b. Severity score per patient (TP53-only + composite TP53+HRD when scar available)
        scored = compute_severity(cohort, maf, hrd_scar=hrd_scar_df)
        band_counts = scored["severity_band"].value_counts().to_dict()
        for band, count in band_counts.items():
            metrics[f"band_{band}_n"] = float(count)

        # 5. Per-patient JSON output
        records = per_patient_records(scored, maf)
        results_path = out_dir / "cohort-15-results.json"
        with results_path.open("w", encoding="utf-8") as fh:
            json.dump(records, fh, indent=2, default=str)

        # 6. Survival analysis — v0.1 univariate TP53 + v0.2 univariate HRD +
        #    bivariate TP53+HRD + interaction. n=15 is severely underpowered
        #    for 3-parameter Cox; results are reported as descriptive only.
        km_3band = kaplan_meier_summary(scored).to_dict(orient="records")
        logrank_3band_p = multivariate_logrank_p(scored)
        two_arm = two_arm_summary(scored)
        cox = cox_severity(scored)

        metrics["logrank_3band_p"] = logrank_3band_p
        metrics["logrank_2arm_p"] = two_arm["logrank_p"]
        metrics["cox_hr"] = cox["hr"]
        metrics["cox_p"] = cox["p_value"]
        metrics["cox_concordance"] = cox["concordance"]

        survival_summary: dict[str, Any] = {
            "km_3band": km_3band,
            "logrank_3band_p": logrank_3band_p,
            "two_arm_high_vs_not": {
                "summary": two_arm["summary"].to_dict(orient="records"),
                "logrank_p": two_arm["logrank_p"],
            },
            "cox_severity_score": cox,
        }

        # v0.2 bivariate Cox if HRD scar data is present
        if hrd_scar_df is not None and "hrd_score" in scored.columns:
            try:
                cox_hrd = cox_severity(scored, covariate="hrd_score")
                survival_summary["cox_hrd_score"] = cox_hrd
                metrics["cox_hrd_hr"] = cox_hrd["hr"]
                metrics["cox_hrd_p"] = cox_hrd["p_value"]
            except Exception as exc:  # noqa: BLE001
                survival_summary["cox_hrd_score"] = {"skipped": str(exc)}
            try:
                from tp53_hrd.survival import cox_bivariate
                cox_biv = cox_bivariate(scored, ["severity_score", "hrd_score"])
                survival_summary["cox_bivariate_tp53_plus_hrd"] = cox_biv
                cox_int = cox_bivariate(
                    scored, ["severity_score", "hrd_score"], interaction=True
                )
                survival_summary["cox_bivariate_with_interaction"] = cox_int
            except Exception as exc:  # noqa: BLE001
                survival_summary["cox_bivariate_skipped"] = str(exc)
            audit.emit(
                action="survival.bivariate_cox.computed",
                job_id=job_id,
                fields={
                    "has_univariate_hrd": "cox_hrd_score" in survival_summary
                                          and "skipped" not in survival_summary["cox_hrd_score"],
                    "has_bivariate": "cox_bivariate_tp53_plus_hrd" in survival_summary,
                    "has_interaction": "cox_bivariate_with_interaction" in survival_summary,
                },
            )

            # ---- v0.3: formal nested-model comparison (LRT + AIC + C-index)
            try:
                from tp53_hrd.model_compare import run_nested_model_suite
                cox_dropna = scored[
                    ["severity_score", "hrd_score", "os_days", "os_event"]
                ].dropna()
                model_comp = run_nested_model_suite(cox_dropna)
                survival_summary["model_comparison_lrt_aic_cindex"] = model_comp
                # Surface key comparison numbers as top-level metrics
                m2m1 = model_comp["m2_vs_m1"]
                m3m2 = model_comp["m3_vs_m2"]
                metrics["m2_vs_m1_lrt_p"] = float(m2m1["lrt_p"])
                metrics["m2_vs_m1_aic_delta"] = float(m2m1["aic_delta"])
                metrics["m3_vs_m2_lrt_p"] = float(m3m2["lrt_p"])
                metrics["m3_vs_m2_aic_delta"] = float(m3m2["aic_delta"])
                audit.emit(
                    action="survival.nested_model_comparison",
                    job_id=job_id,
                    fields={
                        "m2_vs_m1_justified": m2m1["justifies_complex"],
                        "m3_vs_m2_justified": m3m2["justifies_complex"],
                    },
                )
            except Exception as exc:  # noqa: BLE001
                survival_summary["model_comparison_skipped"] = str(exc)

            # ---- v0.3: causal mediation analysis (TP53 -> HRD -> survival)
            try:
                from tp53_hrd.mediation import mediation_to_dict, run_mediation
                med = run_mediation(scored, n_bootstrap=1000, seed=42)
                survival_summary["mediation_tp53_via_hrd"] = mediation_to_dict(med)
                metrics["mediation_indirect_log_hr"] = float(med.indirect_log_hr)
                if med.proportion_mediated is not None:
                    metrics["mediation_proportion_mediated"] = float(med.proportion_mediated)
                audit.emit(
                    action="survival.mediation.computed",
                    job_id=job_id,
                    fields={
                        "indirect_log_hr": float(med.indirect_log_hr),
                        "proportion_mediated": med.proportion_mediated,
                        "ci_low": med.indirect_ci_low,
                        "ci_high": med.indirect_ci_high,
                    },
                )
            except Exception as exc:  # noqa: BLE001
                survival_summary["mediation_skipped"] = str(exc)

        # ---- v0.3: clonal hierarchy approximation from MAF VAF ranking
        try:
            from tp53_hrd.clonal import cohort_clonal_calls, cohort_clonal_summary
            calls_df = cohort_clonal_calls(cohort, maf)
            clonal_summary = cohort_clonal_summary(calls_df)
            (out_dir / "clonal_calls.tsv").write_text(
                calls_df.to_csv(sep="\t", index=True)
            )
            survival_summary["clonal_hierarchy_summary"] = clonal_summary
            if clonal_summary.get("founder_rate_among_tp53_mut") is not None:
                metrics["clonal_founder_rate"] = float(clonal_summary["founder_rate_among_tp53_mut"])
            audit.emit(
                action="clonal.hierarchy.computed",
                job_id=job_id,
                fields=clonal_summary,
            )
        except Exception as exc:  # noqa: BLE001
            survival_summary["clonal_skipped"] = str(exc)

        survival_path = out_dir / "survival_summary.json"
        with survival_path.open("w", encoding="utf-8") as fh:
            json.dump(survival_summary, fh, indent=2, default=str)

        # 7. KM plot
        plot_path = out_dir / "km-severity-bands.png"
        make_km_plot(scored, plot_path)

        elapsed_ms = (time.time() - t0) * 1000.0
        metrics["body_elapsed_ms"] = elapsed_ms

        # --- end body ----------------------------------------------------

        tracking.log_metrics(metrics)

        summary = {
            "cohort": {
                "n_total": cohort_summary.n_total,
                "n_mut": cohort_summary.n_mut,
                "n_wt": cohort_summary.n_wt,
                "tier_counts": cohort_summary.tier_counts,
            },
            "severity_band_counts": band_counts,
            "logrank_3band_p": logrank_3band_p,
            "logrank_2arm_p": two_arm["logrank_p"],
            "cox_hr": cox["hr"],
            "cox_p": cox["p_value"],
            "cox_concordance": cox["concordance"],
            "artifacts": {
                "per_patient_json": str(results_path),
                "survival_summary_json": str(survival_path),
                "km_plot_png": str(plot_path),
            },
        }

    audit.emit(
        action="pipeline_end",
        job_id=job_id,
        fields={"metrics": metrics, "summary": summary},
    )

    return {
        "job_id": job_id,
        "metrics": metrics,
        "summary": summary,
    }


@click.group()
def cli() -> None:
    """tp53_hrd demonstration pipeline."""


@cli.command()
@click.option(
    "--manifest",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=Path("data/manifest.yaml"),
)
@click.option(
    "--out",
    type=click.Path(file_okay=False, path_type=Path),
    default=Path("data"),
)
@click.option(
    "--write-checksums",
    is_flag=True,
    help="Write freshly-computed sha256 values back into the manifest.",
)
def fetch(manifest: Path, out: Path, write_checksums: bool) -> None:
    """Download + checksum-verify the public inputs declared in the manifest."""
    result = fetch_manifest(manifest, out)
    if write_checksums:
        result["checksums_written"] = write_manifest_checksums(manifest, result)
    click.echo(json.dumps(result, indent=2))


@cli.command(name="refresh-manifest")
@click.option(
    "--manifest",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=Path("data/manifest.yaml"),
)
def refresh_manifest_cmd(manifest: Path) -> None:
    """Re-enumerate the GDC TCGA-LAML MAF set into the manifest's inputs block.

    Leaves sha256 blank; follow with ``fetch --write-checksums`` to fill them.
    """
    n = refresh_manifest(manifest)
    click.echo(json.dumps({"manifest": str(manifest), "inputs_written": n}, indent=2))


@cli.command()
@click.option("--name", default="demo")
@click.option(
    "--out",
    type=click.Path(file_okay=False, path_type=Path),
    default=Path("artifacts"),
)
@click.option("--maf-dir", type=click.Path(path_type=Path), default=MAF_DIR)
@click.option("--clinical-json", type=click.Path(path_type=Path), default=CLINICAL_JSON)
@click.option("--seed", type=int, default=DEFAULT_SEED)
def run(name: str, out: Path, maf_dir: Path, clinical_json: Path, seed: int) -> None:
    """Run the end-to-end TP53-HRD severity pipeline."""
    result = run_pipeline(name, out, maf_dir=maf_dir, clinical_json=clinical_json, seed=seed)
    click.echo(json.dumps(result, indent=2, default=str))


if __name__ == "__main__":
    cli()
