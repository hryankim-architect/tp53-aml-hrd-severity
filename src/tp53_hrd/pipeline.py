"""End-to-end TP53-HRD severity pipeline for TCGA-LAML.

The shape mirrors the scaffold template so the audit / tracking / canary
substrate hooks fire in the same shape across every capability-portrait
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
import time
import urllib.request
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import click
import yaml

from tp53_hrd import audit, tracking
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


def _run_id(name: str) -> str:
    stamp = datetime.now(UTC).strftime("%Y%m%d-%H%M%S")
    return f"{name}-{stamp}"


def _checksum(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch_manifest(manifest_path: Path, out_dir: Path) -> dict[str, Any]:
    """Download every entry in the manifest; verify SHA-256 checksums."""
    out_dir.mkdir(parents=True, exist_ok=True)
    with manifest_path.open("r", encoding="utf-8") as fh:
        manifest = yaml.safe_load(fh) or {}

    results: list[dict[str, Any]] = []
    for entry in manifest.get("inputs", []):
        url = entry["url"]
        rel = entry["path"]
        expected = entry.get("sha256")
        size_mb = entry.get("size_mb")
        dest = out_dir / rel
        dest.parent.mkdir(parents=True, exist_ok=True)

        if dest.exists() and expected and _checksum(dest) == expected:
            results.append({"path": str(dest), "status": "cached"})
            continue

        urllib.request.urlretrieve(url, dest)
        actual = _checksum(dest)
        if expected and actual != expected:
            results.append({
                "path": str(dest),
                "status": "checksum_mismatch",
                "expected": expected,
                "actual": actual,
            })
            continue
        results.append({
            "path": str(dest),
            "status": "downloaded",
            "sha256": actual,
            "size_mb": size_mb,
        })

    return {"inputs": results}


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
            ascat_dir = out_dir / "ascat_segments"
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
    """tp53_hrd capability-portrait pipeline."""


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
def fetch(manifest: Path, out: Path) -> None:
    """Download public inputs declared in the manifest."""
    result = fetch_manifest(manifest, out)
    click.echo(json.dumps(result, indent=2))


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
