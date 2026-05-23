"""End-to-end integration tests for the P3 pipeline.

The scaffold-inherited tests checked the demo body (which has been replaced
by the real TP53-HRD pipeline). These tests substitute tiny synthetic
fixtures so the pipeline can be exercised without the 1.5MB MAF tarball.

Three concerns are covered:

* :func:`test_pipeline_writes_three_artifacts` — the wired pipeline produces
  ``cohort-15-results.json`` / ``survival_summary.json`` /
  ``km-severity-bands.png`` from a fixture dataset.
* :func:`test_audit_chain_is_valid_after_pipeline` — the audit-chain
  substrate hook fires correctly when the pipeline runs.
* :func:`test_audit_chain_detects_tamper` — tampering with an audit entry
  invalidates the chain (regression guard for the substrate's tamper
  detection guarantee).
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

from tp53_hrd import audit, pipeline

# Minimal MAF schema — only the columns the pipeline actually reads.
# Two TP53-mutant patients + two non-TP53 patients so the WT pool is non-empty.
_FIXTURE_MAF = """\
#version gdc-1.0.0
#filedate 20220516
Hugo_Symbol\tVariant_Classification\tVariant_Type\tTumor_Sample_Barcode\tHGVSp_Short\tExon_Number\tt_depth\tt_alt_count\tcase_id\thotspot
TP53\tMissense_Mutation\tSNP\tTCGA-AB-0001-03A-01W-0001-08\tp.R248Q\t7/11\t100\t45\tcase-uuid-1\tY
TP53\tMissense_Mutation\tSNP\tTCGA-AB-0002-03A-01W-0001-08\tp.C176Y\t5/11\t100\t80\tcase-uuid-2\tN
DNMT3A\tMissense_Mutation\tSNP\tTCGA-AB-0003-03A-01W-0001-08\tp.R882H\t23/23\t100\t40\tcase-uuid-3\tN
DNMT3A\tMissense_Mutation\tSNP\tTCGA-AB-0004-03A-01W-0001-08\tp.R882H\t23/23\t100\t40\tcase-uuid-4\tN
"""


def _build_fixture_data(root: Path) -> None:
    """Create the smallest MAF + clinical setup the pipeline can run on."""
    maf_dir = root / "data" / "tcga-laml" / "mafs" / "aliquot-fixture"
    maf_dir.mkdir(parents=True, exist_ok=True)
    maf_path = maf_dir / "fixture.maf.gz"
    with gzip.open(maf_path, "wt") as fh:
        fh.write(_FIXTURE_MAF)

    clinical_json = root / "data" / "tcga-laml" / "clinical.json"
    clinical_json.parent.mkdir(parents=True, exist_ok=True)
    body = {
        "data": {
            "hits": [
                {
                    "submitter_id": f"TCGA-AB-{i:04d}",
                    "demographic": {
                        "vital_status": "Dead" if i % 2 == 0 else "Alive",
                        "days_to_death": 200.0 if i % 2 == 0 else None,
                    },
                    "diagnoses": [
                        {
                            "days_to_last_follow_up": None if i % 2 == 0 else 1500.0,
                            "age_at_diagnosis": 20000 + i * 100,
                        }
                    ],
                }
                for i in range(1, 5)
            ]
        }
    }
    clinical_json.write_text(json.dumps(body))


def _run_with_fixture(tmp_path: Path) -> dict:
    return pipeline.run_pipeline(
        "smoke",
        tmp_path / "artifacts",
        maf_dir=tmp_path / "data" / "tcga-laml" / "mafs",
        clinical_json=tmp_path / "data" / "tcga-laml" / "clinical.json",
    )


def test_pipeline_writes_three_artifacts(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("AUDIT_HOST", raising=False)
    monkeypatch.delenv("MLFLOW_TRACKING_URI", raising=False)
    _build_fixture_data(tmp_path)

    result = _run_with_fixture(tmp_path)

    assert "summary" in result
    artifacts = result["summary"]["artifacts"]
    for key in ("per_patient_json", "survival_summary_json", "km_plot_png"):
        path = Path(artifacts[key])
        assert path.exists(), f"missing artifact: {key}"
        assert path.stat().st_size > 0, f"empty artifact: {key}"


def test_audit_chain_is_valid_after_pipeline(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("AUDIT_HOST", raising=False)
    _build_fixture_data(tmp_path)

    _run_with_fixture(tmp_path)

    ok, n_entries, first_bad = audit.verify()
    assert ok, f"audit chain invalid at {first_bad}"
    # pipeline emits at least pipeline_start, cohort_built, pipeline_end
    assert n_entries >= 3


def test_audit_chain_detects_tamper(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("AUDIT_HOST", raising=False)
    _build_fixture_data(tmp_path)

    _run_with_fixture(tmp_path)

    ledger = audit.DEFAULT_LEDGER
    lines = ledger.read_text().splitlines()
    assert len(lines) >= 2

    tampered = json.loads(lines[0])
    tampered["fields"]["out_dir"] = "/etc/evil"
    lines[0] = json.dumps(tampered, sort_keys=True, separators=(",", ":"))
    ledger.write_text("\n".join(lines) + "\n")

    ok, _, first_bad = audit.verify()
    assert not ok
    assert first_bad is not None
