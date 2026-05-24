"""Fetch GDC TCGA-LAML ASCAT2 allele-specific copy-number segment files.

Used by the Arm-3 HRD scar score (see `scar.py`). Open-tier, so no auth
needed. Per-patient lookup goes through the GDC files endpoint, filters
to ASCAT2 (preferred) or ASCAT3 (fallback) allele-specific segments,
picks the latest workflow per patient, downloads, and caches on disk.

GDC schema reference: <https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/CNV_Pipeline/>
"""

from __future__ import annotations

import json
import urllib.parse
import urllib.request
from pathlib import Path

GDC_FILES = "https://api.gdc.cancer.gov/files"
GDC_DATA = "https://api.gdc.cancer.gov/data"


def _gdc_post(url: str, payload: dict) -> dict:
    req = urllib.request.Request(
        url,
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json"},
    )
    with urllib.request.urlopen(req, timeout=30) as fh:
        return json.loads(fh.read().decode("utf-8"))


def find_ascat_file_for_patient(submitter_id: str) -> dict | None:
    """Query GDC for an open-tier ASCAT allele-specific segment file.

    Prefers ASCAT2 (more samples available); falls back to ASCAT3 if
    ASCAT2 is missing for the patient.
    """
    query = {
        "filters": {
            "op": "and",
            "content": [
                {"op": "in", "content": {
                    "field": "cases.submitter_id", "value": [submitter_id]
                }},
                {"op": "in", "content": {
                    "field": "files.data_type",
                    "value": ["Allele-specific Copy Number Segment"],
                }},
                {"op": "in", "content": {
                    "field": "files.access", "value": ["open"]
                }},
            ],
        },
        "fields": "file_id,file_name,analysis.workflow_type,cases.submitter_id",
        "format": "JSON",
        "size": "10",
    }
    response = _gdc_post(GDC_FILES, query)
    hits = response.get("data", {}).get("hits", [])
    if not hits:
        return None
    # Prefer ASCAT2
    for h in hits:
        if h.get("analysis", {}).get("workflow_type") == "ASCAT2":
            return h
    # Fall back to ASCAT3
    for h in hits:
        if h.get("analysis", {}).get("workflow_type") == "ASCAT3":
            return h
    # Otherwise just return the first
    return hits[0]


def download_file(file_id: str, dest: Path) -> Path:
    """Download a GDC file (open tier; no token needed)."""
    url = f"{GDC_DATA}/{urllib.parse.quote(file_id)}"
    dest.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(url, dest)
    return dest


def fetch_cohort_ascat(submitter_ids: list[str], out_dir: Path) -> dict[str, Path]:
    """Fetch ASCAT segments for every patient in the cohort.

    Caches on disk; re-running on top of an already-populated out_dir is
    a no-op for patients whose segment file already exists.

    Returns:
        dict patient_id -> local path. Patients with no ASCAT file
        available are omitted from the result.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    result: dict[str, Path] = {}
    for pid in submitter_ids:
        # Cached?
        cached = sorted(out_dir.glob(f"{pid}__*.seg.txt"))
        if cached:
            result[pid] = cached[0]
            continue
        hit = find_ascat_file_for_patient(pid)
        if hit is None:
            continue
        out_path = out_dir / f"{pid}__{hit['file_name']}"
        try:
            download_file(hit["file_id"], out_path)
        except (OSError, urllib.error.URLError):
            continue
        result[pid] = out_path
    return result
