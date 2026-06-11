"""NIH GDC open-tier discovery for TCGA-LAML public inputs.

This module enumerates the exact, deterministic file set the demo depends on so
``data/manifest.yaml`` can be (re)generated rather than hand-maintained:

* :func:`list_laml_maf_files` — every open-tier WXS "aliquot ensemble masked"
  somatic MAF for the TCGA-LAML project, sorted by ``file_id`` for stable output.
* :func:`build_maf_inputs` — turn that listing into manifest ``inputs`` entries.
* :data:`CLINICAL_QUERY` — the canonical GDC ``/cases`` query for survival fields.

All endpoints are open tier (no token). GDC stores each data file as an immutable
artifact, so a file's bytes — and therefore its sha256 — are stable across
downloads and machines. The MAF *set* can change when GDC reprocesses a project;
that is why regenerating the manifest is an explicit, reviewable step (see
``pipeline.refresh_manifest``) rather than something ``make data`` does silently.
"""

from __future__ import annotations

import json
import urllib.request
from typing import Any

GDC_FILES_URL = "https://api.gdc.cancer.gov/files"
GDC_CASES_URL = "https://api.gdc.cancer.gov/cases"
GDC_DATA_URL = "https://api.gdc.cancer.gov/data"


def _downloads_allowed() -> bool:
    """Reach the GDC API only when AI_ALLOW_DOWNLOAD=1 (default: offline)."""
    import os
    return os.environ.get("AI_ALLOW_DOWNLOAD", "") not in ("", "0", "false", "False")


_OFFLINE_MSG = (
    "GDC access required but downloads are disabled (offline mode): {url}\n"
    "  Seed once with AI_ALLOW_DOWNLOAD=1 (e.g. _offline/bin/seed-assets.sh "
    "tp53-aml-hrd-severity), or place cached inputs under data/tcga-laml/."
)

# Deterministic filter for the TCGA-LAML open-tier somatic MAF set.
MAF_FILTERS: dict[str, Any] = {
    "op": "and",
    "content": [
        {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-LAML"]}},
        {"op": "in", "content": {"field": "files.data_format", "value": ["MAF"]}},
        {"op": "in", "content": {"field": "files.access", "value": ["open"]}},
        {"op": "in", "content": {"field": "files.experimental_strategy", "value": ["WXS"]}},
    ],
}

# Canonical clinical query. Mirrors clinical.CLINICAL_FIELDS; kept here so the
# manifest can record provenance without importing pandas.
CLINICAL_PROJECT = "TCGA-LAML"
CLINICAL_FIELDS = [
    "submitter_id",
    "demographic.vital_status",
    "demographic.days_to_death",
    "diagnoses.days_to_last_follow_up",
    "diagnoses.age_at_diagnosis",
]
CLINICAL_SIZE = 500


def _post(url: str, payload: dict[str, Any], timeout: float = 60.0) -> dict[str, Any]:
    if not _downloads_allowed():
        raise RuntimeError(_OFFLINE_MSG.format(url=url))
    req = urllib.request.Request(
        url,
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode("utf-8"))


def data_url(file_id: str) -> str:
    """Open-tier download URL for a single GDC file."""
    return f"{GDC_DATA_URL}/{file_id}"


def list_laml_maf_files(timeout: float = 60.0) -> list[dict[str, Any]]:
    """Return every TCGA-LAML open WXS masked MAF, sorted by ``file_id``.

    Each item carries ``file_id``, ``file_name``, ``file_size`` and ``md5sum``
    (the GDC-recorded md5, used as an independent integrity check on download).
    Pagination is followed to completion; the sort makes the output stable so a
    regenerated manifest diffs cleanly.
    """
    hits: list[dict[str, Any]] = []
    frm = 0
    page = 100
    while True:
        payload = {
            "filters": MAF_FILTERS,
            "size": page,
            "from": frm,
            "sort": "file_id:asc",
            "fields": "file_id,file_name,file_size,md5sum",
        }
        body = _post(GDC_FILES_URL, payload, timeout=timeout)
        data = body.get("data", {})
        batch = data.get("hits", [])
        hits.extend(batch)
        total = data.get("pagination", {}).get("total", len(hits))
        frm += len(batch)
        if not batch or frm >= total:
            break
    return sorted(hits, key=lambda h: h["file_id"])


def build_maf_inputs(files: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Convert a GDC file listing into manifest ``inputs`` entries.

    The ``sha256`` is left empty here on purpose: it is the hash of the
    *on-disk* file, filled by ``fetch --write-checksums`` after a verified
    download. ``size_bytes`` comes from GDC metadata and is informational.
    """
    inputs: list[dict[str, Any]] = []
    for f in files:
        fid = f["file_id"]
        inputs.append(
            {
                "url": data_url(fid),
                "path": f"tcga-laml/mafs/{fid}.maf.gz",
                "sha256": "",
                "size_bytes": int(f.get("file_size", 0)),
                "license": "NIH GDC Data Use Agreement (open tier)",
                "source": (
                    f"TCGA-LAML WXS aliquot ensemble masked MAF; "
                    f"GDC file {fid} ({f.get('file_name', '')})"
                ),
            }
        )
    return inputs


def clinical_query() -> dict[str, Any]:
    """The canonical GDC ``/cases`` POST body for TCGA-LAML survival fields."""
    return {
        "filters": {
            "op": "in",
            "content": {"field": "cases.project.project_id", "value": [CLINICAL_PROJECT]},
        },
        "size": CLINICAL_SIZE,
        "fields": ",".join(CLINICAL_FIELDS),
    }


def fetch_clinical_raw(timeout: float = 60.0) -> bytes:
    """Fetch the raw clinical JSON bytes (byte-stable for a fixed query)."""
    if not _downloads_allowed():
        raise RuntimeError(_OFFLINE_MSG.format(url=GDC_CASES_URL))
    req = urllib.request.Request(
        GDC_CASES_URL,
        data=json.dumps(clinical_query()).encode("utf-8"),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()
