"""TCGA-LAML clinical metadata fetching and normalization.

GDC returns clinical fields as a deeply nested JSON document — each case has
``demographic`` and ``diagnoses`` sub-objects with mixed presence depending on
the patient's status. This module flattens the response to a tidy DataFrame
with one row per patient and the four columns the survival analysis needs:

* ``patient_id`` — TCGA submitter ID (e.g. ``TCGA-AB-2935``).
* ``vital_status`` — ``"Alive"`` or ``"Dead"``.
* ``os_days`` — overall survival time. For dead patients this is
  ``demographic.days_to_death``; for living patients it is
  ``diagnoses.days_to_last_follow_up``. Either can be missing, in which case
  the row is dropped (no usable survival information).
* ``os_event`` — ``1`` for death, ``0`` for censored.
* ``age_at_diagnosis`` — informational, kept for sanity checks.
"""

from __future__ import annotations

import json
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

import pandas as pd

GDC_CASES_URL = "https://api.gdc.cancer.gov/cases"


def _downloads_allowed() -> bool:
    """Reach the GDC API only when AI_ALLOW_DOWNLOAD=1 (default: offline)."""
    import os
    return os.environ.get("AI_ALLOW_DOWNLOAD", "") not in ("", "0", "false", "False")


CLINICAL_FIELDS = [
    "submitter_id",
    "demographic.vital_status",
    "demographic.days_to_death",
    "diagnoses.days_to_last_follow_up",
    "diagnoses.age_at_diagnosis",
]


def fetch_tcga_laml_clinical(out_path: Path | None = None, timeout: float = 30.0) -> dict[str, Any]:
    """POST to the GDC cases endpoint for the TCGA-LAML project.

    Returns the parsed JSON. If ``out_path`` is supplied, the raw JSON is also
    written there for reproducibility (and to skip re-fetching on later runs).
    """
    if not _downloads_allowed():
        raise RuntimeError(
            "GDC clinical fetch required but downloads are disabled (offline mode): "
            f"{GDC_CASES_URL}\n  Seed once with AI_ALLOW_DOWNLOAD=1, or provide a cached "
            "clinical.json under data/tcga-laml/."
        )
    payload = {
        "filters": {
            "op": "in",
            "content": {
                "field": "cases.project.project_id",
                "value": ["TCGA-LAML"],
            },
        },
        "size": 500,
        "fields": ",".join(CLINICAL_FIELDS),
    }

    req = urllib.request.Request(
        GDC_CASES_URL,
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        raw = resp.read()

    body = json.loads(raw.decode("utf-8"))

    if out_path is not None:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_bytes(raw)

    return body


def parse_clinical(body: dict[str, Any]) -> pd.DataFrame:
    """Flatten the GDC clinical JSON into a tidy DataFrame.

    Drops patients without usable survival information (no days_to_death AND
    no days_to_last_follow_up).
    """
    hits = body.get("data", {}).get("hits", [])
    rows: list[dict[str, Any]] = []

    for case in hits:
        demographic = case.get("demographic") or {}
        # diagnoses is a list (a patient may have multiple); take the first.
        diagnoses_list = case.get("diagnoses") or []
        diagnoses = diagnoses_list[0] if diagnoses_list else {}

        vital_status = demographic.get("vital_status")
        days_to_death = demographic.get("days_to_death")
        days_to_last_follow_up = diagnoses.get("days_to_last_follow_up")

        # OS days: death wins over follow-up
        if vital_status == "Dead" and days_to_death is not None:
            os_days = days_to_death
            os_event = 1
        elif vital_status == "Alive" and days_to_last_follow_up is not None:
            os_days = days_to_last_follow_up
            os_event = 0
        else:
            # Missing survival information — patient unusable for survival.
            continue

        rows.append(
            {
                "patient_id": case.get("submitter_id"),
                "vital_status": vital_status,
                "os_days": float(os_days),
                "os_event": int(os_event),
                "age_at_diagnosis": diagnoses.get("age_at_diagnosis"),
            }
        )

    df = pd.DataFrame(rows)
    if df.empty:
        return df
    return df.sort_values("patient_id").reset_index(drop=True)


def load_or_fetch(json_path: Path) -> pd.DataFrame:
    """Convenience: read cached JSON if present, otherwise fetch + cache."""
    if json_path.exists():
        body = json.loads(json_path.read_text())
    else:
        body = fetch_tcga_laml_clinical(out_path=json_path)
    return parse_clinical(body)
