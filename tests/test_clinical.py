"""Tests for clinical metadata parsing.

The fetch itself (HTTP POST to GDC) is not unit-tested — it's exercised by
the real ``make data`` step. The parser is the testable surface and is
covered here with synthetic GDC-shaped JSON.
"""

from __future__ import annotations

from tp53_hrd.clinical import parse_clinical


def _gdc_like_response(cases: list[dict]) -> dict:
    return {"data": {"hits": cases}}


class TestParseClinical:
    def test_dead_patient_uses_days_to_death(self):
        body = _gdc_like_response(
            [
                {
                    "submitter_id": "TCGA-AB-0001",
                    "demographic": {"vital_status": "Dead", "days_to_death": 365.0},
                    "diagnoses": [{"age_at_diagnosis": 21000}],
                }
            ]
        )
        df = parse_clinical(body)
        assert len(df) == 1
        row = df.iloc[0]
        assert row["os_days"] == 365.0
        assert row["os_event"] == 1
        assert row["vital_status"] == "Dead"

    def test_alive_patient_uses_days_to_last_follow_up(self):
        body = _gdc_like_response(
            [
                {
                    "submitter_id": "TCGA-AB-0002",
                    "demographic": {"vital_status": "Alive"},
                    "diagnoses": [
                        {"days_to_last_follow_up": 1800.0, "age_at_diagnosis": 18000}
                    ],
                }
            ]
        )
        df = parse_clinical(body)
        assert len(df) == 1
        row = df.iloc[0]
        assert row["os_days"] == 1800.0
        assert row["os_event"] == 0

    def test_drops_patient_missing_survival(self):
        body = _gdc_like_response(
            [
                {
                    "submitter_id": "TCGA-AB-0003",
                    "demographic": {"vital_status": "Alive"},
                    "diagnoses": [{"age_at_diagnosis": 22000}],
                },
                {
                    "submitter_id": "TCGA-AB-0004",
                    "demographic": {"vital_status": "Dead", "days_to_death": 100.0},
                    "diagnoses": [],
                },
            ]
        )
        df = parse_clinical(body)
        # 0003 has no follow-up days, 0004 is Dead with valid days_to_death
        assert list(df["patient_id"]) == ["TCGA-AB-0004"]

    def test_handles_empty_response(self):
        df = parse_clinical(_gdc_like_response([]))
        assert df.empty

    def test_sorted_by_patient_id(self):
        body = _gdc_like_response(
            [
                {
                    "submitter_id": "TCGA-AB-0099",
                    "demographic": {"vital_status": "Dead", "days_to_death": 100.0},
                    "diagnoses": [{}],
                },
                {
                    "submitter_id": "TCGA-AB-0001",
                    "demographic": {"vital_status": "Dead", "days_to_death": 200.0},
                    "diagnoses": [{}],
                },
            ]
        )
        df = parse_clinical(body)
        assert list(df["patient_id"]) == ["TCGA-AB-0001", "TCGA-AB-0099"]
