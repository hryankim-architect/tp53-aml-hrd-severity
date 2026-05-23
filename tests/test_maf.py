"""Tests for MAF loading utilities."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from tp53_hrd.maf import (
    EXPECTED_COLUMNS,
    load_aliquot_maf,
    load_combined_maf,
    patient_from_barcode,
)


class TestPatientFromBarcode:
    def test_full_tcga_barcode(self):
        assert patient_from_barcode("TCGA-AB-2935-03A-01W-0745-08") == "TCGA-AB-2935"

    def test_short_tcga_barcode(self):
        assert patient_from_barcode("TCGA-AB-2935") == "TCGA-AB-2935"

    def test_two_aliquots_same_patient(self):
        # The Saturday-discovered duplicate
        a = patient_from_barcode("TCGA-AB-2935-03A-01W-0745-08")
        b = patient_from_barcode("TCGA-AB-2935-03A-01W-0755-09")
        assert a == b == "TCGA-AB-2935"

    def test_non_tcga_returns_unchanged(self):
        assert patient_from_barcode("UUID-1234") == "UUID-1234"

    def test_handles_none(self):
        assert patient_from_barcode(None) == ""

    def test_handles_non_string(self):
        assert patient_from_barcode(42) == ""


# Minimal real-shaped MAF fixture — 6 columns of the 140 in production MAFs.
# The full schema is asserted by tests that touch real GDC data.
_FIXTURE_MAF = """\
#version gdc-1.0.0
#annotation.spec gdc-2.0.0-aliquot-merged-masked
#filedate 20220516
Hugo_Symbol\tVariant_Classification\tVariant_Type\tTumor_Sample_Barcode\tHGVSp_Short\tExon_Number\tt_depth\tt_alt_count\tcase_id\thotspot
TP53\tMissense_Mutation\tSNP\tTCGA-AB-2935-03A-01W-0745-08\tp.R248Q\t7/11\t103\t45\tcase-uuid-1\tY
TP53\tMissense_Mutation\tSNP\tTCGA-AB-2943-03A-01W-0745-08\tp.R273C\t8/11\t34\t29\tcase-uuid-2\tY
DNMT3A\tMissense_Mutation\tSNP\tTCGA-AB-2935-03A-01W-0745-08\tp.R882H\t23/23\t89\t38\tcase-uuid-1\tY
"""


@pytest.fixture
def fixture_maf_path(tmp_path: Path) -> Path:
    p = tmp_path / "fixture.maf.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(_FIXTURE_MAF)
    return p


class TestLoadAliquotMaf:
    def test_loads_three_rows(self, fixture_maf_path):
        df = load_aliquot_maf(fixture_maf_path)
        assert len(df) == 3

    def test_skips_comment_lines(self, fixture_maf_path):
        df = load_aliquot_maf(fixture_maf_path)
        # No row should start with "#" — those are comment lines
        assert not any(str(s).startswith("#") for s in df["Hugo_Symbol"])

    def test_has_expected_subset_of_columns(self, fixture_maf_path):
        df = load_aliquot_maf(fixture_maf_path)
        for c in EXPECTED_COLUMNS:
            assert c in df.columns, f"missing column: {c}"


class TestLoadCombinedMaf:
    def test_concatenates_multiple_files(self, tmp_path):
        # Create two fixtures
        d = tmp_path / "mafs" / "aliquot-1"
        d.mkdir(parents=True)
        for sub in ["a/x.maf.gz", "b/y.maf.gz"]:
            (tmp_path / "mafs" / sub).parent.mkdir(parents=True, exist_ok=True)
            with gzip.open(tmp_path / "mafs" / sub, "wt") as fh:
                fh.write(_FIXTURE_MAF)

        combined = load_combined_maf(tmp_path / "mafs")
        assert len(combined) == 6  # 3 rows × 2 files

    def test_raises_when_no_mafs(self, tmp_path):
        empty = tmp_path / "empty"
        empty.mkdir()
        with pytest.raises(FileNotFoundError, match="no .maf.gz files"):
            load_combined_maf(empty)
