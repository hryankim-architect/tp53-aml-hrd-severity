"""Unit tests for src/tp53_hrd/scar.py (HRD genomic-scar score, Telli 2016)."""

from __future__ import annotations

import pandas as pd

from tp53_hrd import scar


def _seg(chrom, start, end, major=1, minor=1) -> dict:
    """Synthetic segment row helper."""
    return {
        "Chromosome": chrom,
        "Start": int(start),
        "End": int(end),
        "Copy_Number": int(major + minor),
        "Major_Copy_Number": int(major),
        "Minor_Copy_Number": int(minor),
        "length_bp": int(end - start + 1),
    }


def _segments(rows: list[dict]) -> pd.DataFrame:
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# LOH
# ---------------------------------------------------------------------------
def test_loh_requires_minor_zero():
    # 20 Mb segment with minor=1 -> NOT LOH
    segs = _segments([_seg("chr1", 1, 20_000_000, major=2, minor=1)])
    assert scar.count_loh(segs) == 0


def test_loh_requires_length_above_15mb():
    # 10 Mb segment with minor=0 -> NOT LOH (too short)
    segs = _segments([_seg("chr1", 1, 10_000_000, major=2, minor=0)])
    assert scar.count_loh(segs) == 0
    # 20 Mb segment with minor=0 -> LOH
    segs = _segments([_seg("chr1", 1, 20_000_000, major=2, minor=0)])
    assert scar.count_loh(segs) == 1


def test_loh_excludes_whole_chromosome():
    # Segment covering essentially all of chr22 with minor=0 -> NOT LOH
    chr22_len = scar.CHROM_SIZES_HG38["chr22"]
    segs = _segments([_seg("chr22", 1, chr22_len, major=2, minor=0)])
    assert scar.count_loh(segs) == 0


def test_loh_excludes_sex_chromosomes():
    segs = _segments([_seg("chrX", 1, 20_000_000, major=2, minor=0)])
    assert scar.count_loh(segs) == 0


# ---------------------------------------------------------------------------
# TAI
# ---------------------------------------------------------------------------
def test_tai_requires_allelic_imbalance():
    # AI requires major != minor. Balanced 2:2 segment -> NOT TAI
    segs = _segments([_seg("chr1", 1, 30_000_000, major=2, minor=2)])
    assert scar.count_tai(segs) == 0


def test_tai_requires_telomere_proximity():
    # 30 Mb mid-chromosome AI -> NOT TAI (doesn't touch telomere)
    segs = _segments([_seg("chr1", 60_000_000, 90_000_000, major=2, minor=1)])
    assert scar.count_tai(segs) == 0
    # Same shape but anchored at chr start (p-telomere) -> TAI
    segs = _segments([_seg("chr1", 1, 30_000_000, major=2, minor=1)])
    assert scar.count_tai(segs) == 1


def test_tai_rejects_segments_spanning_centromere():
    # Spans chr1 centromere (~123 Mb) start->end touches p-telomere but
    # crosses centromere -> NOT TAI per Telli/Birkbak
    segs = _segments([_seg("chr1", 1, 200_000_000, major=2, minor=1)])
    assert scar.count_tai(segs) == 0


# ---------------------------------------------------------------------------
# LST
# ---------------------------------------------------------------------------
def test_lst_two_large_adjacent_segments_with_different_state():
    # Two 20 Mb adjacent segments on chr1 p-arm with different allele states
    # -> 1 LST transition
    segs = _segments([
        _seg("chr1", 1, 20_000_000, major=2, minor=2),
        _seg("chr1", 20_000_001, 40_000_000, major=3, minor=1),
    ])
    assert scar.count_lst(segs) == 1


def test_lst_zero_when_both_segments_short():
    # Two 5 Mb segments -> NOT LST (each <10 Mb)
    segs = _segments([
        _seg("chr1", 1, 5_000_000, major=2, minor=2),
        _seg("chr1", 5_000_001, 10_000_000, major=3, minor=1),
    ])
    assert scar.count_lst(segs) == 0


def test_lst_zero_when_states_are_identical():
    # Two 20 Mb segments with the same allele state -> NOT LST
    # (smoothing merges them into one)
    segs = _segments([
        _seg("chr1", 1, 20_000_000, major=2, minor=2),
        _seg("chr1", 20_000_001, 40_000_000, major=2, minor=2),
    ])
    assert scar.count_lst(segs) == 0


# ---------------------------------------------------------------------------
# Aggregate HRD score
# ---------------------------------------------------------------------------
def test_compute_hrd_score_known_synthetic_example():
    """Hand-built segment list with known LOH=1, TAI=1, LST=1 -> hrd_score=3."""
    segs = _segments([
        # LOH event on chr1 p-arm (start before centromere): 20 Mb, minor=0
        _seg("chr1", 1, 20_000_000, major=2, minor=0),
        # TAI: starts at chr2 telomere (Start <= 3 Mb), 30 Mb,
        # imbalanced 3:1, does NOT cross centromere (~93 Mb)
        _seg("chr2", 1, 30_000_000, major=3, minor=1),
        # LST: 2 adjacent 20 Mb segments on chr3 with different states
        _seg("chr3", 50_000_000, 70_000_000, major=2, minor=2),
        _seg("chr3", 70_000_001, 90_000_000, major=3, minor=1),
    ])
    result = scar.compute_hrd_score(segs, "TEST-001")
    assert result.loh == 1
    assert result.tai == 1
    assert result.lst == 1
    assert result.hrd_score == 3
    assert result.hrd_positive is False  # below threshold of 42


def test_compute_hrd_score_threshold_call():
    # Manually engineer 42 LOH events on different chromosomes to flip
    # hrd_positive.
    rows = []
    chroms = [f"chr{i}" for i in range(1, 23)]
    # 2 LOH segments per chromosome (p arm + q arm) on first 21 chromosomes
    for c in chroms[:21]:
        # p-arm LOH ~20 Mb
        rows.append(_seg(c, 1_000_000, 21_000_000, major=2, minor=0))
        # q-arm LOH ~20 Mb (well away from centromere)
        cen = scar.CENTROMERE_MID_HG38[c]
        rows.append(_seg(c, cen + 5_000_000, cen + 25_000_000, major=2, minor=0))
    segs = _segments(rows)
    result = scar.compute_hrd_score(segs, "TEST-HRD-POS")
    assert result.loh == 42
    assert result.hrd_score == 42
    assert result.hrd_positive is True  # exactly at threshold


def test_cohort_hrd_scores_skips_missing_files(tmp_path):
    out = scar.cohort_hrd_scores({"P1": tmp_path / "missing.txt"})
    assert out.empty
