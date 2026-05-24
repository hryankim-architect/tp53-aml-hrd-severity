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
    """Hand-built segment list with known LOH=1, TAI=1, LST=1 -> hrd_score=3.

    Anchored carefully so each segment matches *only* its intended scar
    category (e.g. the LOH segment is moved away from the p-telomere
    margin so it does not double-count as a TAI event).
    """
    segs = _segments([
        # LOH only: chr1 p-arm, away from p-telomere (start > 3 Mb margin),
        # away from centromere, 20 Mb, minor=0
        _seg("chr1", 4_000_000, 24_000_000, major=2, minor=0),
        # TAI only: chr2 p-telomere-anchored (start=1), 30 Mb, AI 3:1,
        # does NOT cross chr2 centromere (~93.9 Mb), minor != 0 so NOT LOH
        _seg("chr2", 1, 30_000_000, major=3, minor=1),
        # LST: 2 adjacent 20 Mb segments on chr3 with different states,
        # well inside the chromosome so neither end touches a telomere
        # and neither segment is LOH (minor >= 1) nor TAI (not telomere-adjacent)
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
    """Engineer exactly 42 LOH events with NO TAI / LST contamination.

    Each segment is placed away from both telomeres (start > 3 Mb,
    end < chrom_len - 3 Mb) so it cannot count as TAI. p-arm segments
    are placed entirely before the centromere for metacentric chroms
    or span the centromere for acrocentric chroms (chr13-15, 21, 22) —
    in either case NOT contributing to TAI (TAI requires extending to
    a telomere). Each pair lives on a different arm so they cannot
    form an LST transition with each other.
    """
    rows = []
    chroms = [f"chr{i}" for i in range(1, 22)]  # 21 autosomes (chr22 has very short q-arm)
    for c in chroms:
        cen = scar.CENTROMERE_MID_HG38[c]
        chrom_len = scar.CHROM_SIZES_HG38[c]
        # p-arm LOH ~20 Mb, away from p-telomere (start >> 3 Mb margin)
        rows.append(_seg(c, 4_000_000, 24_000_000, major=2, minor=0))
        # q-arm LOH ~20 Mb, away from centromere AND away from q-telomere
        q_start = cen + 5_000_000
        q_end = min(q_start + 20_000_000, chrom_len - 5_000_000)
        rows.append(_seg(c, q_start, q_end, major=2, minor=0))
    segs = _segments(rows)
    result = scar.compute_hrd_score(segs, "TEST-HRD-POS")
    assert result.loh == 42
    assert result.tai == 0, f"unexpected TAI contamination: {result.tai}"
    assert result.lst == 0, f"unexpected LST contamination: {result.lst}"
    assert result.hrd_score == 42
    assert result.hrd_positive is True  # exactly at threshold


def test_cohort_hrd_scores_skips_missing_files(tmp_path):
    out = scar.cohort_hrd_scores({"P1": tmp_path / "missing.txt"})
    assert out.empty
