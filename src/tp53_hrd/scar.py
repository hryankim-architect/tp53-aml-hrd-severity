"""HRD genomic-scar score per Telli et al. 2016 (PMC6773427).

This module implements the Myriad MyChoice HRD score:

    HRD_score = LOH_count + TAI_count + LST_count

where the three components are defined as follows:

* **LOH (Loss of Heterozygosity)** — Abkevich et al. 2012, Br J Cancer:
  count of genomic segments where Minor_Copy_Number == 0, segment length
  is greater than 15 Mb, and the segment does NOT span the entire
  chromosome (a whole-chromosome LOH event reflects whole-arm loss /
  uniparental disomy, not HR repair failure).

* **TAI (Telomeric Allelic Imbalance)** — Birkbak et al. 2012, Cancer Discov:
  count of subchromosomal regions with allelic imbalance
  (Major_Copy_Number != Minor_Copy_Number) that extend to a telomere AND
  do not cross the centromere of the chromosome.

* **LST (Large-scale State Transitions)** — Popova et al. 2012, Cancer Res:
  number of chromosomal breaks between adjacent segments that are both
  at least 10 Mb long. Adjacent here means consecutive on the same
  chromosome arm; we smooth across gaps < 3 Mb before counting.

Cutoff per Telli 2016 Triple-Negative Breast Cancer validation:

    HRD-positive if HRD_score >= 42

We carry this cutoff verbatim in the implementation but flag that the
TCGA-LAML cohort was not part of the original Telli validation; per-cohort
thresholds in AML are not yet established. The README documents this
honestly.

Inputs:
    GDC TCGA-LAML ASCAT2 / ASCAT3 allele-specific segment files
    (open-tier, columns: GDC_Aliquot / Chromosome / Start / End /
    Copy_Number / Major_Copy_Number / Minor_Copy_Number).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

# ---------------------------------------------------------------------------
# Telli 2016 constants
# ---------------------------------------------------------------------------
LOH_MIN_LENGTH_BP: int = 15_000_000   # 15 Mb (Abkevich 2012)
LST_MIN_LENGTH_BP: int = 10_000_000   # 10 Mb (Popova 2012)
LST_SMOOTHING_GAP_BP: int = 3_000_000  # smooth across <3 Mb gaps (Popova 2012)
TELOMERE_MARGIN_BP: int = 3_000_000    # within 3 Mb of telomere = telomeric (Birkbak 2012)
HRD_POSITIVE_THRESHOLD: int = 42       # Telli 2016 TNBC cutoff

# ---------------------------------------------------------------------------
# GRCh38 chromosome lengths (UCSC hg38, primary assembly contigs)
# ---------------------------------------------------------------------------
CHROM_SIZES_HG38: dict[str, int] = {
    "chr1": 248_956_422, "chr2": 242_193_529, "chr3": 198_295_559,
    "chr4": 190_214_555, "chr5": 181_538_259, "chr6": 170_805_979,
    "chr7": 159_345_973, "chr8": 145_138_636, "chr9": 138_394_717,
    "chr10": 133_797_422, "chr11": 135_086_622, "chr12": 133_275_309,
    "chr13": 114_364_328, "chr14": 107_043_718, "chr15": 101_991_189,
    "chr16": 90_338_345, "chr17": 83_257_441, "chr18": 80_373_285,
    "chr19": 58_617_616, "chr20": 64_444_167, "chr21": 46_709_983,
    "chr22": 50_818_468, "chrX": 156_040_895, "chrY": 57_227_415,
}

# Centromere midpoints (GRCh38, UCSC gap/centromere track), used to split
# chromosomes into p-arm vs q-arm so TAI's "does not cross centromere" rule
# and LST's per-arm scoring are correct.
CENTROMERE_MID_HG38: dict[str, int] = {
    "chr1": 123_400_000, "chr2": 93_900_000, "chr3": 90_900_000,
    "chr4": 50_000_000, "chr5": 48_750_000, "chr6": 59_550_000,
    "chr7": 60_100_000, "chr8": 45_200_000, "chr9": 43_000_000,
    "chr10": 39_800_000, "chr11": 53_400_000, "chr12": 35_500_000,
    "chr13": 17_700_000, "chr14": 17_150_000, "chr15": 19_000_000,
    "chr16": 36_850_000, "chr17": 25_050_000, "chr18": 18_450_000,
    "chr19": 26_150_000, "chr20": 28_050_000, "chr21": 12_000_000,
    "chr22": 15_550_000, "chrX": 61_000_000, "chrY": 10_400_000,
}

# Telli's HRD count is autosome-only; sex chromosomes are excluded because
# normal allelic imbalance there is sex-dependent and confounds the metric.
AUTOSOMES: tuple[str, ...] = tuple(f"chr{i}" for i in range(1, 23))


# ---------------------------------------------------------------------------
# Loader
# ---------------------------------------------------------------------------
def load_ascat_segments(path: Path) -> pd.DataFrame:
    """Load a GDC ASCAT2/3 allele-specific segment TSV.

    Expected columns (GDC schema):
        GDC_Aliquot / Chromosome / Start / End / Copy_Number /
        Major_Copy_Number / Minor_Copy_Number

    Returns a DataFrame with normalised column names + an added
    ``length_bp`` column.
    """
    df = pd.read_csv(path, sep="\t")
    expected = {
        "Chromosome", "Start", "End", "Copy_Number",
        "Major_Copy_Number", "Minor_Copy_Number",
    }
    missing = expected - set(df.columns)
    if missing:
        raise ValueError(f"ASCAT file {path} missing columns: {sorted(missing)}")
    df = df.copy()
    df["length_bp"] = df["End"] - df["Start"] + 1
    return df


# ---------------------------------------------------------------------------
# LOH (Loss of Heterozygosity)
# ---------------------------------------------------------------------------
def is_loh_segment(row: pd.Series, chrom_sizes: dict[str, int] | None = None) -> bool:
    """Telli 2016 LOH segment: minor=0, length>15 Mb, NOT whole chromosome."""
    chrom_sizes = chrom_sizes or CHROM_SIZES_HG38
    if int(row["Minor_Copy_Number"]) != 0:
        return False
    length = int(row["length_bp"])
    if length <= LOH_MIN_LENGTH_BP:
        return False
    chrom = str(row["Chromosome"])
    if chrom not in AUTOSOMES:
        return False
    # Reject whole-chromosome events: segment covers >= 95% of chromosome length
    chrom_len = chrom_sizes.get(chrom, 0)
    return not (chrom_len and length >= 0.95 * chrom_len)


def count_loh(segments: pd.DataFrame, chrom_sizes: dict[str, int] | None = None) -> int:
    """Telli 2016 LOH count for a single sample."""
    mask = segments.apply(lambda r: is_loh_segment(r, chrom_sizes), axis=1)
    return int(mask.sum())


# ---------------------------------------------------------------------------
# TAI (Telomeric Allelic Imbalance)
# ---------------------------------------------------------------------------
def _extends_to_telomere(
    start: int, end: int, chrom: str, chrom_sizes: dict[str, int]
) -> bool:
    """True if segment touches (within margin of) either telomere."""
    margin = TELOMERE_MARGIN_BP
    if start <= margin:
        return True
    chrom_len = chrom_sizes.get(chrom, 0)
    return bool(chrom_len and end >= chrom_len - margin)


def _crosses_centromere(
    start: int, end: int, chrom: str, centromere_mid: dict[str, int]
) -> bool:
    """True if segment spans the centromere of its chromosome."""
    cen = centromere_mid.get(chrom)
    if cen is None:
        return False
    return start < cen < end


def is_tai_segment(
    row: pd.Series,
    chrom_sizes: dict[str, int] | None = None,
    centromere_mid: dict[str, int] | None = None,
) -> bool:
    """Telli 2016 TAI: AI (major != minor), extends to telomere, NOT crossing centromere."""
    chrom_sizes = chrom_sizes or CHROM_SIZES_HG38
    centromere_mid = centromere_mid or CENTROMERE_MID_HG38
    chrom = str(row["Chromosome"])
    if chrom not in AUTOSOMES:
        return False
    if int(row["Major_Copy_Number"]) == int(row["Minor_Copy_Number"]):
        return False
    start, end = int(row["Start"]), int(row["End"])
    if _crosses_centromere(start, end, chrom, centromere_mid):
        return False
    return _extends_to_telomere(start, end, chrom, chrom_sizes)


def count_tai(
    segments: pd.DataFrame,
    chrom_sizes: dict[str, int] | None = None,
    centromere_mid: dict[str, int] | None = None,
) -> int:
    """Telli 2016 TAI count for a single sample."""
    mask = segments.apply(
        lambda r: is_tai_segment(r, chrom_sizes, centromere_mid), axis=1
    )
    return int(mask.sum())


# ---------------------------------------------------------------------------
# LST (Large-scale State Transitions)
# ---------------------------------------------------------------------------
def _segment_arm(row: pd.Series, centromere_mid: dict[str, int]) -> str:
    """Return 'p' if segment is entirely before centromere, 'q' if after, else 'spanning'."""
    chrom = str(row["Chromosome"])
    cen = centromere_mid.get(chrom)
    if cen is None:
        return "unknown"
    start, end = int(row["Start"]), int(row["End"])
    if end <= cen:
        return "p"
    if start >= cen:
        return "q"
    return "spanning"


def count_lst(
    segments: pd.DataFrame,
    centromere_mid: dict[str, int] | None = None,
) -> int:
    """Popova 2012 / Telli 2016 LST count for a single sample.

    Walk each chromosome arm in genomic order. After smoothing across
    gaps < 3 Mb, count each transition between two consecutive segments
    where both are >= 10 Mb long.
    """
    centromere_mid = centromere_mid or CENTROMERE_MID_HG38
    autosome_segs = segments[segments["Chromosome"].isin(AUTOSOMES)].copy()
    if autosome_segs.empty:
        return 0
    autosome_segs["arm"] = autosome_segs.apply(
        lambda r: _segment_arm(r, centromere_mid), axis=1
    )
    autosome_segs = autosome_segs[autosome_segs["arm"].isin(("p", "q"))]

    lst_count = 0
    for (_chrom, _arm), arm_segs in autosome_segs.groupby(["Chromosome", "arm"]):
        arm_segs = arm_segs.sort_values("Start").reset_index(drop=True)
        # Smoothing: drop tiny gaps (<3 Mb) by merging adjacent segments
        # that share the same allele state.
        smoothed = _smooth_segments(arm_segs)
        # Count transitions where both adjacent segments are >= 10 Mb.
        for i in range(len(smoothed) - 1):
            left = smoothed[i]
            right = smoothed[i + 1]
            if (
                left["length_bp"] >= LST_MIN_LENGTH_BP
                and right["length_bp"] >= LST_MIN_LENGTH_BP
                and (
                    left["Major_Copy_Number"] != right["Major_Copy_Number"]
                    or left["Minor_Copy_Number"] != right["Minor_Copy_Number"]
                )
            ):
                lst_count += 1
    return lst_count


def _smooth_segments(arm_segs: pd.DataFrame) -> list[dict]:
    """Merge consecutive segments separated by <3 Mb of un-segmented sequence
    AND sharing the same (major, minor) copy state."""
    out: list[dict] = []
    for _, row in arm_segs.iterrows():
        new = {
            "Start": int(row["Start"]),
            "End": int(row["End"]),
            "length_bp": int(row["length_bp"]),
            "Major_Copy_Number": int(row["Major_Copy_Number"]),
            "Minor_Copy_Number": int(row["Minor_Copy_Number"]),
        }
        if not out:
            out.append(new)
            continue
        prev = out[-1]
        gap = new["Start"] - prev["End"]
        same_state = (
            prev["Major_Copy_Number"] == new["Major_Copy_Number"]
            and prev["Minor_Copy_Number"] == new["Minor_Copy_Number"]
        )
        if gap < LST_SMOOTHING_GAP_BP and same_state:
            # Merge into the prior segment
            prev["End"] = new["End"]
            prev["length_bp"] = prev["End"] - prev["Start"] + 1
        else:
            out.append(new)
    return out


# ---------------------------------------------------------------------------
# HRD score (Telli 2016)
# ---------------------------------------------------------------------------
@dataclass
class HRDScarResult:
    patient_id: str
    loh: int
    tai: int
    lst: int
    hrd_score: int
    hrd_positive: bool
    n_segments_input: int


def compute_hrd_score(segments: pd.DataFrame, patient_id: str) -> HRDScarResult:
    """Top-level: compute Telli 2016 HRD score for one sample."""
    loh = count_loh(segments)
    tai = count_tai(segments)
    lst = count_lst(segments)
    hrd_score = loh + tai + lst
    return HRDScarResult(
        patient_id=patient_id,
        loh=loh,
        tai=tai,
        lst=lst,
        hrd_score=hrd_score,
        hrd_positive=hrd_score >= HRD_POSITIVE_THRESHOLD,
        n_segments_input=int(len(segments)),
    )


def cohort_hrd_scores(
    patient_segment_paths: dict[str, Path],
) -> pd.DataFrame:
    """Compute HRD scores for a whole cohort.

    Args:
        patient_segment_paths: dict patient_id -> ASCAT segment file path

    Returns:
        DataFrame indexed by patient_id with columns
        [loh, tai, lst, hrd_score, hrd_positive, n_segments_input].
    """
    rows = []
    for pid, path in patient_segment_paths.items():
        try:
            segs = load_ascat_segments(Path(path))
        except (FileNotFoundError, ValueError):
            continue
        result = compute_hrd_score(segs, pid)
        rows.append({
            "patient_id": result.patient_id,
            "loh": result.loh,
            "tai": result.tai,
            "lst": result.lst,
            "hrd_score": result.hrd_score,
            "hrd_positive": result.hrd_positive,
            "n_segments_input": result.n_segments_input,
        })
    if not rows:
        return pd.DataFrame(
            columns=["loh", "tai", "lst", "hrd_score", "hrd_positive", "n_segments_input"]
        )
    return pd.DataFrame(rows).set_index("patient_id")
