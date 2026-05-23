"""MAF loading utilities for TCGA-LAML aliquot ensemble masked MAFs.

The GDC project-level MAF was deprecated in favor of per-aliquot files. A single
download from the GDC bulk endpoint produces a tarball with one MAF per aliquot::

    data/tcga-laml/mafs/<file_id>/<aliquot_uuid>.wxs.aliquot_ensemble_masked.maf.gz

Each MAF file shares the same 140-column schema (see ``EXPECTED_COLUMNS`` below)
and is preceded by 4-5 ``#``-prefixed comment lines (version, contigs, sort
order, filedate, annotation spec).

This module provides:

* :func:`load_aliquot_maf` — read one MAF.gz file into a pandas DataFrame.
* :func:`load_combined_maf` — walk a directory and concatenate all aliquot MAFs.
* :func:`patient_from_barcode` — collapse TCGA aliquot barcodes to patient IDs
  (TCGA-AB-2935-03A-01W-0745-08 → TCGA-AB-2935). Patient-level grouping is
  required for survival analysis to avoid double-counting technical replicates.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

# Columns the rest of the package depends on. Other columns are kept but
# unused by v0.1. This list is asserted in tests to catch upstream schema drift.
EXPECTED_COLUMNS = (
    "Hugo_Symbol",
    "Variant_Classification",
    "Variant_Type",
    "Tumor_Sample_Barcode",
    "HGVSp_Short",
    "Exon_Number",
    "t_depth",
    "t_alt_count",
    "case_id",
    "hotspot",
)


def load_aliquot_maf(path: Path) -> pd.DataFrame:
    """Read one aliquot MAF (.maf.gz) into a DataFrame.

    The ``#`` comment lines at the top of GDC MAFs are skipped automatically by
    ``pd.read_csv(comment="#")``. The next line is the column header.
    """
    return pd.read_csv(
        path,
        sep="\t",
        comment="#",
        compression="gzip",
        low_memory=False,
        dtype={"Chromosome": str},  # avoid mixed-int issues on chrX/chrY/chrM
    )


def load_combined_maf(maf_dir: Path) -> pd.DataFrame:
    """Concatenate every aliquot MAF in ``maf_dir`` (recursively)."""
    maf_dir = Path(maf_dir)
    paths = sorted(maf_dir.rglob("*.maf.gz"))
    if not paths:
        raise FileNotFoundError(f"no .maf.gz files found under {maf_dir}")

    frames = []
    for p in paths:
        df = load_aliquot_maf(p)
        if df.empty:
            # Aliquot with zero somatic calls — skip silently. Many TCGA-LAML
            # aliquots are sparse (median ~10 variants per sample).
            continue
        frames.append(df)

    combined = pd.concat(frames, ignore_index=True)

    missing = [c for c in EXPECTED_COLUMNS if c not in combined.columns]
    if missing:
        raise ValueError(
            f"MAF schema drift detected — missing expected columns: {missing}"
        )

    return combined


def patient_from_barcode(barcode: str) -> str:
    """Collapse a TCGA aliquot barcode to its patient ID.

    Example::

        >>> patient_from_barcode("TCGA-AB-2935-03A-01W-0745-08")
        'TCGA-AB-2935'

    Returns the input unchanged if it does not look like a TCGA barcode (e.g.
    a UUID or non-TCGA sample identifier).
    """
    if not isinstance(barcode, str):
        return ""
    parts = barcode.split("-")
    if len(parts) >= 3 and parts[0] == "TCGA":
        return "-".join(parts[:3])
    return barcode
