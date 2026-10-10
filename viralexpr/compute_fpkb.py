"""
compute_fpkb.py
Compute FPKB, FPKM, and TPM normalization from featureCounts output.

Denominators (v2026-10-09 fix):
  FPKB: total sequenced fragments (library size; fixed per sequencing run)
        -> count * 1e12 / (CDS length_bp * total sequenced fragments)
  FPKM: mapped fragments, taken as the sum of CDS-assigned fragments per
        sample (i.e., the column sum of the count matrix). Because the
        reference is the viral genome only, this sum is the total number of
        viral reads and grows as the virus replicates.
        -> count * 1e9 / (CDS length_bp * mapped fragments)
  TPM:  sum of length-normalized rates across genes
        -> (count / length_kb) / sum(count / length_kb) * 1e6

Note: FPKB and FPKM intentionally use DIFFERENT denominators. FPKB is scaled
by the sequencing library (fixed), FPKM by the viral mapped reads (dynamic);
the contrast between the two is the object of the accompanying manuscript.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Compute FPKB, FPKM, and TPM")
    parser.add_argument("--counts", "-c", required=True, help="Count matrix (genes x samples)")
    parser.add_argument("--lengths", "-l", required=True, help="Gene lengths (gene_id, length_bp)")
    parser.add_argument("--reads", "-r", required=True, help="Reads summary (sample_id, total_reads)")
    parser.add_argument("--output-prefix", "-o", required=True, help="Output prefix")
    return parser.parse_args()


def compute_fpkb(counts, lengths, total_reads):
    """FPKB = count / (length_kb * reads_billion) = count * 1e12 / (length_bp * total_reads)

    Denominator: total sequenced fragments (library size, fixed per run).
    """
    length_kb = lengths / 1000.0
    reads_billion = total_reads / 1e9
    fpkb = counts.div(length_kb, axis=0).div(reads_billion, axis=1)
    return fpkb.round(4)


def compute_fpkm(counts, lengths):
    """FPKM = count / (length_kb * mapped_million)

    Denominator: mapped fragments, defined as the sum of CDS-assigned
    fragments in each sample (column sum of the count matrix). Since the
    reference is the viral genome alone, this sum equals the total number of
    viral reads and therefore grows with viral replication -- which is
    precisely the behaviour described in the manuscript.
    """
    length_kb = lengths / 1000.0
    mapped = counts.sum(axis=0).astype(float)
    # Guard against samples with zero assigned fragments (no viral signal):
    # such samples yield NaN -> reported as 0.0 instead of crashing.
    mapped = mapped.replace(0, np.nan)
    mapped_million = mapped / 1e6
    fpkm = counts.div(length_kb, axis=0).div(mapped_million, axis=1)
    return fpkm.fillna(0.0).round(4)


def compute_tpm(counts, lengths):
    """TPM = (count/length_kb) / sum(count/length_kb) * 1e6"""
    length_kb = lengths / 1000.0
    rate = counts.div(length_kb, axis=0)
    tpm = rate.div(rate.sum(axis=0), axis=1) * 1e6
    return tpm.round(4)


def main():
    args = parse_args()

    counts = pd.read_csv(args.counts, sep="\t", index_col=0)
    lengths = pd.read_csv(args.lengths, sep="\t", skiprows=1, header=None, names=["gene_id", "length_bp"])
    lengths = lengths.drop_duplicates("gene_id").set_index("gene_id")["length_bp"].astype(float)
    reads = pd.read_csv(args.reads, sep="\t", skiprows=1, header=None, names=["sample_id", "total_reads"])
    reads = reads.drop_duplicates("sample_id").set_index("sample_id")["total_reads"].astype(float)

    # Validate inputs
    if counts.empty:
        raise ValueError(f"Count matrix is empty: {args.counts}")
    if counts.columns.empty:
        raise ValueError(
            "Count matrix has no sample columns. "
            f"First row of {args.counts} looks like: {counts.index[0] if len(counts.index) > 0 else 'N/A'}. "
            "Check that the file has a proper header row with sample names."
        )

    common_samples = set(counts.columns) & set(reads.index)
    if not common_samples:
        raise ValueError(
            f"No common samples between count matrix and reads summary.\n"
            f"  Samples in count matrix columns: {list(counts.columns)}\n"
            f"  Samples in reads summary:        {list(reads.index)}\n"
            "Most likely cause: count_matrix.txt is missing a header row. "
            "Delete the quantification checkpoint and re-run the pipeline."
        )

    counts = counts[list(common_samples)]
    reads = reads[list(common_samples)]
    common_genes = counts.index.intersection(lengths.index)
    counts = counts.loc[common_genes]
    lengths = lengths.loc[common_genes]

    if counts.empty or lengths.empty:
        raise ValueError(
            f"No common genes between count matrix ({len(counts.index)} genes) "
            f"and gene lengths ({len(lengths.index)} genes)."
        )

    out_prefix = Path(args.output_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    # FPKB: library-size denominator (total sequenced fragments)
    compute_fpkb(counts, lengths, reads).to_csv(f"{args.output_prefix}_FPKB.txt", sep="\t")
    # FPKM: mapped-read denominator (sum of CDS-assigned fragments per sample)
    compute_fpkm(counts, lengths).to_csv(f"{args.output_prefix}_FPKM.txt", sep="\t")
    # TPM: sum of length-normalized rates
    compute_tpm(counts, lengths).to_csv(f"{args.output_prefix}_TPM.txt", sep="\t")

    print("Done: {}".format(out_prefix))


if __name__ == "__main__":
    main()
