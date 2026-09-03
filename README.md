# ViralExpr

ViralExpr is a universal framework for viral transcript quantification from infected plant RNA-seq data. By providing a viral genome FASTA and an annotation file (GFF3/GTF) as input, it automatically completes alignment, quantification, and normalization from RNA-seq data of infected tissues, outputting three expression matrices: FPKB, FPKM, and TPM.

## Features

- **Virus-agnostic**: Analyze different viruses by simply replacing the genome and annotation files; no code modification required
- **Index reuse**: HISAT2 index is automatically built on first run and directly reused afterwards; custom shared directories are supported
- **Annotation validation**: Automatically validates GFF3/GTF format, checks CDS presence and gene identifiers; GFF3 is automatically converted to GTF
- **Resume capability**: Checkpoints are generated at each step; re-running after interruption automatically skips completed steps
- **One-command workflow**: A single command completes the entire pipeline from raw data to expression matrices

## Dependencies

All dependencies are managed uniformly through conda; no manual installation of individual packages is required.

| Software | Purpose |
|------|------|
| bash | Workflow control |
| hisat2 | Sequence alignment |
| samtools | BAM processing |
| featureCounts | Gene quantification |
| fastp | Quality control (optional) |
| python3 + pandas/numpy/matplotlib/seaborn/scipy | Scripts and visualization |

## Installation

### Prerequisites

- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Mamba](https://mamba.readthedocs.io/en/latest/installation.html)
- bash >= 4.0

### Install (to control environment location)

```bash
# Create environment (use -p /path/to/env to change installation location)
conda env create -f viralexpr.yaml
# Or faster: mamba env create -f viralexpr.yaml

# Activate environment
conda activate viralexpr

# Install ViralExpress
pip install -e .
```

### Verify Installation

After installation, the `viralexpr` command is available. Activate the environment before each use:

```bash
conda activate viralexpr
viralexpr -h
```

## Quick Start

### Input Files

**Sample Information Table** `sample_info.txt` (no header, tab-delimited, 3 columns):

```
sample1	Sample1_R1.fq.gz	Sample1_R2.fq.gz
sample2	Sample2_R1.fq.gz	Sample2_R2.fq.gz
```

**Plotting Metadata** `metadata.txt` (with header, tab-delimited):

```
sample_id	timepoint	replicate
sample1	5d	rep1
sample2	11d	rep1
```

### Run

```bash
viralexpr \
    -r ${ViralExpr_PATH}/SRBSDV/genome.fa \
    -g ${ViralExpr_PATH}/SRBSDV/annotation.gtf \
    -c ./clean_fastq \
    -s sample_info.txt \
    -t 8 \
    -o ./results
```

### Visualization

```bash
python viralexpr/plot_comparison.py \
    --fpkb results/04_fpkb/viralexprion_FPKB.txt \
    --fpkm results/04_fpkb/viralexprion_FPKM.txt \
    --tpm results/04_fpkb/viralexprion_TPM.txt \
    --counts results/03_quantification/count_matrix.txt \
    --metadata metadata.txt \
    --output results/05_summary/plots
```

## Parameters

### viralexpr

| Parameter | Required | Description |
|------|------|------|
| `-r` | Yes | Viral genome FASTA |
| `-g` | Yes | Viral gene annotation (GFF3 or GTF) |
| `-c` | Yes | Clean FASTQ directory |
| `-s` | Yes | Sample information table |
| `-t` | Yes | Number of threads |
| `-o` | No | Output directory (default: `viralexpr_results`) |
| `--index-dir` | No | HISAT2 index directory (default: same directory as FASTA) |
| `-q` | No | Enable fastp QC/adapter trimming (default: off) |
| `--overwrite` | No | Force overwrite existing checkpoints |
| `--cleanup` | No | Clean up temporary files after completion |

### Index Sharing Example

```bash
viralexpr \
    -r ${ViralExpr_PATH}/SRBSDV/genome.fa \
    -g ${ViralExpr_PATH}/SRBSDV/annotation.gtf \
    --index-dir /shared/indices/virus_hisat2 \
    -c ./data -s sample_info.txt -t 8
```

The index is automatically built and stored in `--index-dir` on the first run; all subsequent runs pointing to the same `--index-dir` directly reuse it.

### Resume from Interruption

A `.ok` checkpoint is generated in `logs/` after each step completes. Re-running the same command after interruption automatically skips completed steps:

```bash
# First run interrupted (e.g., power outage)
viralexpr -r genome.fa -g anno.gtf -c data -s info.txt -t 8

# After power is restored, re-run the same command
viralexpr -r genome.fa -g anno.gtf -c data -s info.txt -t 8
# -> Automatically recognizes existing checkpoints and resumes from the breakpoint
```

To force re-running a specific step, add `--overwrite`:

```bash
viralexpr -r genome.fa -g anno.gtf -c data -s info.txt -t 8 --overwrite
```

## Output

```
results/
├── logs/
│   ├── pipeline.log              # Run log
│   ├── index.ok                  # Index checkpoint
│   ├── annotation.ok             # Annotation checkpoint
│   ├── qc.ok                     # QC checkpoint (if enabled)
│   ├── alignment.ok              # Alignment checkpoint
│   ├── quantification.ok         # Quantification checkpoint
│   └── fpkb.ok                   # Normalization checkpoint
├── 01_qc/                        # fastp reports (when -q is enabled)
├── 02_mapping/                   # BAM + index + HISAT2 summary
├── 03_quantification/            # featureCounts output + cleaned matrix + validated GTF
├── 04_fpkb/                      # FPKB / FPKM / TPM matrices + statistics summary
└── 05_summary/
    └── pipeline_summary.txt      # Run summary
```

## Annotation Format Requirements

**GFF3**: Standard 9 columns, containing `CDS` features; `ID=` or `Name=` attributes are used for gene identification.

**GTF**: Standard 9 columns, containing `CDS` features; `gene_id` attribute is used for gene identification.

When GFF3 is provided as input, it is automatically converted to GTF and stored in `03_quantification/annotation_validated.gtf`. If validation fails, the process exits and reports the specific error line number in the log.

## FPKB Formula

```
FPKB =(CDS mapped fragments ×10^12)/(total sequenced fragments × CDS length (bp) )
```
"CDS mapped fragments" is the number of paired-end fragments mapped to the CDS of a given gene, where one fragment corresponds to one read pair. "Total sequenced fragments" is the total number of fragments in the sequencing library. "CDS length" is in base pairs. The scaling factor 10^12 standardizes the output by kilobase of CDS length and billions of sequenced (not mapped) fragments.

## FAQ

| Problem | Cause | Solution |
|------|------|------|
| HISAT2 index build failure | Insufficient memory | Reduce `-t` thread count |
| Annotation validation failure | Missing CDS or gene_id / ID attributes | Check the error line number in `logs/pipeline.log` |
| Step skipped but results corrupted | Previous abnormal interruption | Add `--overwrite` to re-run that step |
| FASTQ files not found | Sample info table filenames do not match actual files | Verify columns 2 and 3 of `sample_info.txt` match filenames in the `-c` directory |

## Citation

If you use ViralExpress, please cite:

## License

MIT License
