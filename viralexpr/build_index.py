"""
HISAT2 index manager: check existence, auto-build if missing, reuse if present.
"""

import os
import sys
from pathlib import Path

from .MPRUN import logger, run_cmd
from .Dtools import mkdirs


def get_index_prefix(fasta_path, index_dir=None):
    """Determine index prefix path from FASTA and optional index directory."""
    fasta_path = Path(fasta_path)
    if index_dir:
        index_path = Path(index_dir)
    else:
        index_path = fasta_path.parent / f"{fasta_path.stem}_hisat2_index"
    mkdirs(index_path)
    return index_path / fasta_path.stem


def index_exists(prefix):
    """Check if all expected HISAT2 index files exist."""
    prefix = Path(prefix)
    expected = [prefix.with_name(f"{prefix.name}.{i}.ht2") for i in range(1, 9)]
    expected_alt = [prefix.with_name(f"{prefix.name}.{i}.ht2l") for i in range(1, 9)]
    return all(p.exists() for p in expected) or all(p.exists() for p in expected_alt)


def build_index(fasta_path, index_prefix, threads=1):
    """Run hisat2-build to create the index."""
    cmd = [
        "hisat2-build",
        "-p", str(threads),
        str(fasta_path),
        str(index_prefix)
    ]
    logger.info('Building HISAT2 index: {}'.format(index_prefix))
    run_cmd(cmd)
    logger.info('Index build completed')


def ensure_index(fasta_path, index_dir=None, threads=1):
    """Ensure index exists, build if necessary. Returns index prefix path."""
    fasta_path = Path(fasta_path)
    if not fasta_path.exists():
        logger.error('FASTA file not found: {}'.format(fasta_path))
        sys.exit(1)

    index_prefix = get_index_prefix(fasta_path, index_dir)

    if index_exists(index_prefix):
        logger.info('HISAT2 index already exists: {}'.format(index_prefix))
    else:
        build_index(fasta_path, index_prefix, threads)
        if not index_exists(index_prefix):
            logger.error('Index build reported success but files not found')
            sys.exit(1)

    return str(index_prefix)
