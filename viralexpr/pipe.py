"""
pipe.py
Main pipeline for ViralExpress: viral gene expression quantification from RNA-seq data.

Inspired by Centromics pipeline architecture:
  - Pipeline class encapsulates all steps
  - Checkpoint system for resumable execution
  - Unified logger with timestamps
  - Modular design for each analysis stage

Usage:
    viralexpr -r genome.fa -g annotation.gff -c data_dir -s sample_info.txt -t 8
    # Or after pip install:
    python -m viralexpr.pipe -r genome.fa -g annotation.gff -c data_dir -s sample_info.txt -t 8
"""

import argparse
import os
import re
import sys
from pathlib import Path

from .__version__ import version
from .MPRUN import logger, run_cmd
from .Dtools import mkdirs, check_ckp, mk_ckp
from .build_index import ensure_index
from .validate_annotation import ensure_valid_gtf


def makeArgparse():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='ViralExpress: viral gene expression quantification from infected host RNA-seq data',
        epilog='''
Example:
    viralexpr -r virus_genome.fa -g virus_annotation.gtf -c ./clean_fastq -s sample_info.txt -t 8

Input Files:
    -r: Viral genome FASTA file
    -g: Viral gene annotation (GFF3 or GTF)
    -c: Directory containing paired-end FASTQ files
    -s: Sample info file (3 columns, no header: sample_id R1 R2)

Output:
    FPKB, FPKM, and TPM normalized expression matrices
    Comparison plots (matplotlib/seaborn)
        '''
    )

    # Input
    group_in = parser.add_argument_group('Input')
    group_in.add_argument('-r', '--ref', required=True, metavar='FILE',
                          dest='genome_fasta',
                          help='Viral genome FASTA file [required]')
    group_in.add_argument('-g', '--gff', required=True, metavar='FILE',
                          dest='annotation',
                          help='Viral gene annotation (GFF3 or GTF) [required]')
    group_in.add_argument('-c', '--clean-data', required=True, metavar='DIR',
                          dest='clean_data_dir',
                          help='Directory with clean FASTQ files (*.fq.gz) [required]')
    group_in.add_argument('-s', '--sample-info', required=True, metavar='FILE',
                          dest='sample_info',
                          help='Sample info (3 columns, no header: sample_id R1 R2) [required]')

    # Output
    group_out = parser.add_argument_group('Output')
    group_out.add_argument('-o', '--outdir', default='viralexpr_results', metavar='DIR',
                           dest='outdir',
                           help='Output directory [default: %(default)s]')
    group_out.add_argument('--index-dir', default='', metavar='DIR',
                           dest='index_dir',
                           help='Directory for HISAT2 index (default: auto-generated beside FASTA)')

    # Options
    group_opt = parser.add_argument_group('Options')
    group_opt.add_argument('-t', '--threads', type=int, default=1, metavar='INT',
                           dest='threads',
                           help='Number of threads per sample [default: %(default)s]')
    group_opt.add_argument('-p', '--parallel', type=int, default=1, metavar='INT',
                           dest='parallel',
                           help='Number of samples to align in parallel [default: %(default)s]')
    group_opt.add_argument('-q', '--qc', action='store_true', default=False,
                           dest='run_qc',
                           help='Run fastp QC/trimming (default: off, input assumed clean)')
    group_opt.add_argument('--overwrite', action='store_true', default=False,
                           help='Overwrite existing results even if checkpoints exist')
    group_opt.add_argument('--cleanup', action='store_true', default=False,
                           help='Remove temporary files after completion')
    group_opt.add_argument('-v', '--version', action='version', version=version)

    args = parser.parse_args()
    return args


class Pipeline:
    def __init__(self, **kwargs):
        self.__dict__.update(**kwargs)
        # Resolve absolute paths
        self.outdir = os.path.realpath(self.outdir)
        self.clean_data_dir = os.path.realpath(self.clean_data_dir)
        self.genome_fasta = os.path.realpath(self.genome_fasta)
        self.annotation = os.path.realpath(self.annotation)
        self.sample_info = os.path.realpath(self.sample_info)

        # Set up directory structure
        self.dirs = {
            'logs': os.path.join(self.outdir, 'logs'),
            'qc': os.path.join(self.outdir, '01_qc'),
            'mapping': os.path.join(self.outdir, '02_mapping'),
            'quant': os.path.join(self.outdir, '03_quantification'),
            'fpkb': os.path.join(self.outdir, '04_fpkb'),
            'summary': os.path.join(self.outdir, '05_summary'),
        }
        for d in self.dirs.values():
            mkdirs(d)

        # Script directory
        self.bindir = os.path.dirname(os.path.realpath(__file__))

        # Set up logging to file
        log_file = os.path.join(self.dirs['logs'], 'pipeline.log')
        from logging import FileHandler
        fh = FileHandler(log_file)
        fh.setFormatter(logger.handlers[0].formatter if logger.handlers else None)
        logger.addHandler(fh)

        logger.info('Command: {}'.format(' '.join(sys.argv)))
        logger.info('Version: {}'.format(version))
        logger.info('Output directory: {}'.format(self.outdir))

    def mk_ckp(self, name):
        """Generate checkpoint file path."""
        return os.path.join(self.dirs['logs'], '{}.ok'.format(name))

    def run(self):
        """Execute the full pipeline with checkpoint support."""
        # Step 0: Build or check HISAT2 index
        self.step_index()

        # Step 0b: Validate annotation
        self.step_validate_annotation()

        # Step 1: Quality control (optional)
        self.step_qc()

        # Step 2: Alignment
        self.step_alignment()

        # Step 3: Quantification
        self.step_quantification()

        # Step 4: Compute FPKB/FPKM/TPM
        self.step_fpkb()

        # Step 5: Summary
        self.step_summary()

        # Cleanup
        if self.cleanup:
            self.step_cleanup()

        logger.info('Pipeline completed.')

    def step_index(self):
        """Step 0: Build or check HISAT2 index."""
        ckp = self.mk_ckp('index')
        if check_ckp(ckp, overwrite=self.overwrite):
            logger.info('Step 0 (Index): checkpoint exists, skipping')
            # Read index prefix from checkpoint
            with open(ckp) as f:
                self.index_prefix = f.read().strip()
            return

        logger.info('## Step 0: Building/checking HISAT2 index')
        self.index_prefix = ensure_index(
            self.genome_fasta,
            index_dir=self.index_dir if self.index_dir else None,
            threads=self.threads
        )
        with open(ckp, 'w') as f:
            f.write(self.index_prefix)
        mk_ckp(ckp)

    def step_validate_annotation(self):
        """Step 0b: Validate annotation file."""
        ckp = self.mk_ckp('annotation')
        if check_ckp(ckp, overwrite=self.overwrite):
            logger.info('Step 0b (Annotation): checkpoint exists, skipping')
            with open(ckp) as f:
                self.validated_gtf = f.read().strip()
            return

        logger.info('## Step 0b: Validating annotation')
        output_gtf = os.path.join(self.dirs['quant'], 'annotation_validated.gtf')
        self.validated_gtf = ensure_valid_gtf(self.annotation, output_gtf)
        with open(ckp, 'w') as f:
            f.write(self.validated_gtf)
        mk_ckp(ckp)

    def step_qc(self):
        """Step 1: Optional fastp QC."""
        if not self.run_qc:
            logger.info('Step 1 (QC): skipped (use -q to enable)')
            return

        ckp = self.mk_ckp('qc')
        if check_ckp(ckp, overwrite=self.overwrite):
            logger.info('Step 1 (QC): checkpoint exists, skipping')
            return

        logger.info('## Step 1: Quality Control (fastp)')
        sample_count = 0
        with open(self.sample_info) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split()
                if len(parts) != 3:
                    continue
                sid, r1, r2 = parts
                sample_count += 1
                in1 = os.path.join(self.clean_data_dir, r1)
                in2 = os.path.join(self.clean_data_dir, r2)
                out1 = os.path.join(self.dirs['qc'], '{}_R1.fq.gz'.format(sid))
                out2 = os.path.join(self.dirs['qc'], '{}_R2.fq.gz'.format(sid))

                cmd = [
                    'fastp',
                    '-i', in1, '-I', in2,
                    '-o', out1, '-O', out2,
                    '-j', os.path.join(self.dirs['qc'], '{}_fastp.json'.format(sid)),
                    '-h', os.path.join(self.dirs['qc'], '{}_fastp.html'.format(sid)),
                    '--detect_adapter_for_pe',
                    '--cut_front', '--cut_tail',
                    '--cut_window_size', '4', '--cut_mean_quality', '20',
                    '--length_required', '20',
                    '--thread', str(self.threads)
                ]
                run_cmd(cmd)

        logger.info('QC completed for {} samples'.format(sample_count))
        mk_ckp(ckp)

    def _get_fastq_dir(self):
        """Return directory containing FASTQ files (QC output or original)."""
        if self.run_qc:
            return self.dirs['qc']
        return self.clean_data_dir

    def _align_sample(self, sample):
        """Align a single sample. Called by multiprocessing pool."""
        sid, r1, r2 = sample
        input_dir = self._get_fastq_dir()

        if self.run_qc:
            in1 = os.path.join(input_dir, '{}_R1.fq.gz'.format(sid))
            in2 = os.path.join(input_dir, '{}_R2.fq.gz'.format(sid))
        else:
            in1 = os.path.join(input_dir, r1)
            in2 = os.path.join(input_dir, r2)

        bam = os.path.join(self.dirs['mapping'], '{}.bam'.format(sid))
        summary = os.path.join(self.dirs['mapping'], '{}.summary'.format(sid))

        logger.info('Aligning: {}'.format(sid))

        cmd1 = [
            'hisat2',
            '-p', str(self.threads),
            '--dta', '--dta-cufflinks',
            '-x', self.index_prefix,
            '-1', in1, '-2', in2,
            '--summary-file', summary,
            '--rg-id', sid,
            '--rg', 'SM:{}'.format(sid),
            '--rg', 'PL:ILLUMINA'
        ]
        cmd2 = [
            'samtools', 'sort',
            '--threads', str(self.threads),
            '-o', bam, '-'
        ]

        import subprocess
        p1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p2 = subprocess.Popen(cmd2, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p1.stdout.close()
        stdout, stderr = p2.communicate()
        if p2.returncode != 0:
            logger.error('Alignment failed for {}'.format(sid))
            raise RuntimeError('Alignment failed for {}'.format(sid))

        run_cmd(['samtools', 'index', bam])
        logger.info('Alignment done: {}'.format(sid))
        return sid

    def step_alignment(self):
        """Step 2: HISAT2 alignment with optional parallel processing."""
        ckp = self.mk_ckp('alignment')
        if check_ckp(ckp, overwrite=self.overwrite):
            logger.info('Step 2 (Alignment): checkpoint exists, skipping')
            return

        logger.info('## Step 2: Alignment (HISAT2)')

        samples = []
        with open(self.sample_info) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split()
                if len(parts) != 3:
                    continue
                samples.append(parts)

        logger.info('Total samples: {}, parallel: {}, threads per sample: {}'.format(
            len(samples), self.parallel, self.threads))

        if self.parallel <= 1:
            for sample in samples:
                self._align_sample(sample)
        else:
            from multiprocessing import Pool
            with Pool(processes=self.parallel) as pool:
                results = pool.map(self._align_sample, samples)
            logger.info('All {} samples aligned'.format(len(results)))

        logger.info('Alignment completed for {} samples'.format(len(samples)))
        mk_ckp(ckp)

    def step_quantification(self):
        """Step 3: featureCounts quantification."""
        ckp = self.mk_ckp('quantification')
        if check_ckp(ckp, overwrite=self.overwrite):
            logger.info('Step 3 (Quantification): checkpoint exists, skipping')
            return

        logger.info('## Step 3: Quantification (featureCounts)')

        bam_list = os.path.join(self.dirs['quant'], 'bam_list.txt')
        with open(bam_list, 'w') as fout:
            for bam in sorted(Path(self.dirs['mapping']).glob('*.bam')):
                fout.write(str(bam) + '\n')

        cmd = [
            'featureCounts',
            '-p', '-Q', '10',
            '-t', 'CDS',
            '-g', 'gene_id',
            '-s', '0',
            '-T', str(self.threads),
            '-a', self.validated_gtf,
            '-o', os.path.join(self.dirs['quant'], 'raw_counts.txt')
        ]
        with open(bam_list) as f:
            bams = [line.strip() for line in f if line.strip()]
        cmd.extend(bams)

        run_cmd(cmd)

        # Extract clean count matrix
        count_matrix = os.path.join(self.dirs['quant'], 'count_matrix.txt')
        with open(os.path.join(self.dirs['quant'], 'raw_counts.txt')) as fin, \
             open(count_matrix, 'w') as fout:
            for i, line in enumerate(fin):
                if i == 0:
                    continue
                cols = line.rstrip('\n').split('\t')
                # Columns: Geneid Chr Start End Strand Length [sample_bams...]
                out_cols = [cols[0]] + cols[6:]
                fout.write('\t'.join(out_cols) + '\n')

        # Fix header (remove .bam suffix and path)
        with open(count_matrix) as f:
            lines = f.readlines()
        header_parts = lines[0].strip().split('\t')
        clean_header = [header_parts[0]]
        for part in header_parts[1:]:
            clean_name = os.path.basename(part).replace('.bam', '')
            clean_header.append(clean_name)
        lines[0] = '\t'.join(clean_header) + '\n'
        with open(count_matrix, 'w') as f:
            f.writelines(lines)

        logger.info('Quantification completed')
        mk_ckp(ckp)

    def step_fpkb(self):
        """Step 4: Compute FPKB, FPKM, TPM."""
        ckp = self.mk_ckp('fpkb')
        if check_ckp(ckp, overwrite=self.overwrite):
            logger.info('Step 4 (FPKB): checkpoint exists, skipping')
            return

        logger.info('## Step 4: Computing FPKB, FPKM, TPM')

        # Extract gene lengths
        gene_lengths = os.path.join(self.dirs['quant'], 'gene_lengths.txt')
        with open(os.path.join(self.dirs['quant'], 'raw_counts.txt')) as fin, \
             open(gene_lengths, 'w') as fout:
            fout.write('gene_id\tlength\n')
            for i, line in enumerate(fin):
                if i <= 1:
                    continue
                cols = line.strip().split('\t')
                fout.write('{}\t{}\n'.format(cols[0], cols[5]))

        # Extract total reads per sample from HISAT2 summaries
        reads_summary = os.path.join(self.dirs['mapping'], 'reads_summary.txt')
        with open(reads_summary, 'w') as fout:
            fout.write('sample_id\ttotal_reads\n')
            with open(self.sample_info) as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split()
                    if len(parts) != 3:
                        continue
                    sid = parts[0]
                    summary_file = os.path.join(self.dirs['mapping'], '{}.summary'.format(sid))

                    total_reads = None
                    if os.path.exists(summary_file):
                        with open(summary_file) as sf:
                            first_line = sf.readline().strip()
                            if first_line.isdigit():
                                total_reads = first_line
                            else:
                                # Try second line
                                second_line = sf.readline().strip()
                                if second_line.isdigit():
                                    total_reads = second_line

                    if total_reads is None:
                        bam_file = os.path.join(self.dirs['mapping'], '{}.bam'.format(sid))
                        import subprocess
                        result = subprocess.run(
                            ['samtools', 'view', '-c', '-f', '1', bam_file],
                            capture_output=True, text=True
                        )
                        if result.returncode == 0:
                            total_reads = result.stdout.strip()
                        else:
                            result = subprocess.run(
                                ['samtools', 'view', '-c', bam_file],
                                capture_output=True, text=True
                            )
                            total_reads = result.stdout.strip()

                    fout.write('{}\t{}\n'.format(sid, total_reads))

        # Run compute_fpkb.py
        count_matrix = os.path.join(self.dirs['quant'], 'count_matrix.txt')
        fpkb_script = os.path.join(self.bindir, 'compute_fpkb.py')
        cmd = [
            'python3', fpkb_script,
            '--counts', count_matrix,
            '--lengths', gene_lengths,
            '--reads', reads_summary,
            '--output-prefix', os.path.join(self.dirs['fpkb'], 'viralexprion')
        ]
        run_cmd(cmd)

        logger.info('FPKB computation completed')
        mk_ckp(ckp)

    def step_summary(self):
        """Step 5: Generate summary report."""
        logger.info('## Step 5: Summary Report')

        summary_file = os.path.join(self.dirs['summary'], 'pipeline_summary.txt')
        with open(summary_file, 'w') as f:
            f.write('=' * 80 + '\n')
            f.write('ViralExpress v{} - Pipeline Summary\n'.format(version))
            f.write('=' * 80 + '\n')
            f.write('Run Date: {}\n'.format(__import__('datetime').datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
            f.write('\nInput:\n')
            f.write('  Reference: {}\n'.format(self.genome_fasta))
            f.write('  Annotation: {}\n'.format(self.annotation))
            f.write('  Samples: {}\n'.format(self._count_samples()))
            f.write('  Threads: {}\n'.format(self.threads))
            f.write('\nOutput:\n')
            f.write('  Directory: {}\n'.format(self.outdir))
            f.write('\nResults:\n')
            f.write('  QC reports: {}\n'.format(self.dirs['qc']))
            f.write('  Alignments: {}\n'.format(self.dirs['mapping']))
            f.write('  Raw counts: {}\n'.format(os.path.join(self.dirs['quant'], 'count_matrix.txt')))
            f.write('  FPKB matrix: {}\n'.format(os.path.join(self.dirs['fpkb'], 'viralexprion_FPKB.txt')))
            f.write('  FPKM matrix: {}\n'.format(os.path.join(self.dirs['fpkb'], 'viralexprion_FPKM.txt')))
            f.write('  TPM matrix: {}\n'.format(os.path.join(self.dirs['fpkb'], 'viralexprion_TPM.txt')))
            f.write('\nLogs:\n')
            f.write('  Pipeline log: {}\n'.format(os.path.join(self.dirs['logs'], 'pipeline.log')))
            f.write('=' * 80 + '\n')

        logger.info('Summary written to: {}'.format(summary_file))

    def step_cleanup(self):
        """Remove temporary files."""
        logger.info('Cleaning up temporary files...')
        # Remove intermediate files but keep outputs
        for tmp_file in ['bam_list.txt']:
            path = os.path.join(self.dirs['quant'], tmp_file)
            if os.path.exists(path):
                os.remove(path)

    def _count_samples(self):
        """Count number of samples in sample info file."""
        count = 0
        with open(self.sample_info) as f:
            for line in f:
                if line.strip():
                    count += 1
        return count


def main():
    args = makeArgparse()
    pipeline = Pipeline(**vars(args))
    pipeline.run()


if __name__ == '__main__':
    main()