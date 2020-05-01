"""
Copyright 2020 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Trycycler

This file is part of Trycycler. Trycycler is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Trycycler is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Trycycler.
If not, see <http://www.gnu.org/licenses/>.
"""

import collections
import sys

from .alignment import align_reads_to_seq
from .intrange import IntRange
from .log import log, section_header, explanation
from .misc import get_sequence_file_type, load_fasta, iterate_fastq, get_fastq_stats
from .software import check_minimap2


def partition(args):
    welcome_message()
    check_inputs_and_requirements(args)
    best_clusters = align_reads(args.cluster_dirs, args.reads, args.threads, args.min_aligned_len,
                                args.min_read_cov)
    save_reads_per_cluster(args.cluster_dirs, args.reads, best_clusters)


def welcome_message():
    section_header('Starting Trycycler partitioning')
    explanation('Trycycler partition is a tool for partitioning reads by cluster. I.e. each read '
                'will be assigned to one cluster and saved in a file for that cluster.')


def check_inputs_and_requirements(args):
    check_input_reads(args.reads)
    check_input_clusters(args.cluster_dirs)
    check_required_software()


def check_input_reads(filename):
    read_type = get_sequence_file_type(filename)
    if read_type != 'FASTQ':
        sys.exit(f'\nError: input reads ({filename}) are not in FASTQ format')
    log(f'Input reads: {filename}')
    read_count, total_size, n50 = get_fastq_stats(filename)
    log(f'  {read_count:,} reads ({total_size:,} bp)')
    log(f'  N50 = {n50:,} bp')
    log()


def check_input_clusters(cluster_dirs):
    if len(cluster_dirs) < 1:
        sys.exit('\nError: one or more input cluster directories are required')
    log(f'Input clusters:')
    for d in cluster_dirs:
        contigs = sorted(d.glob('2_all_seqs.fasta'))
        if not contigs:
            sys.exit(f'\nError: there is not 2_all_seqs.fasta file in {d}')
        assert len(contigs) == 1
        f = contigs[0]
        contig_type = get_sequence_file_type(f)
        if contig_type != 'FASTA':
            sys.exit(f'\nError: input contig file ({f}) is not in FASTA format')
        seqs = load_fasta(f)
        if len(seqs) == 0:
            sys.exit(f'\nError: contig file ({f}) contains no sequences')
        noun = 'contig' if len(seqs) == 1 else 'contigs'
        mean_len = sum(len(s[1]) for s in seqs) // len(seqs)
        log(f'  {f}: {len(seqs)} {noun}, mean length = {mean_len:,} bp')
    log()


def check_required_software():
    log('Checking required software:')
    check_minimap2()
    log()


def align_reads(cluster_dirs, reads, threads, min_aligned_len, min_read_cov):
    section_header('Aligning reads to each contig')
    explanation('The reads are independently aligned to each of the contigs and Trycycler will '
                'remember the single best alignment for each read.')
    best_clusters = {}
    best_matching_bases = collections.defaultdict(int)
    for d in cluster_dirs:
        contigs = sorted(d.glob('2_all_seqs.fasta'))
        assert len(contigs) == 1
        seqs = load_fasta(contigs[0])

        for name, seq in seqs:
            seq_len = len(seq)
            log(f'{name} ({len(seq):,} bp)', end=': ')
            doubled_seq = seq + seq
            alignments = align_reads_to_seq(reads, doubled_seq, threads, include_cigar=False)

            # Toss out alignments entirely in the second half of the doubled sequence.
            alignments = [a for a in alignments if a.ref_start < seq_len]
            log(f'{len(alignments):,} alignments')

            # Group alignments by read name.
            alignments_by_read = collections.defaultdict(list)
            for a in alignments:
                alignments_by_read[a.query_name].append(a)

            # Only consider reads for which enough sequence aligned and a large enough fraction
            # of the read aligned. Note that it's okay for this to come in multiple alignments,
            # so a read which glitches out in the middle can still pass.
            for read_name, read_alignments in alignments_by_read.items():
                r = IntRange()
                read_length = None
                for a in read_alignments:
                    read_length = a.query_length
                    r.add_range(a.query_start, a.query_end)
                total_len = r.total_length()
                read_coverage = 100.0 * total_len / read_length
                if total_len < min_aligned_len or read_coverage < min_read_cov:
                    continue

                # If the read passed the checks, then we remember its best single alignment.
                for a in read_alignments:
                    if a.matching_bases > best_matching_bases[read_name]:
                        best_clusters[read_name] = d
                        best_matching_bases[read_name] = a.matching_bases
    log()

    return best_clusters


def save_reads_per_cluster(cluster_dirs, reads, best_clusters):
    section_header('Saving reads to cluster directories')
    explanation('Reads are now saved in a file in the cluster directory to which they aligned '
                'best.')
    for d in cluster_dirs:
        save_reads_one_cluster(d, reads, best_clusters)


def save_reads_one_cluster(cluster_dir, reads, best_clusters):
    cluster_reads = cluster_dir / '4_reads.fastq'
    log(f'{cluster_reads}:')
    total_read_count, total_read_bases = 0, 0
    cluster_read_count, cluster_read_bases = 0, 0
    with open(cluster_reads, 'wt') as f:
        for read_name, header, seq, qual in iterate_fastq(reads):
            total_read_count += 1
            total_read_bases += len(seq)
            if read_name in best_clusters:
                best_cluster = best_clusters[read_name]
                if best_cluster == cluster_dir:
                    cluster_read_count += 1
                    cluster_read_bases += len(seq)
                    f.write(f'{header}\n{seq}\n+\n{qual}\n')
    count_percent = 100 * cluster_read_count / total_read_count
    bases_percent = 100 * cluster_read_bases / total_read_bases
    log(f'  {cluster_read_count:,} reads ({count_percent:.2f}%)')
    log(f'  {cluster_read_bases:,} bases ({bases_percent:.2f}%)')
    log()
