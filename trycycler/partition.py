"""
Copyright 2019 Ryan Wick (rrwick@gmail.com)
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
from .consensus import check_input_reads
from .log import log, section_header, explanation
from .misc import get_sequence_file_type, load_fasta, iterate_fastq


def partition(args):
    welcome_message()
    check_inputs_and_requirements(args)


def welcome_message():
    section_header('Starting Trycycler partitioning')
    explanation('Trycycler partition is a tool for partitioning reads by cluster. I.e. each read '
                'will be assigned to the cluster it best aligns to.')


def check_inputs_and_requirements(args):
    check_input_reads(args.reads)
    check_input_clusters(args.cluster_dirs)
    check_required_software()
    best_clusters = align_reads(args.cluster_dirs, args.reads, args.threads, args.coverage)
    save_reads_per_cluster(args.cluster_dirs, args.reads, best_clusters)


def check_input_clusters(cluster_dirs):
    if len(cluster_dirs) < 1:
        sys.exit('Error: one or more input cluster directories are required')
    log(f'Input clusters:')
    for d in cluster_dirs:
        log(f'  {d}:')
        contigs = sorted(d.glob("*.fasta"))
        if not contigs:
            sys.exit(f'Error: there are no FASTA files in {d}')
        for f in contigs:
            contig_type = get_sequence_file_type(f)
            if contig_type != 'FASTA':
                sys.exit(f'Error: input contig file ({f}) is not in FASTA format')
            seqs = load_fasta(f)
            if len(seqs) == 0:
                sys.exit(f'Error: contig file ({f}) contains no sequences')
            if len(seqs) > 1:
                sys.exit(f'Error: contig file ({f}) contains multiple sequences')
            contig_len = len(seqs[0][1])
            log(f'    {f.name} ({contig_len:,} bp)')
    log()


def check_required_software():
    pass
    # TODO
    # TODO
    # TODO
    # TODO
    # TODO


def align_reads(cluster_dirs, reads, threads, coverage_threshold):
    section_header('Aligning reads to each contig')
    explanation('The reads are independently aligned to each of the contigs. Alignments are only '
                'kept if they cover a sufficient fraction of the read, which helps to eliminate '
                'some chimeric reads.')
    best_clusters = {}
    best_matching_bases = collections.defaultdict(int)
    for d in cluster_dirs:
        contigs = sorted(d.glob("*.fasta"))
        for f in contigs:
            seqs = load_fasta(f)
            seq = seqs[0][1]
            seq_len = len(seq)
            log(f'{f} ({len(seq):,} bp)', end=': ')
            doubled_seq = seq + seq
            alignments = align_reads_to_seq(reads, doubled_seq, threads, include_cigar=False)
            alignments = [a for a in alignments if a.ref_start < seq_len]
            alignments = [a for a in alignments if a.query_cov >= coverage_threshold]
            log(f'{len(alignments):,} alignments')

            for a in alignments:
                read_name = a.query_name
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
    cluster_reads = cluster_dir / 'reads.fastq'
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
