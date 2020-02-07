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

import random
import string
import sys

from .circularisation import circularise
from .initial_check import initial_sanity_check
from .log import log, section_header, explanation
from .misc import get_sequence_file_type, load_fasta, get_fastq_stats
from .pairwise import get_pairwise_alignments
from .starting_seq import get_starting_seq, rotate_to_starting_seq
from . import settings


def align(args):
    random.seed(0)
    welcome_message()
    check_inputs_and_requirements(args)
    seqs, fasta_names = load_contig_sequences(args.cluster_dir)
    initial_sanity_check(seqs, args.max_mash_dist, args.max_length_diff)
    starting_seq, seqs = get_starting_seq(seqs, args.threads)
    circular = not args.linear
    if circular:
        seqs = circularise(seqs, args.reads, args.threads)
        seqs = rotate_to_starting_seq(seqs, starting_seq)
    save_seqs_to_fasta(seqs, args.cluster_dir / '2_all_seqs.fasta')
    pairwise_alignments = get_pairwise_alignments(seqs)
    save_pairwise_alignments(pairwise_alignments, args.cluster_dir / '3_pairwise_alignments')


def welcome_message():
    section_header('Starting Trycycler align')
    explanation('Trycycler align is a tool for reconciling multiple alternative contigs with each '
                'other.')


def check_inputs_and_requirements(args):
    check_input_reads(args.reads)
    check_cluster_directory(args.cluster_dir)
    check_input_contigs(args.cluster_dir)
    check_required_software()


def load_contig_sequences(cluster_dir):
    filenames = get_contigs_from_cluster_dir(cluster_dir)
    contig_seqs, fasta_names = {}, {}
    for i, f in enumerate(filenames):
        letter = string.ascii_uppercase[i]
        seqs = load_fasta(f)
        assert len(seqs) == 1
        contig_seqs[letter] = seqs[0][1]
        fasta_names[letter] = f
    return contig_seqs, fasta_names


def check_input_reads(filename):
    read_type = get_sequence_file_type(filename)
    if read_type != 'FASTQ':
        sys.exit(f'Error: input reads ({filename}) are not in FASTQ format')
    log(f'Input reads: {filename}')
    read_count, total_size, n50 = get_fastq_stats(filename)
    log(f'  {read_count:,} reads ({total_size:,} bp)')
    log(f'  N50 = {n50:,} bp')
    log()


def check_input_contigs(cluster_dir):
    filenames = get_contigs_from_cluster_dir(cluster_dir)
    if len(filenames) < 2:
        sys.exit('Error: two or more input contigs are required')
    if len(filenames) > settings.MAX_INPUT_CONTIGS:
        sys.exit(f'Error: you cannot have more than {settings.MAX_INPUT_CONTIGS} input contigs')
    log(f'Input contigs:')
    contig_names = set()
    for i, f in enumerate(filenames):
        contig_type = get_sequence_file_type(f)
        if contig_type != 'FASTA':
            sys.exit(f'Error: input contig file ({f}) is not in FASTA format')
        seqs = load_fasta(f)
        if len(seqs) == 0:
            sys.exit(f'Error: input contig file ({f}) contains no sequences')
        if len(seqs) > 1:
            sys.exit(f'Error: input contig file ({f}) contains multiple sequences')
        contig_name = seqs[0][0]
        if contig_name in contig_names:
            sys.exit(f'Error: duplicate contig name: {contig_name}')
        contig_names.add(contig_name)
        contig_len = len(seqs[0][1])
        letter = string.ascii_uppercase[i]
        log(f'  {letter}: {f} ({contig_name}: {contig_len:,} bp)')
    log()


def get_contigs_from_cluster_dir(cluster_dir):
    contig_dir = cluster_dir / '1_contigs'
    if not contig_dir.is_dir():
        sys.exit(f'Error: contig directory ({contig_dir}) does not exist')
    return sorted(contig_dir.glob('*.fasta'))


def check_cluster_directory(directory):
    if directory.is_file():
        sys.exit(f'Error: output directory ({directory}) already exists as a file')
    if not directory.is_dir():
        sys.exit(f'Error: output directory ({directory}) does not exist')


def check_required_software():
    pass
    # TODO
    # TODO
    # TODO
    # TODO
    # TODO


def save_seqs_to_fasta(seqs, filename):
    seq_word = 'sequence' if len(seqs) == 1 else 'sequences'
    log(f'Saving {seq_word} to file: {filename}')
    with open(filename, 'wt') as fasta:
        for name, seq in seqs.items():
            fasta.write(f'>{name}\n')
            fasta.write(f'{seq}\n')
    log()


def get_consensus_seq(seqs, per_base_scores, pairwise_alignments):
    section_header('Consensus sequence')
    explanation('Trycycler now builds the consensus sequence by switching between contigs '
                'to stay on whichever has the highest local score.')

    seq_names = list(seqs.keys())
    other_seq_names = get_other_seq_names(seq_names)

    current_seq_name = choose_starting_sequence(seqs, per_base_scores)
    current_seq = seqs[current_seq_name]
    current_pos = 0

    consensus_seq = []
    counts = {n: 0 for n in seq_names}
    log('Consensus sequence composition:')
    while True:
        current_score = per_base_scores[current_seq_name][current_pos]
        best_other_score, best_other_seq_name, best_other_pos = 0, None, None
        for other_seq_name in other_seq_names[current_seq_name]:
            pairwise = pairwise_alignments[(current_seq_name, other_seq_name)]
            other_pos = pairwise[current_pos]
            if other_pos is not None:
                other_score = per_base_scores[other_seq_name][other_pos]
                if other_score > best_other_score:
                    best_other_score = other_score
                    best_other_seq_name = other_seq_name
                    best_other_pos = other_pos
        if best_other_score > current_score:
            current_seq_name = best_other_seq_name
            current_seq = seqs[current_seq_name]
            current_pos = best_other_pos
            log_proportion(counts)
        consensus_seq.append(current_seq[current_pos])
        counts[current_seq_name] += 1

        current_pos += 1
        if current_pos >= len(current_seq):
            break

    log_proportion(counts)
    log('\n')
    return ''.join(consensus_seq)


def choose_starting_sequence(seqs, per_base_scores):
    best_seq_name, best_score = None, 0
    for seq_name in seqs.keys():
        starting_score = per_base_scores[seq_name][0]
        if best_seq_name is None or starting_score > best_score:
            best_seq_name = seq_name
            best_score = starting_score
    return best_seq_name


def get_other_seq_names(seq_names):
    """
    Builds a dictionary where each seq name gives a list of the other seq names.
    """
    other_seq_names = {}
    for seq_name in seq_names:
        other_seq_names[seq_name] = [n for n in seq_names if n != seq_name]
    return other_seq_names


def log_proportion(counts):
    total = sum(counts.values())
    proportions = []
    for seq_name, count in counts.items():
        proportion = 100.0 * count / total
        proportions.append(f'{seq_name}: {proportion:.2f}%')
    log('\r  ' + ', '.join(proportions), end='    ')


def save_pairwise_alignments(pairwise_alignments, filename):
    log(f'Saving pairwise alignments to file: {filename}')
    with open(filename, 'wt') as f:
        for seq_names, cigar in pairwise_alignments.items():
            a, b = seq_names
            f.write(f'{a}\t{b}\t')
            f.write(cigar)
            f.write('\n')
    log()
