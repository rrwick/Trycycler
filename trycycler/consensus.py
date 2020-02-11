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

import string
import sys

from .base_scores import get_per_base_scores
from .log import log, section_header, explanation
from .misc import get_sequence_file_type, load_fasta, get_fastq_stats
from .pairwise import get_all_pairwise_coordinates


def consensus(args):
    welcome_message()
    check_inputs_and_requirements(args)
    seqs, seq_names, seq_lengths = load_seqs(args.cluster_dir)
    cluster_name = args.cluster_dir.name
    pairwise_alignments = load_pairwise_alignments(args.cluster_dir, seq_names, seq_lengths)
    circular = not args.linear
    per_base_scores, plot_max = get_per_base_scores(seqs, args.cluster_dir, circular, args.threads,
                                                    args.plot_qual)
    consensus_seq = get_consensus_seq(seqs, per_base_scores, pairwise_alignments)
    save_seqs_to_fasta({cluster_name + '_consensus': consensus_seq},
                       args.cluster_dir / '5_consensus.fasta')
    get_per_base_scores({'consensus': consensus_seq}, args.cluster_dir, circular, args.threads,
                        args.plot_qual, plot_max, consensus=True)


def welcome_message():
    section_header('Starting Trycycler consensus')
    explanation('Trycycler consensus is a tool for combining multiple contigs from the same '
                'long-read set (e.g. assemblies from different assemblers) into a consensus '
                'contig that takes the best parts of each.')


def check_inputs_and_requirements(args):
    check_input_reads(args.cluster_dir)
    check_cluster_directory(args.cluster_dir)
    check_seqs(args.cluster_dir)
    check_required_software()


def load_contig_sequences(filenames):
    contig_seqs, fasta_names = {}, {}
    for i, f in enumerate(filenames):
        letter = string.ascii_uppercase[i]
        seqs = load_fasta(f)
        assert len(seqs) == 1
        contig_seqs[letter] = seqs[0][1]
        fasta_names[letter] = f
    return contig_seqs, fasta_names


def check_input_reads(cluster_dir):
    filename = cluster_dir / '4_reads.fastq'
    read_type = get_sequence_file_type(filename)
    if read_type != 'FASTQ':
        sys.exit(f'Error: input reads ({filename}) are not in FASTQ format')
    log(f'Input reads: {filename}')
    read_count, total_size, n50 = get_fastq_stats(filename)
    log(f'  {read_count:,} reads ({total_size:,} bp)')
    log(f'  N50 = {n50:,} bp')
    log()


def check_seqs(cluster_dir):
    filename = cluster_dir / '2_all_seqs.fasta'
    log(f'Input contigs: {filename}')
    contig_type = get_sequence_file_type(filename)
    if contig_type != 'FASTA':
        sys.exit(f'Error: input contig file ({filename}) is not in FASTA format')
    seqs = load_fasta(filename)
    if len(seqs) == 0:
        sys.exit(f'Error: input contig file ({filename}) contains no sequences')
    contig_names = set()
    for contig_name, seq in seqs:
        if contig_name in contig_names:
            sys.exit(f'Error: duplicate contig name: {contig_name}')
        contig_names.add(contig_name)
        log(f'  {contig_name}: {len(seq):,} bp')
    log()


def check_cluster_directory(directory):
    if directory.is_file():
        sys.exit(f'Error: output directory ({directory}) already exists as a file')
    if not directory.is_dir():
        sys.exit(f'Error: output directory ({directory}) does not exist')

    seq_file = directory / '2_all_seqs.fasta'
    if not seq_file.is_file():
        sys.exit(f'Error: output directory ({directory}) does not contain2_all_seqs.fasta')

    pairwise_file = directory / '3_pairwise_alignments'
    if not pairwise_file.is_file():
        sys.exit(f'Error: output directory ({directory}) does not contain 3_pairwise_alignments')

    reads_file = directory / '4_reads.fastq'
    if not reads_file.is_file():
        sys.exit(f'Error: output directory ({directory}) does not contain 4_reads.fastq')


def check_required_software():
    pass
    # TODO
    # TODO
    # TODO
    # TODO
    # TODO


def load_seqs(cluster_dir):
    filename = cluster_dir / '2_all_seqs.fasta'
    seqs = dict(load_fasta(filename))
    seq_names = sorted(seqs.keys())
    seq_lengths = {name: len(seq) for name, seq in seqs.items()}
    return seqs, seq_names, seq_lengths


def load_pairwise_alignments(cluster_dir, seq_names, seq_lengths):
    filename = cluster_dir / '3_pairwise_alignments'
    pairwise_cigars = {}
    with open(filename, 'rt') as f:
        for line in f:
            parts = line.strip().split('\t')
            assert len(parts) == 3
            a, b, cigar = parts
            pairwise_cigars[(a, b)] = cigar
    return get_all_pairwise_coordinates(seq_names, pairwise_cigars, seq_lengths)


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
    double_longest = 2 * max(len(seq) for seq in seqs.values())

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
        if len(consensus_seq) > double_longest:
            sys.exit('\nError: consensus sequence has grown too long - there is a cycle in the '
                     'pairwise alignments')

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
