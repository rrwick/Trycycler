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
import random
import string
import sys

from .base_scores import get_per_base_scores, add_indels_to_per_base_scores
from .log import log, section_header, explanation
from .misc import get_sequence_file_type, load_fasta, get_fastq_stats
from .software import check_minimap2


def consensus(args):
    welcome_message()
    check_inputs_and_requirements(args)

    seqs, seq_names, seq_lengths = load_seqs(args.cluster_dir)
    msa_seqs, msa_names, msa_lengths = load_msa(args.cluster_dir)
    sanity_check_msa(seqs, seq_names, seq_lengths, msa_seqs, msa_names, msa_lengths)

    circular = not args.linear
    per_base_scores, plot_max = get_per_base_scores(seqs, args.cluster_dir, circular,
                                                    args.threads, args.plot_qual)
    per_base_scores = add_indels_to_per_base_scores(msa_seqs, per_base_scores)
    consensus_seq_with_gaps, consensus_seq_without_gaps = \
        get_consensus_seq(msa_seqs, per_base_scores)
    save_seqs_to_fasta({args.cluster_dir.name + '_consensus': consensus_seq_with_gaps},
                       args.cluster_dir / '5_consensus_with_gaps.fasta')
    save_seqs_to_fasta({args.cluster_dir.name + '_consensus': consensus_seq_without_gaps},
                       args.cluster_dir / '6_consensus.fasta')
    get_per_base_scores({'consensus': consensus_seq_without_gaps}, args.cluster_dir, circular,
                        args.threads, args.plot_qual, plot_max, consensus=True)


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
        sys.exit(f'\nError: input reads ({filename}) are not in FASTQ format')
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
        sys.exit(f'\nError: input contig file ({filename}) is not in FASTA format')
    seqs = load_fasta(filename)
    if len(seqs) == 0:
        sys.exit(f'\nError: input contig file ({filename}) contains no sequences')
    contig_names = set()
    for contig_name, seq in seqs:
        if contig_name in contig_names:
            sys.exit(f'\nError: duplicate contig name: {contig_name}')
        contig_names.add(contig_name)
        log(f'  {contig_name}: {len(seq):,} bp')
    log()


def check_cluster_directory(directory):
    if directory.is_file():
        sys.exit(f'\nError: output directory ({directory}) already exists as a file')
    if not directory.is_dir():
        sys.exit(f'\nError: output directory ({directory}) does not exist')

    seq_file = directory / '2_all_seqs.fasta'
    if not seq_file.is_file():
        sys.exit(f'\nError: output directory ({directory}) does not contain2_all_seqs.fasta')

    pairwise_file = directory / '3_msa.fasta'
    if not pairwise_file.is_file():
        sys.exit(f'\nError: output directory ({directory}) does not contain 3_msa.fasta')

    reads_file = directory / '4_reads.fastq'
    if not reads_file.is_file():
        sys.exit(f'\nError: output directory ({directory}) does not contain 4_reads.fastq')


def check_required_software():
    log('Checking required software:')
    check_minimap2()
    log()


def load_seqs(cluster_dir):
    filename = cluster_dir / '2_all_seqs.fasta'
    seqs = dict(load_fasta(filename))
    seq_names = sorted(seqs.keys())
    seq_lengths = {name: len(seq) for name, seq in seqs.items()}
    return seqs, seq_names, seq_lengths


def load_msa(cluster_dir):
    filename = cluster_dir / '3_msa.fasta'
    seqs = dict(load_fasta(filename))
    seq_names = sorted(seqs.keys())
    seq_lengths = {name: len(seq) for name, seq in seqs.items()}
    return seqs, seq_names, seq_lengths


def sanity_check_msa(seqs, seq_names, seq_lengths, msa_seqs, msa_names, msa_lengths):
    assert seq_names == msa_names
    msa_length = msa_lengths[seq_names[0]]
    for n in seq_names:
        assert msa_lengths[n] == msa_length
        assert seq_lengths[n] <= msa_lengths[n]
        assert seqs[n] == msa_seqs[n].replace('-', '')


def save_seqs_to_fasta(seqs, filename):
    seq_word = 'sequence' if len(seqs) == 1 else 'sequences'
    log(f'Saving {seq_word} to file: {filename}')
    with open(filename, 'wt') as fasta:
        for name, seq in seqs.items():
            fasta.write(f'>{name}\n')
            fasta.write(f'{seq}\n')
    log()


def get_consensus_seq(msa_seqs, per_base_scores):
    section_header('Consensus sequence')
    explanation('Trycycler now builds the consensus sequence by switching between contigs '
                'to stay on whichever has the highest local score.')

    seq_names = list(msa_seqs.keys())
    msa_length = len(msa_seqs[seq_names[0]])

    # Sanity check!
    for n in seq_names:
        assert len(msa_seqs[n]) == msa_length == len(per_base_scores[n])

    counts = {n: 0 for n in seq_names}  # number of bases from each sequence in the consensus
    consensus_seq = []

    log('Consensus sequence composition:')
    for i in range(msa_length):
        best_base, best_seq_name = get_best_base(msa_seqs, per_base_scores, seq_names, i)
        consensus_seq.append(best_base)
        counts[best_seq_name] += 1
        if i % 1000 == 0:
            log_proportion(counts)
    log_proportion(counts)
    log('\n')

    # Sanity check: each base in the consensus sequence should contribute to the counts.
    assert sum(counts.values()) == len(consensus_seq)

    consensus_seq_with_gaps = ''.join(consensus_seq)
    consensus_seq_without_gaps = consensus_seq_with_gaps.replace('-', '')

    log(f'Consensus sequence length: {len(consensus_seq_without_gaps):,} bp')
    log()
    return consensus_seq_with_gaps, consensus_seq_without_gaps


def log_proportion(counts):
    total = sum(counts.values())
    proportions = []
    for seq_name, count in counts.items():
        proportion = 100.0 * count / total
        proportions.append(f'{seq_name}: {proportion:.2f}%')
    log('\r  ' + ', '.join(proportions), end='    ')


def get_best_base(msa_seqs, per_base_scores, seq_names, i):
    # Get the base with the highest total score. This will probably be the most common base at that
    # position in the MSA.
    base_scores = collections.defaultdict(int)
    for name in seq_names:
        base = msa_seqs[name][i]
        score = per_base_scores[name][i]
        base_scores[base] += score
    best_base = max(base_scores, key=base_scores.get)

    # Get the sequence name with the best base and the best score. If multiple sequences tie, then
    # we choose one at random.
    best_seq_names, best_score = [], None
    for name in seq_names:
        base = msa_seqs[name][i]
        if base == best_base:
            score = per_base_scores[name][i]
            if best_score is None or score > best_score:
                best_seq_names = [name]
                best_score = score
            elif score == best_score:
                best_seq_names.append(name)
    assert best_score is not None
    assert len(best_seq_names) >= 1
    best_seq_name = random.choice(best_seq_names)

    return best_base, best_seq_name
