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

import random
import sys

from .circularisation import circularise
from .initial_check import initial_check
from .log import log, section_header, explanation, dim, red, quit_with_error
from .misc import get_sequence_file_type, load_fasta, check_input_reads
from .pairwise import get_pairwise_alignments
from .software import check_minimap2
from .starting_seq import normalise_strands, get_starting_seq, rotate_to_starting_seq
from . import settings


def reconcile(args):
    random.seed(0)
    welcome_message()
    check_inputs_and_requirements(args)
    seqs, fasta_names = load_contig_sequences(args.cluster_dir)
    initial_check(seqs, args.max_mash_dist, args.max_length_diff)
    seqs = normalise_strands(seqs)
    circular = not args.linear
    if circular:
        seqs = circularise(seqs, args)
        seqs, starting_seq = get_starting_seq(seqs, args.threads)
        seqs = rotate_to_starting_seq(seqs, starting_seq)
    pairwise_cigars, percent_identities, worst_1kbp_identities = get_pairwise_alignments(seqs)
    print_identity_matrix(seqs, percent_identities, args.min_identity)
    print_worst_1kbp_matrix(seqs, worst_1kbp_identities, args.min_1kbp_identity)
    finished_message()
    save_seqs_to_fasta(seqs, args.cluster_dir / '2_all_seqs.fasta')


def welcome_message():
    section_header('Starting Trycycler reconcile')
    explanation('Trycycler reconcile is a tool for reconciling multiple alternative contigs with '
                'each other.')


def finished_message():
    section_header('Finished!')
    explanation('All contig sequences are now reconciled and ready for the next step in the '
                'pipeline: trycycler msa.')


def check_inputs_and_requirements(args):
    check_input_reads(args.reads, file_size_only=True)
    check_cluster_directory(args.cluster_dir)
    check_input_contigs(args.cluster_dir)
    check_required_software()


def load_contig_sequences(cluster_dir):
    filenames = get_contigs_from_cluster_dir(cluster_dir)
    contig_seqs, fasta_names = {}, {}
    for f in filenames:
        seqs = load_fasta(f)
        assert len(seqs) == 1
        name, seq = seqs[0]
        contig_seqs[name] = seq
        fasta_names[name] = f
    return contig_seqs, fasta_names


def check_input_contigs(cluster_dir, require_two_or_more=True):
    filenames = get_contigs_from_cluster_dir(cluster_dir)
    if require_two_or_more and len(filenames) < 2:
        sys.exit('\nError: two or more input contigs are required')
    if len(filenames) > settings.MAX_INPUT_CONTIGS:
        sys.exit(f'\nError: you cannot have more than {settings.MAX_INPUT_CONTIGS} input contigs')
    log(f'Input contigs:')
    contig_names = set()
    for f in filenames:
        contig_type = get_sequence_file_type(f)
        if contig_type != 'FASTA':
            sys.exit(f'\nError: input contig file ({f}) is not in FASTA format')
        seqs = load_fasta(f)
        if len(seqs) == 0:
            sys.exit(f'\nError: input contig file ({f}) contains no sequences')
        if len(seqs) > 1:
            sys.exit(f'\nError: input contig file ({f}) contains multiple sequences')
        contig_name = seqs[0][0]
        if contig_name in contig_names:
            sys.exit(f'\nError: duplicate contig name: {contig_name}')
        contig_names.add(contig_name)
        contig_len = len(seqs[0][1])
        log(f'  {f} ({contig_len:,} bp)')
    log()


def get_contigs_from_cluster_dir(cluster_dir):
    contig_dir = cluster_dir / '1_contigs'
    if not contig_dir.is_dir():
        sys.exit(f'\nError: contig directory ({contig_dir}) does not exist')
    return sorted(contig_dir.glob('*.fasta'))


def check_cluster_directory(directory):
    if directory.is_file():
        sys.exit(f'\nError: output directory ({directory}) already exists as a file')
    if not directory.is_dir():
        sys.exit(f'\nError: output directory ({directory}) does not exist')


def check_required_software():
    log('Checking required software:')
    check_minimap2()
    log()


def save_seqs_to_fasta(seqs, filename):
    seq_word = 'sequence' if len(seqs) == 1 else 'sequences'
    log(f'Saving {seq_word} to file: {filename}')
    with open(filename, 'wt') as fasta:
        for name, seq in seqs.items():
            fasta.write(f'>{name}\n')
            fasta.write(f'{seq.upper()}\n')
    log()


def print_identity_matrix(seqs, percent_identities, min_allowed_identity):
    log('Overall pairwise identities:')
    seq_names = sorted(seqs.keys())
    max_seq_name_len = max(len(x) for x in seq_names)
    failed = False
    for a in seq_names:
        log('  ' + a, end=':')
        log(' ' * (max_seq_name_len - len(a)), end=' ')
        for b in seq_names:
            if a == b:
                identity = 100.0
            else:
                identity = percent_identities[(a, b)]
            identity_str = f'{identity:.3f}%'.rjust(8)
            if a == b:
                log(dim(identity_str), end='')
            elif identity < min_allowed_identity:
                log(red(identity_str), end='')
                failed = True
            else:
                log(identity_str, end='')
            if b != seq_names[-1]:  # if not the last one in the row
                log('  ', end='')
        log()
    log()
    if failed:
        quit_with_error(f'Error: some pairwise identities are below the minimum allowed value '
                        f'of {min_allowed_identity}%. Please remove offending sequences or lower '
                        f'the --min_identity threshold and try again.')


def print_worst_1kbp_matrix(seqs, worst_1kbp_identities, min_allowed_1kbp_identity):
    log('Worst-1kbp pairwise identities:')
    seq_names = sorted(seqs.keys())
    max_seq_name_len = max(len(x) for x in seq_names)
    failed = False
    for a in seq_names:
        log('  ' + a, end=':')
        log(' ' * (max_seq_name_len - len(a)), end=' ')
        for b in seq_names:
            if a == b:
                identity = 100.0
            else:
                identity = worst_1kbp_identities[(a, b)]
            identity_str = f'{identity:.1f}%'.rjust(6)
            if a == b:
                log(dim(identity_str), end='')
            elif identity < min_allowed_1kbp_identity:
                log(red(identity_str), end='')
                failed = True
            else:
                log(identity_str, end='')
            if b != seq_names[-1]:  # if not the last one in the row
                log('  ', end='')
        log()
    log()
    if failed:
        quit_with_error(f'Error: some pairwise worst-1kbp identities are below the minimum '
                        f'allowed value of {min_allowed_1kbp_identity}%. Please remove offending '
                        f'sequences or lower the --min_1kbp_identity threshold and try again.')
