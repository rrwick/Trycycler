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

import multiprocessing
import pathlib
import statistics
import subprocess
import sys
import tempfile

from .log import log, section_header, explanation
from .misc import load_fasta, count_substrings, get_sequence_file_type
from .software import check_muscle, get_muscle_version


def msa(args):
    welcome_message()
    check_inputs_and_requirements(args)
    seqs = dict(load_fasta(args.cluster_dir / '2_all_seqs.fasta'))

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = pathlib.Path(temp_dir)
        piece_count = partition_sequences(seqs, args.kmer, args.step, args.lookahead, temp_dir)
        run_muscle_all_pieces(temp_dir, args.threads)
        check_muscle_results(temp_dir, piece_count)
        merge_pieces(temp_dir, args.cluster_dir, seqs)


def welcome_message():
    section_header('Starting Trycycler MSA')
    explanation('Trycycler MSA is a tool for conducting global multiple sequence alignment of '
                'contig sequences.')


def partition_sequences(seqs, kmer, step, lookahead, temp_dir):
    section_header('Partitioning sequences')
    explanation('The sequences are now partitioned into smaller chunks to make the multiple '
                'sequence alignment more tractable.')

    seq_names = sorted(seqs.keys())
    first_seq_name = seq_names[0]
    first_seq = seqs[first_seq_name]
    seq_positions = {n: 0 for n in seq_names}

    chunk_sizes = []
    i = 0
    while True:
        log(f'\rpieces: {i+1}', end='')
        fasta_filename = temp_dir / f'{i:012d}.fasta'
        new_positions = find_next_cutoff_positions(seqs, seq_names, first_seq, seq_positions, kmer,
                                                   step, lookahead)
        chunk_sizes.append(new_positions[first_seq_name] - seq_positions[first_seq_name])

        with open(fasta_filename, 'wt') as fasta:
            for n in seq_names:
                fasta.write(f'>{n} {seq_positions[n]}-{new_positions[n]}\n')
                fasta.write(seqs[n][seq_positions[n]:new_positions[n]])
                fasta.write('\n')

        if new_positions[first_seq_name] == len(first_seq):
            break

        seq_positions = new_positions
        i += 1
    log('\n')
    log(f'median piece size: {int(statistics.median(chunk_sizes)):,} bp')
    log(f'max piece size:    {max(chunk_sizes):,} bp')
    log()
    return len(chunk_sizes)


def run_muscle_all_pieces(temp_dir: pathlib.Path, threads):
    section_header('Running Muscle')
    explanation('Trycycler now runs Muscle on each of the pieces to turn them into multiple '
                'sequence alignments.')
    fasta_files = sorted(temp_dir.glob('*.fasta'))
    muscle_version = get_muscle_version()

    parameters = []
    for f in fasta_files:
        input_fasta = str(f)
        output_filename = input_fasta.replace('.fasta', '_msa.fasta')
        parameters.append((input_fasta, output_filename, muscle_version))
    i = 0
    if threads == 1:
        for f in parameters:
            run_muscle_one_piece(f)
            i += 1
            log(f'\rpieces: {i}', end='')
    else:
        with multiprocessing.Pool(threads) as pool:
            for _ in pool.imap_unordered(run_muscle_one_piece, parameters):
                i += 1
                log(f'\rpieces: {i}', end='')
    log('\n')


def run_muscle_one_piece(parameters):
    input_filename, output_filename, muscle_version = parameters
    muscle_output_filename = output_filename.replace('_msa.fasta', '_muscle.out')
    if muscle_version.startswith('3'):
        muscle_command = ['muscle', '-in', input_filename, '-out', output_filename]
    elif muscle_version.startswith('5'):
        muscle_command = ['muscle', '-align', input_filename, '-output', output_filename]
    else:
        assert False
    with open(muscle_output_filename, 'wt') as muscle_output:
        subprocess.run(muscle_command, stderr=muscle_output)


def check_muscle_results(temp_dir: pathlib.Path, piece_count):
    missing_files = []
    for i in range(piece_count):
        muscle_out_filename = temp_dir / f'{i:012d}_msa.fasta'
        if not muscle_out_filename.is_file():
            missing_files.append(muscle_out_filename)
    if missing_files:
        sys.exit(f'Error: MUSCLE failed to complete on {len(missing_files)} of the {piece_count} '
                 f'pieces. Please remove the most divergent sequences from this cluster and then '
                 f'try again.')


def find_next_cutoff_positions(seqs, seq_names, first_seq, seq_positions, kmer_size, step,
                               lookahead):
    lookahead_size = lookahead
    first_seq_pos = seq_positions[seq_names[0]] + step - kmer_size

    while True:
        # If we've left the end of the first sequence, then take each sequence to its end.
        if first_seq_pos + kmer_size > len(first_seq):
            new_positions = {n: len(seqs[n]) for n in seq_names}
            return new_positions

        # Get the k-mer from the first sequence.
        k = first_seq[first_seq_pos:first_seq_pos + kmer_size]

        # Get the next lookahead chunk of each of the sequences.
        lookahead_seqs = {n: seqs[n][seq_positions[n]:seq_positions[n] + lookahead_size]
                          for n in seq_names}

        # Check that the trial k-mer occurs exactly once in each of the lookahead sequences. If
        # that is not the case, then we'll step forward and try with another k-mer.
        kmer_counts = [count_substrings(s, k) for s in lookahead_seqs.values()]
        if not all(c == 1 for c in kmer_counts):
            first_seq_pos += step
            lookahead_size += step * 2
            continue

        # If we got here, then the trial k-mer is good to use! We get its position in each of the
        # lookahead sequences and use that as the new position.
        new_positions = {n: seq_positions[n] + lookahead_seqs[n].find(k) + kmer_size
                         for n in seq_names}
        return new_positions


def merge_pieces(temp_dir: pathlib.Path, cluster_dir, seqs):
    section_header('Merging MSA')
    explanation('Each of the MSA pieces are now merged together and saved to file.')
    msa_fasta_files = sorted(temp_dir.glob('*_msa.fasta'))
    seq_names = sorted(seqs.keys())
    aligned_seq_parts = {n: [] for n in seq_names}
    for f in msa_fasta_files:
        parts = dict(load_fasta(f))
        for n in seq_names:
            aligned_seq_parts[n].append(parts[n].upper())
    aligned_seqs = {}
    for n in seq_names:
        aligned_seqs[n] = ''.join(aligned_seq_parts[n])

    final_seqs = {}
    for n in seq_names:
        final_seqs[n] = ''.join(aligned_seq_parts[n])

    # Sanity check: the MSA sequences should all be the same length
    msa_length = len(final_seqs[seq_names[0]])
    log(f'MSA length: {msa_length:,} bp')
    for n in seq_names:
        assert len(final_seqs[n]) == msa_length

    # Sanity check: the MSA sequences should match the original sequences.
    for n in seq_names:
        msa_minus_dashes = final_seqs[n].replace('-', '')
        assert seqs[n].upper() == msa_minus_dashes

    # Save the full MSA to file.
    final_msa_fasta_filename = cluster_dir / '3_msa.fasta'
    with open(final_msa_fasta_filename, 'wt') as final_msa_fasta_file:
        for n in seq_names:
            final_msa_fasta_file.write(f'>{n}\n')
            final_msa_fasta_file.write(final_seqs[n])
            final_msa_fasta_file.write('\n')
    log(f'Saving to: {final_msa_fasta_filename}')
    log()


def check_inputs_and_requirements(args):
    check_cluster_directory(args.cluster_dir)
    check_input_sequences(args.cluster_dir)
    check_required_software()


def check_cluster_directory(directory: pathlib.Path):
    if directory.is_file():
        sys.exit(f'\nError: output directory ({directory}) already exists as a file')
    if not directory.is_dir():
        sys.exit(f'\nError: output directory ({directory}) does not exist')


def check_required_software():
    log('Checking required software:')
    check_muscle()
    log()


def check_input_sequences(cluster_dir):
    f = cluster_dir / '2_all_seqs.fasta'

    contig_type = get_sequence_file_type(f)
    if contig_type != 'FASTA':
        sys.exit(f'\nError: input sequence file ({f}) is not in FASTA format')

    seqs = load_fasta(f)
    if len(seqs) < 2:
        sys.exit(f'\nError: input file ({f}) must have two or more sequences')

    log('Input sequences:')
    for name, seq in seqs:
        log(f'  {name}: {len(seq):,} bp')
    log()
