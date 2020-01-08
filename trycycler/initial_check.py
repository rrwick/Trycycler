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

import os
import pathlib
import subprocess
import sys
import tempfile

from .log import log, section_header, explanation, dim, red
from .misc import write_seq_to_fasta, reverse_complement
from . import settings


def initial_sanity_check(seqs):
    section_header('Checking closeness of contigs')
    explanation('Before proceeding, Trycycler ensures that the input contigs appear sufficiently '
                'close to each other to make a consensus. If not, the program will quit and the '
                'user must fix the input contigs (make them more similar to each other) or exclude '
                'some before trying again.')

    seq_names = list(seqs.keys())

    log('Relative sequence lengths:')
    length_matrix = get_length_ratio_matrix(seq_names, seqs)
    check_length_ratios(length_matrix)

    log('Mash distances:')
    mash_matrix = get_mash_dist_matrix(seq_names, seqs)
    check_mash_distances(mash_matrix)

    log('Contigs have passed the initial check - they seem sufficiently close to reconcile.')
    log()


def get_length_ratio_matrix(seq_names, seqs):
    min_threshold, max_threshold = get_length_thresholds()
    length_matrix = {}
    for a in seq_names:
        log(f'  {a}: ', end='')
        for b in seq_names:
            ratio = len(seqs[a]) / len(seqs[b])
            ratio_str = f'{ratio:.3f}'
            if a == b:
                log(dim(ratio_str), end='')
            elif ratio < min_threshold or ratio > max_threshold:
                log(red(ratio_str), end='')
            else:
                log(ratio_str, end='')
            if b != seq_names[-1]:  # if not the last one in the row
                log('  ', end='')
            length_matrix[(a, b)] = ratio
        log()
    log()
    return length_matrix


def get_mash_dist_matrix(seq_names, seqs):
    # TODO: I could make this faster by pre-sketching all the sequences and running mash on the
    #       sketches.
    mash_matrix = {}
    for a in seq_names:
        log(f'  {a}: ', end='')
        for b in seq_names:
            if a == b:
                distance = 0.0
            else:
                distance = min(get_mash_dist(seqs[a], seqs[b], '+'),
                               get_mash_dist(seqs[a], seqs[b], '-'))
                mash_matrix[(a, b)] = distance
            distance_str = f'{distance:.3f}'

            if a == b:
                log(dim(distance_str), end='')
            elif distance > settings.MASH_DISTANCE_THRESHOLD:
                log(red(distance_str), end='')
            else:
                log(distance_str, end='')
            if b != seq_names[-1]:  # if not the last one in the row
                log('  ', end='')
            mash_matrix[(b, a)] = distance
        log()
    log()
    return mash_matrix


def print_matrix(seq_names, matrix, min_threshold, max_threshold):
    longest_name_len = max(len(n) for n in seq_names)
    for a in seq_names:
        log('  ', end='')
        log(a.ljust(longest_name_len), end='  ')
        val_strings = []
        for b in seq_names:
            val = matrix[(a, b)]
            val_str = f'{val:.3f}'
            if a == b:
                val_strings.append(dim(val_str))
            elif val < min_threshold or val > max_threshold:
                val_strings.append(red(val_str))
            else:
                val_strings.append(val_str)
        log('  '.join(val_strings))
    log()


def check_length_ratios(length_matrix):
    """
    This function looks at the whole length ratio matrix and quits the program if any value is too
    out of bounds.
    """
    min_pair = min(length_matrix, key=length_matrix.get)
    max_pair = max(length_matrix, key=length_matrix.get)
    min_ratio = length_matrix[min_pair]
    max_ratio = length_matrix[max_pair]
    min_threshold, max_threshold = get_length_thresholds()
    if min_ratio < min_threshold or max_ratio > max_threshold:
        sys.exit('Error: there is too much length difference between contigs')


def get_length_thresholds():
    min_threshold = settings.LENGTH_DIFFERENCE_THRESHOLD
    max_threshold = 1.0 / min_threshold
    assert min_threshold < 1.0
    assert max_threshold > 1.0
    return min_threshold, max_threshold


def check_mash_distances(mash_matrix):
    """
    This function looks at the whole Mash distance matrix and quits the program if any value is too
    out of bounds.
    """
    max_pair = max(mash_matrix, key=mash_matrix.get)
    max_dist = mash_matrix[max_pair]
    if max_dist > settings.MASH_DISTANCE_THRESHOLD:
        sys.exit(f'Error: there is too much Mash distance between contigs')


def get_mash_dist(seq_a, seq_b, strand):
    assert strand == '+' or strand == '-'
    if strand == '-':
        seq_b = reverse_complement(seq_b)
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = pathlib.Path(temp_dir)
        temp_a = temp_dir / 'a.fasta'
        temp_b = temp_dir / 'b.fasta'
        write_seq_to_fasta(seq_a, 'A', temp_a)
        write_seq_to_fasta(seq_b, 'B', temp_b)
        with open(os.devnull, 'w') as dev_null:
            out = subprocess.check_output(['mash', 'dist', '-n', str(temp_a), str(temp_b)],
                                          stderr=dev_null)
        out = out.decode()
        parts = out.split('\t')
        return float(parts[2])
