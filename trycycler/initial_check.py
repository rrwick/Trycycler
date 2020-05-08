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

import sys

from .log import log, section_header, explanation, dim, red, quit_with_error
from .mash import get_mash_dist_matrix


def initial_check(seqs, max_mash_dist, max_length_diff):
    section_header('Initial check of contigs')
    explanation('Before proceeding, Trycycler ensures that the input contigs appear sufficiently '
                'close to each other to make a consensus. If not, the program will quit and the '
                'user must fix the input contigs (make them more similar to each other) or exclude '
                'some before trying again.')

    seq_names = list(seqs.keys())

    log('Relative sequence lengths:')
    length_matrix = get_length_ratio_matrix(seq_names, seqs, max_length_diff)
    check_length_ratios(length_matrix, max_length_diff)

    log('Mash distances:')
    mash_matrix = get_mash_dist_matrix(seq_names, seqs, max_mash_dist)
    check_mash_distances(mash_matrix, max_mash_dist)

    log('Contigs have passed the initial check - they seem sufficiently close to reconcile.')
    log()


def get_length_ratio_matrix(seq_names, seqs, max_length_diff):
    max_seq_name_len = max(len(x) for x in seq_names)
    min_threshold, max_threshold = get_length_thresholds(max_length_diff)
    length_matrix = {}
    for a in seq_names:
        log('  ' + a, end=':')
        log(' ' * (max_seq_name_len - len(a)), end=' ')
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


def check_length_ratios(length_matrix, max_length_diff):
    """
    This function looks at the whole length ratio matrix and quits the program if any value is too
    out of bounds.
    """
    min_pair = min(length_matrix, key=length_matrix.get)
    max_pair = max(length_matrix, key=length_matrix.get)
    min_ratio = length_matrix[min_pair]
    max_ratio = length_matrix[max_pair]
    min_threshold, max_threshold = get_length_thresholds(max_length_diff)
    if min_ratio < min_threshold or max_ratio > max_threshold:
        quit_with_error('Error: there is too much length difference between contigs. You must '
                        'either exclude or repair the offending contig sequences and then try '
                        'running trycycler reconcile again. If one of the sequences is too long, '
                        'it could be due to excessive circularisation overlap, and trimming that '
                        'overlap may allow trycycler reconcile to continue.')


def get_length_thresholds(max_length_diff):
    max_threshold = max_length_diff
    min_threshold = 1.0 / max_threshold
    assert min_threshold < 1.0
    assert max_threshold > 1.0
    return min_threshold, max_threshold


def check_mash_distances(mash_matrix, max_mash_dist):
    """
    This function looks at the whole Mash distance matrix and quits the program if any value is too
    out of bounds.
    """
    max_pair = max(mash_matrix, key=mash_matrix.get)
    max_dist = mash_matrix[max_pair]
    if max_dist > max_mash_dist:
        quit_with_error('Error: there is too much Mash distance between contigs. You should '
                        'exclude the offending contig sequences and then try running trycycler '
                        'reconcile again. Alternatively, if you want to accept the high Mash '
                        'distances, you can increase the value of the --max_mash_dist option.')
