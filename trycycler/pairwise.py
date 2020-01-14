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

import edlib
import re
import sys

from .log import log, section_header, explanation
from . import settings


def get_pairwise_alignments(seqs):
    section_header('Pairwise global alignments')
    explanation('Trycycler uses the edlib aligner to get global alignments between all pairs of '
                'sequences. This will allow it to \'jump\' from any sequence to the corresponding '
                'position of any other sequence.')
    seq_names = list(seqs.keys())
    seq_lengths = {name: len(seqs[name]) for name in seq_names}
    pairwise_cigars = {}
    worst_identity = 100.0
    for i, a in enumerate(seq_names):
        seq_a = seqs[a]
        for j in range(i+1, len(seq_names)):
            b = seq_names[j]
            seq_b = seqs[b]
            log(f'{a} vs {b}... ', end='')

            result = edlib.align(seq_a, seq_b, mode="NW", task="path")
            cigar = result['cigar']
            percent_identity, max_indel = identity_and_max_indel_from_cigar(cigar)
            log(f'{percent_identity:.2f}% identity, max indel = {max_indel}')
            worst_identity = min(percent_identity, worst_identity)
            pairwise_cigars[(a, b)] = cigar

    if worst_identity < settings.MIN_ALLOWED_PAIRWISE_IDENTITY:
        sys.exit(f'Error: some pairwise identities are below the minimum allowed'
                 f'({settings.MIN_ALLOWED_PAIRWISE_IDENTITY}%). Please remove offending '
                 f'sequences and try again.')

    return get_all_pairwise_coordinates(seq_names, pairwise_cigars, seq_lengths)


def get_all_pairwise_coordinates(seq_names, pairwise_cigars, seq_lengths):
    pairwise_coordinates = {}
    for i, a in enumerate(seq_names):
        for j in range(i+1, len(seq_names)):
            b = seq_names[j]
            cigar = pairwise_cigars[(a, b)]
            coords_1, coords_2 = get_pairwise_coordinates(a, b, cigar, seq_lengths)
            pairwise_coordinates[(a, b)] = coords_1
            pairwise_coordinates[(b, a)] = coords_2
    return pairwise_coordinates


def get_pairwise_coordinates(a, b, cigar, seq_lengths):
    """
    Builds two lists of coordinates for translating seq A positions to seq B positions and vice
    versa.
    """
    a_to_b = [None] * seq_lengths[a]
    b_to_a = [None] * seq_lengths[b]
    cigar_parts = re.findall(r'\d+[IDX=]', cigar)

    i, j = 0, 0
    for p in cigar_parts:
        size = int(p[:-1])
        letter = p[-1]
        if letter == '=' or letter == 'X':
            for _ in range(size):
                a_to_b[i] = j
                b_to_a[j] = i
                i += 1
                j += 1
        elif letter == 'I':  # insertion means in seq A but not in seq B
            for _ in range(size):
                a_to_b[i] = None
                i += 1
        elif letter == 'D':  # deletion means in seq B but not in seq A
            for _ in range(size):
                b_to_a[j] = None
                j += 1
        else:
            assert False

    assert i == seq_lengths[a]
    assert j == seq_lengths[b]

    return a_to_b, b_to_a


def identity_and_max_indel_from_cigar(cigar):
    cigar_parts = re.findall(r'\d+[IDX=]', cigar)
    total, matches, max_indel = 0, 0, 0
    for p in cigar_parts:
        size = int(p[:-1])
        letter = p[-1]
        total += size
        if letter == '=':
            matches += size
        if letter == 'I' or letter == 'D':
            max_indel = max(max_indel, size)
    percent_identity = 100.0 * matches / total
    return percent_identity, max_indel

