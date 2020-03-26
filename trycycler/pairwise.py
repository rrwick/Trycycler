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

from .log import log, section_header, explanation


def get_pairwise_alignments(seqs):
    section_header('Pairwise global alignments')
    explanation('Trycycler uses the edlib aligner to get global alignments between all pairs of '
                'sequences. This will allow it to \'jump\' from any sequence to the corresponding '
                'position of any other sequence.')
    seq_names = list(seqs.keys())
    max_seq_name_len = max(len(x) for x in seq_names)
    pairwise_cigars, percent_identities = {}, {}

    for i, a in enumerate(seq_names):
        seq_a = seqs[a]
        for j in range(i+1, len(seq_names)):
            b = seq_names[j]
            seq_b = seqs[b]
            log(' ' * (max_seq_name_len - len(a)) + a, end='')
            log(' vs ', end='')
            log(b + '...' + ' ' * (max_seq_name_len - len(b)), end=' ')

            result = edlib.align(seq_a, seq_b, mode="NW", task="path")
            cigar = result['cigar']
            percent_identity, max_indel = identity_and_max_indel_from_cigar(cigar)
            log(f'{percent_identity:.2f}% identity, max indel = {max_indel}')

            pairwise_cigars[(a, b)] = cigar
            percent_identities[(a, b)] = percent_identity
            percent_identities[(b, a)] = percent_identity
    log()

    return pairwise_cigars, percent_identities


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

