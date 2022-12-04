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

import edlib
import re

from .log import log, section_header, explanation


def get_pairwise_alignments(seqs):
    section_header('Pairwise global alignments')
    explanation('Trycycler uses the edlib aligner to get global alignments between all pairs of '
                'sequences. This can help you to spot any problematic sequences that should be '
                'excluded before continuing. If you see any sequences with notably worse '
                'identities, you can remove them (delete the contig\'s FASTA) and run this '
                'command again.')
    seq_names = list(seqs.keys())
    max_seq_name_len = max(len(x) for x in seq_names)
    pairwise_cigars, percent_identities, worst_1kbp_identities = {}, {}, {}

    for i, a in enumerate(seq_names):
        seq_a = seqs[a]
        for j in range(i+1, len(seq_names)):
            b = seq_names[j]
            seq_b = seqs[b]
            log(' ' * (max_seq_name_len - len(a)) + a, end='')
            log(' vs ', end='')
            log(b + '...' + ' ' * (max_seq_name_len - len(b)), end=' ')

            result = edlib.align(seq_a, seq_b, mode='NW', task='path')
            cigar = result['cigar']
            percent_identity, worst_1kbp = identity_and_worst_1kbp_from_cigar(cigar)
            log(f'{percent_identity:.3f}% overall identity, '
                f'{worst_1kbp:.1f}% worst-1kbp identity')

            pairwise_cigars[(a, b)] = cigar
            percent_identities[(a, b)] = percent_identity
            percent_identities[(b, a)] = percent_identity
            worst_1kbp_identities[(a, b)] = worst_1kbp
            worst_1kbp_identities[(b, a)] = worst_1kbp
    log()

    return pairwise_cigars, percent_identities, worst_1kbp_identities


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


def identity_and_worst_1kbp_from_cigar(cigar):
    percent_identity, _ = identity_and_max_indel_from_cigar(cigar)
    expanded_cigar = get_expanded_cigar(cigar)
    worst_1kbp = worst_window_identity(expanded_cigar, 1000)
    return percent_identity, worst_1kbp


def get_expanded_cigar(cigar):
    expanded_cigar = []
    cigar_parts = re.findall(r'\d+[IDX=]', cigar)
    for p in cigar_parts:
        size = int(p[:-1])
        letter = p[-1]
        expanded_cigar.append(letter * size)
    return ''.join(expanded_cigar)


def worst_window_identity(expanded_cigar, window_size):
    """
    This function returns the worst percent identity in a sliding window across an expanded CIGAR
    string.
    """
    cigar_len = len(expanded_cigar)
    if cigar_len <= window_size:
        return 100.0 * expanded_cigar.count('=') / cigar_len
    start, end = 0, window_size
    window_match_count = expanded_cigar[start:end].count('=')
    min_match_count = window_match_count
    while end < cigar_len:
        if expanded_cigar[start] == '=':
            window_match_count -= 1
        if expanded_cigar[end] == '=':
            window_match_count += 1
        if window_match_count < min_match_count:
            min_match_count = window_match_count
        start += 1
        end += 1
    return 100.0 * min_match_count / window_size
