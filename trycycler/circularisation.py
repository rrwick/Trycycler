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

import sys

from .alignment import align_a_to_b
from .log import log, section_header, explanation
from . import settings


def circularise(seqs, reads):
    section_header('Circularisation')
    explanation('Trycycler now compares the contigs to each other to repair any circularisation '
                'issues. After this step, each sequence should be cleanly circularised - i.e. the '
                'first base in the contig immediately follows the last base. Each contig will '
                'be circularised by looking for the position of its start and end in the other '
                'contigs. If necessary, additional sequence will be added or duplicated sequence '
                'will be removed. If there are multiple possible ways to fix a contig\'s '
                'circularisation, then Trycycler will use read alignments to choose the best one.')
    circularised_seqs = {}
    for name in seqs.keys():
        circularised_seq = circularise_one_seq_with_all_others(name, seqs, reads)
        circularised_seqs[name] = circularised_seq
    return circularised_seqs


def circularise_one_seq_with_all_others(name_a, seqs, reads):
    log(f'Circularising {name_a}')
    seq_a = seqs[name_a]
    candidate_seqs = set()
    for name_b, seq_b in seqs.items():
        if name_a == name_b:
            continue
        circularised_seq = circularise_one_seq_with_one_other(seq_a, seq_b, name_a, name_b)
        if circularised_seq is not None:
            candidate_seqs.add(circularised_seq)
    candidate_seqs = list(candidate_seqs)

    if len(candidate_seqs) == 0:
        log()
        sys.exit(f'Error: failed to circularise sequence {name_a}')

    if len(candidate_seqs) == 1:
        log('  circularisation complete')
        circularised_seq = candidate_seqs[0]
    else:  # more than one:
        circularised_seq = choose_best_circularisation(candidate_seqs, reads)
    log()
    return circularised_seq


def circularise_one_seq_with_one_other(seq_a, seq_b, name_a, name_b):
    log(f'  using {name_b}')

    end_alignment, start_alignment = find_end_and_start_in_other_seq(seq_a, seq_b, name_a, name_b)
    if end_alignment is None or start_alignment is None:
        log('    cannot circularise')
        return None

    # If we found them with no gap at all (end followed by start), that implies seq A is already
    # cleanly circularised.
    if end_alignment.ref_end == start_alignment.ref_start:
        log('    no adjustment needed (already circular)')
        return seq_a

    # If we found them with a small gap (end then missing seq then start), that implies seq A has a
    # gapped circularisation and needs some sequence added.
    if start_alignment.ref_start > end_alignment.ref_end and \
            start_alignment.ref_start - end_alignment.ref_end < settings.START_END_GAP_SIZE:
        missing_seq = seq_b[end_alignment.ref_end:start_alignment.ref_start]
        log(f'    circularising by adding {len(missing_seq)} bp of sequence')
        return seq_a + missing_seq


    # If we found them overlapping (end alignment and start alignment covering each other), that
    # implies that seq A has a slightly overlapping circularisation.


    # If we found them in the wrong orientation (start then missing seq then end), that implies
    # that seq A has a big overlapping circularisation.


    log('    unable to circularise')  # TEMP
    return seq_a  # TEMP



def find_end_and_start_in_other_seq(seq_a, seq_b, name_a, name_b):
    start_seq = seq_a[:settings.START_END_SIZE]
    end_seq = seq_a[-settings.START_END_SIZE:]

    # Look for seq A's end sequence in seq B. We should hopefully find one instance.
    log(f'    looking for {name_a}\'s end in {name_b}...   ', end='')
    end_alignments = align_a_to_b(end_seq, seq_b)
    end_alignments = [a for a in end_alignments if a.strand == '+'
                      and a.percent_identity >= settings.START_END_IDENTITY_THRESHOLD
                      and a.query_cov >= settings.START_END_COV_THRESHOLD]
    if len(end_alignments) == 0:
        log('not found')
        return None, None
    elif len(end_alignments) > 1:
        log('multiple hits')
        return None, None
    assert len(end_alignments) == 1
    end_alignment = end_alignments[0]
    log(f'found at position {end_alignment.ref_start}-{end_alignment.ref_end}')

    # Look for seq A's start sequence in seq B. We should hopefully find one instance.
    log(f'    looking for {name_a}\'s start in {name_b}... ', end='')
    start_alignments = align_a_to_b(start_seq, seq_b)
    start_alignments = [a for a in start_alignments if a.strand == '+'
                        and a.percent_identity >= settings.START_END_IDENTITY_THRESHOLD
                        and a.query_cov >= settings.START_END_COV_THRESHOLD]
    if len(start_alignments) == 0:
        log('not found')
        return None, None
    elif len(start_alignments) > 1:
        log('multiple hits')
        return None, None
    assert len(start_alignments) == 1
    start_alignment = start_alignments[0]
    log(f'found at position {start_alignment.ref_start}-{start_alignment.ref_end}')

    return end_alignment, start_alignment


def choose_best_circularisation(candidate_seqs, reads):
    log(f'  choosing best circularisation of {len(candidate_seqs)} alternatives')
    # TODO
    # TODO
    # TODO
    # TODO
    return None  # TEMP
