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

from .alignment import align_a_to_b, align_reads_to_seq, get_best_alignment_per_read
from .log import log, section_header, explanation, quit_with_error
from .misc import remove_duplicates
from . import settings


def circularise(seqs, args):
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
        circularised_seq = circularise_seq_with_others(name, seqs, args)
        circularised_seqs[name] = circularised_seq
    return circularised_seqs


def circularise_seq_with_others(name_a, seqs, args):
    log(f'Circularising {name_a}:')
    seq_a = seqs[name_a]

    candidate_seqs = []
    for name_b, seq_b in seqs.items():
        if name_a != name_b:
            candidate_seqs.append(circularise_seq_with_another(seq_a, seq_b, name_a, name_b, args))
    candidate_seqs = [s for s in candidate_seqs if s is not None]
    candidate_seqs = remove_duplicates(candidate_seqs)

    if len(candidate_seqs) == 0:
        circularised_seq = None
        # TODO: make this message more verbose with suggestions on what to do.
        quit_with_error(f'Error: failed to circularise sequence {name_a}.')

    elif len(candidate_seqs) == 1:
        circularised_seq = candidate_seqs[0]

    else:  # more than one:
        circularised_seq = choose_best_circularisation(candidate_seqs, args.reads, args.threads)
    log(f'  circularisation complete ({len(circularised_seq):,} bp)')
    log()
    return circularised_seq


def circularise_seq_with_another(seq_a, seq_b, name_a, name_b, args):
    log(f'  using {name_b}:')

    max_add_seq_relative = int(round(args.max_add_seq_percent * len(seq_a) / 100.0))
    max_add_seq = min(args.max_add_seq, max_add_seq_relative)

    max_trim_seq_relative = int(round(args.max_trim_seq_percent * len(seq_a) / 100.0))
    max_trim_seq = min(args.max_trim_seq, max_trim_seq_relative)

    end_alignment, start_alignment = find_end_and_start(seq_a, seq_b, name_a, name_b, args)
    if end_alignment is None or start_alignment is None:
        log('    cannot circularise')
        return None

    # If we found them with no gap at all (end immediately followed by start), that implies seq A
    # is already cleanly circularised.
    if end_alignment.ref_end == start_alignment.ref_start:
        log(f'    no adjustment needed ({name_a} is already circular)')
        return seq_a

    # If we found them with a small gap (end then missing seq then start), that implies seq A has a
    # gapped circularisation and needs some sequence added.
    if start_alignment.ref_start > end_alignment.ref_end:
        missing_seq = seq_b[end_alignment.ref_end:start_alignment.ref_start]
        if len(missing_seq) < max_add_seq:
            log(f'    circularising {name_a} by adding {len(missing_seq)} bp of sequence from'
                f' {name_b} ({end_alignment.ref_end}-{start_alignment.ref_start})')
            return seq_a + missing_seq
        else:
            log(f'    unable to circularise: {name_a} requires {len(missing_seq)} bp to be added ' 
                f'but settings only allow {max_add_seq} bp')
            return None

    # If we got here, then it seems seq A has overlapping circularisation and some sequence will
    # need to be removed. To figure out how much to trim, we take the part of seq B which precedes
    # the start alignment and align it back to seq A.
    pre_start_alignment = find_pre_start_alignment(seq_a, seq_b, start_alignment)
    if pre_start_alignment is not None:
        trim_point = pre_start_alignment.ref_end
        trim_amount = len(seq_a) - trim_point
        if trim_amount < max_trim_seq:
            log(f'    circularising {name_a} by trimming {trim_amount} bp of sequence from the end')
            return seq_a[:trim_point]
        else:
            log(f'    unable to circularise: {name_a} requires {trim_amount} bp to be trimmed ' 
                f'but settings only allow {max_trim_seq} bp')
            return None

    log('    unable to circularise')
    return None


def find_end_and_start(seq_a, seq_b, name_a, name_b, args):
    seq_a_start_start, seq_a_start_end = 0, settings.START_END_SIZE
    seq_a_end_start, seq_a_end_end = len(seq_a) - settings.START_END_SIZE, len(seq_a)

    start_seq = seq_a[seq_a_start_start:seq_a_start_end]
    end_seq = seq_a[seq_a_end_start:seq_a_end_end]

    # Look for seq A's end sequence in seq B. We should hopefully find one instance.
    if args.verbose:
        log(f'    looking for {name_a}\'s end in {name_b}...   ', end='')
    end_alignments = align_a_to_b(end_seq, seq_b)
    end_alignments = [a for a in end_alignments if a.strand == '+'
                      and a.percent_identity >= settings.START_END_IDENTITY_THRESHOLD
                      and a.query_cov >= settings.START_END_COV_THRESHOLD]
    if len(end_alignments) == 0:
        if args.verbose:
            log('not found')
        return None, None
    elif len(end_alignments) > 1:
        if args.verbose:
            log('multiple hits')
        return None, None
    assert len(end_alignments) == 1
    end_alignment = end_alignments[0]
    if args.verbose:
        log(f'found at {end_alignment.ref_start}-{end_alignment.ref_end}')

    # Look for seq A's start sequence in seq B. We should hopefully find one instance.
    if args.verbose:
        log(f'    looking for {name_a}\'s start in {name_b}... ', end='')
    start_alignments = align_a_to_b(start_seq, seq_b)
    start_alignments = [a for a in start_alignments if a.strand == '+'
                        and a.percent_identity >= settings.START_END_IDENTITY_THRESHOLD
                        and a.query_cov >= settings.START_END_COV_THRESHOLD]
    if len(start_alignments) == 0:
        if args.verbose:
            log('not found')
        return None, None
    elif len(start_alignments) > 1:
        if args.verbose:
            log('multiple hits')
        return None, None
    assert len(start_alignments) == 1
    start_alignment = start_alignments[0]
    if args.verbose:
        log(f'found at {start_alignment.ref_start}-{start_alignment.ref_end}')

    return end_alignment, start_alignment


def find_pre_start_alignment(seq_a, seq_b, start_alignment):
    pre_start_seq = seq_b[start_alignment.ref_start - settings.START_END_SIZE:
                          start_alignment.ref_start]
    if len(pre_start_seq) != settings.START_END_SIZE:
        return None
    pre_start_alignments = align_a_to_b(pre_start_seq, seq_a)
    pre_start_alignments = [a for a in pre_start_alignments if a.strand == '+'
                            and a.percent_identity >= settings.START_END_IDENTITY_THRESHOLD
                            and a.query_cov >= settings.START_END_COV_THRESHOLD]
    if len(pre_start_alignments) == 0:
        return None
    elif len(pre_start_alignments) > 1:
        return None
    else:
        return pre_start_alignments[0]


def choose_best_circularisation(candidate_seqs, reads, threads):
    """
    This function chooses between multiple alternative circularisations. It does so by aligning
    reads to the junction and choosing whichever circularisation has the highest total alignment
    score.
    """
    log(f'  choosing best circularisation of {len(candidate_seqs)} alternatives')
    best_seq, best_score, best_i = None, 0.0, 0
    for i, candidate_seq in enumerate(candidate_seqs):
        log(f'    alternative {i+1} ({len(candidate_seq):,} bp): score = ', end='')
        loop_seq = (candidate_seq[:settings.CIRCULARISATION_CHOICE_ALIGNMENT_SIZE] +
                    candidate_seq[-settings.CIRCULARISATION_CHOICE_ALIGNMENT_SIZE:])
        alignments = align_reads_to_seq(reads, loop_seq, threads)
        alignments = get_best_alignment_per_read(alignments)
        score = sum(a.alignment_score for a in alignments)
        log(f'{score:,}')
        if score > best_score:
            best_seq = candidate_seq
            best_score = score
            best_i = i
    log(f'    best alternative: {best_i+1}')
    return best_seq
