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

import collections

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

    candidate_seqs, fail_reasons = [], []
    candidate_seq_counts = collections.defaultdict(int)
    for name_b, seq_b in seqs.items():
        if name_a != name_b:
            candidate_seq, fail_reason = \
                circularise_seq_with_another(seq_a, seq_b, name_a, name_b, args)
            candidate_seqs.append(candidate_seq)
            candidate_seq_counts[candidate_seq] += 1
            fail_reasons.append(fail_reason)
    candidate_seqs = [s for s in candidate_seqs if s is not None]
    candidate_seqs = remove_duplicates(candidate_seqs)

    if len(candidate_seqs) == 0:
        circularised_seq = None
        quit_with_error(get_fail_message(name_a, fail_reasons))

    elif len(candidate_seqs) == 1:
        circularised_seq = candidate_seqs[0]

    else:  # more than one:
        circularised_seq = choose_best_circularisation(candidate_seqs, candidate_seq_counts,
                                                       args.reads, args.threads)

    log(f'  circularisation complete ({len(circularised_seq):,} bp)')
    log()
    return circularised_seq


def get_fail_message(name_a, fail_reasons):
    fail_message = f'Error: failed to circularise sequence {name_a}.'
    fail_reasons = set(x for x in fail_reasons if x is not None)
    if len(fail_reasons) > 1:
        return f'Error: failed to circularise sequence {name_a} for multiple reasons. You must ' \
               f'either repair this sequence or exclude it and then try running trycycler ' \
               f'reconcile again.'
    elif 'end not found' in fail_reasons:
        return f'Error: failed to circularise sequence {name_a} because its end could not be ' \
               f'found in other sequences. You can either trim some sequence off the ' \
               f'end of {name_a} or exclude the sequence altogether and try again.'
    elif 'start not found' in fail_reasons:
        return f'Error: failed to circularise sequence {name_a} because its start could not be ' \
               f'found in other sequences. You can either trim some sequence off the ' \
               f'start of {name_a} or exclude the sequence altogether and try again.'
    elif 'same start/end' in fail_reasons:
        return f'Error: failed to circularise sequence {name_a} because it had the same ' \
               f'start/end as the other sequences. I.e. the circularisations were not different ' \
               f'and therefore could not reconcile each other.'
    elif 'multiple possibilities' in fail_reasons:
        return f'Error: failed to circularise sequence {name_a} because its start/end sequences ' \
               f'were found in multiple ambiguous places in other sequences. This is likely ' \
               f'because {name_a} starts/ends in a repetitive region. You can either manually ' \
               f'repair its circularisation (and ensure it does not start/end in a repetitive ' \
               f'region) or exclude the sequence altogether and try again.'
    elif 'too much extra' in fail_reasons:
        return f'Error: failed to circularise sequence {name_a} because it contained too much ' \
               f'start-end overlap. You can either manually trim the contig sequence, ' \
               f'increase the value of the --max_trim_seq/--max_trim_seq_percent parameters or ' \
               f'exclude the sequence altogether and try again.'
    elif 'too much missing' in fail_reasons:
        return f'Error: failed to circularise sequence {name_a} because it contained too large ' \
               f'of a start-end gap. You can either increase the value of the --max_add_seq/' \
               f'--max_add_seq_percent parameters or exclude the sequence altogether and try again.'
    return fail_message


def circularise_seq_with_another(seq_a, seq_b, name_a, name_b, args):
    log(f'  using {name_b}:')

    max_add_seq_relative = int(round(args.max_add_seq_percent * len(seq_a) / 100.0))
    max_add_seq = min(args.max_add_seq, max_add_seq_relative)

    max_trim_seq_relative = int(round(args.max_trim_seq_percent * len(seq_a) / 100.0))
    max_trim_seq = min(args.max_trim_seq, max_trim_seq_relative)

    end_alignment, start_alignment, fail_reason = \
        find_end_and_start(seq_a, seq_b, name_a, name_b, args)
    if end_alignment is None or start_alignment is None:
        start_or_end = 'start' if fail_reason == 'start not found' else 'end'
        log(f'    unable to circularise: {name_a}\'s {start_or_end} could not be found in {name_b}')
        return None, fail_reason

    # If we found them with no gap at all (end immediately followed by start), that implies seq A
    # is already cleanly circularised.
    if end_alignment.ref_end == start_alignment.ref_start:
        log(f'    no adjustment needed ({name_a} is already circular)')
        return seq_a, None

    # If the starts or ends of the two sequences are the same, then we can't circularise because
    # there's no difference.
    if start_alignment.ref_start == 0 or end_alignment.ref_end == len(seq_b):
        log(f'    unable to circularise: {name_a}\'s start/end is the same as {name_b}\'s '
            f'start/end')
        return None, 'same start/end'

    # If we found them with a small gap (end then missing seq then start), that implies seq A has a
    # gapped circularisation and needs some sequence added.
    if start_alignment.ref_start > end_alignment.ref_end:
        missing_seq = seq_b[end_alignment.ref_end:start_alignment.ref_start]
        if len(missing_seq) < max_add_seq:
            log(f'    circularising {name_a} by adding {len(missing_seq)} bp of sequence from'
                f' {name_b} ({end_alignment.ref_end}-{start_alignment.ref_start})')
            return seq_a + missing_seq, None
        else:
            log(f'    unable to circularise: {name_a} requires {len(missing_seq)} bp to be added ' 
                f'but settings only allow {max_add_seq} bp')
            return None, 'too much missing'

    # If we got here, then it seems seq A has overlapping circularisation and some sequence will
    # need to be removed. To figure out how much to trim, we take the part of seq B which precedes
    # the start alignment and align it back to seq A.
    pre_start_alignment = find_pre_start_alignment(seq_a, seq_b, name_a, name_b, start_alignment,
                                                   args.verbose)
    if pre_start_alignment is not None:
        trim_point = pre_start_alignment.ref_end
        trim_amount = len(seq_a) - trim_point
        if trim_amount < max_trim_seq:
            log(f'    circularising {name_a} by trimming {trim_amount} bp of sequence from the end')
            return seq_a[:trim_point], None
        else:
            log(f'    unable to circularise: {name_a} requires {trim_amount} bp to be trimmed ' 
                f'but settings only allow {max_trim_seq} bp')
            return None, 'too much extra'

    log('    unable to circularise: cannot determine trim amount')
    return None, 'other'


def get_start_end_size(seq):
    """
    Don't use a start/end size of more than 10% of the sequence length.
    """
    start_end_size = settings.START_END_SIZE
    if start_end_size > len(seq) / 10.0:
        start_end_size = int(round(len(seq) / 10.0))
    return start_end_size


def find_end_and_start(seq_a, seq_b, name_a, name_b, args):
    start_end_size = get_start_end_size(seq_a)
    seq_a_start_start, seq_a_start_end = 0, start_end_size
    seq_a_end_start, seq_a_end_end = len(seq_a) - start_end_size, len(seq_a)

    start_seq = seq_a[seq_a_start_start:seq_a_start_end]
    end_seq = seq_a[seq_a_end_start:seq_a_end_end]

    # Look for seq A's end sequence in seq B.
    if args.verbose:
        log(f'    looking for {name_a}\'s end ({seq_a_end_start}-{seq_a_end_end}) in {name_b}...')
    end_alignments = align_a_to_b(end_seq, seq_b)
    end_alignments = [a for a in end_alignments if a.strand == '+'
                      and a.percent_identity >= settings.START_END_IDENTITY_THRESHOLD
                      and a.query_cov >= settings.START_END_COV_THRESHOLD]
    if len(end_alignments) == 0:
        if args.verbose:
            log('      not found')
        return None, None, 'end not found'
    if args.verbose:
        for end_alignment in end_alignments:
            log(f'      found at {end_alignment.ref_start}-{end_alignment.ref_end}')

    # Look for seq A's start sequence in seq B.
    if args.verbose:
        log(f'    looking for {name_a}\'s start ({seq_a_start_start}-{seq_a_start_end}) in '
            f'{name_b}... ')
    start_alignments = align_a_to_b(start_seq, seq_b)
    start_alignments = [a for a in start_alignments if a.strand == '+'
                        and a.percent_identity >= settings.START_END_IDENTITY_THRESHOLD
                        and a.query_cov >= settings.START_END_COV_THRESHOLD]
    if len(start_alignments) == 0:
        if args.verbose:
            log('      not found')
        return None, None, 'start not found'
    if args.verbose:
        for start_alignment in start_alignments:
            log(f'      found at {start_alignment.ref_start}-{start_alignment.ref_end}')

    if len(end_alignments) == 1 and len(start_alignments) == 1:
        return end_alignments[0], start_alignments[0], None

    # Since we have multiple possibilities, we will try all start/end combinations.
    if args.verbose:
        log(f'    checking different end-start alignment combinations...')
    assert len(end_alignments) > 1 or len(start_alignments) > 1
    combinations = []
    for end_alignment in end_alignments:
        end_1, end_2 = end_alignment.ref_start, end_alignment.ref_end
        for start_alignment in start_alignments:
            start_1, start_2 = start_alignment.ref_start, start_alignment.ref_end
            distance = abs(end_alignment.ref_end - start_alignment.ref_start)
            if args.verbose:
                log(f'      {end_1}-{end_2}, {start_1}-{start_2} distance={distance}')
            combinations.append((distance, end_alignment, start_alignment))
    combinations = sorted(combinations, key=lambda x: x[0])  # sort from closest to furthest
    assert len(combinations) > 1

    best_second_best_diff = combinations[1][0] - combinations[0][0]
    assert best_second_best_diff >= 0

    if best_second_best_diff >= settings.START_END_MIN_COMBINATION_DIFF:
        _, end_alignment, start_alignment = combinations[0]
        return end_alignment, start_alignment, None
    else:
        return None, None, 'multiple possibilities'


def find_pre_start_alignment(seq_a, seq_b, name_a, name_b, start_alignment, verbose):
    start_end_size = get_start_end_size(seq_a)
    start, end = start_alignment.ref_start - start_end_size, start_alignment.ref_start
    if start < 0:
        if verbose:
            log(f'    cannot search for {name_b}\'s before-start ({start}-{end}) in {name_a}')
        return None
    pre_start_seq = seq_b[start:end]
    if verbose:
        log(f'    looking for {name_b}\'s before-start ({start}-{end}) in {name_a}... ')

    if len(pre_start_seq) != start_end_size:
        if verbose:
            log('      not found')
        return None
    pre_start_alignments = align_a_to_b(pre_start_seq, seq_a)
    pre_start_alignments = [a for a in pre_start_alignments if a.strand == '+'
                            and a.percent_identity >= settings.START_END_IDENTITY_THRESHOLD
                            and a.query_cov >= settings.START_END_COV_THRESHOLD]
    if len(pre_start_alignments) == 0:
        if verbose:
            log('      not found')
        return None
    elif len(pre_start_alignments) == 1:
        pre_start_alignment = pre_start_alignments[0]
    else:
        # If we got here, then there are multiple alignments. We choose the one closest to the end.
        pre_start_alignments = sorted(pre_start_alignments, key=lambda a: a.ref_end)
        pre_start_alignment = pre_start_alignments[-1]

    if verbose:
        log(f'      found at {pre_start_alignment.ref_start}-{pre_start_alignment.ref_end}')
    return pre_start_alignment


def choose_best_circularisation(candidate_seqs, candidate_seq_counts, reads, threads):
    """
    This function chooses between multiple alternative circularisations. It first tries to choose
    the a most-common option: it excludes all just-once options, and if a single option remains,
    that's the one. If that fails to find a winner, then it aligns reads to the option and chooses
    whichever circularisation has the highest total alignment score.
    """
    more_than_one = [seq for seq, count in candidate_seq_counts.items()
                     if count > 1 and seq is not None]
    if len(more_than_one) == 1:
        log('  choosing most common circularisation')
        return more_than_one[0]

    # If there are multiple options each with more than one vote, then we will only consider those
    # for the read-based decision.
    elif len(more_than_one) > 1:
        candidate_seqs = more_than_one

    log(f'  choosing best circularisation of {len(candidate_seqs)} alternatives')
    best_seq, best_score, best_i = None, 0.0, 0
    halfway_point = len(candidate_seqs[0]) // 2
    for i, candidate_seq in enumerate(candidate_seqs):
        log(f'    alternative {i+1} ({len(candidate_seq):,} bp): score = ', end='')
        loop_seq = candidate_seq[:halfway_point] + candidate_seq[halfway_point:]
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
