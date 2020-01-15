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

from .alignment import align_reads_to_seq
from .log import log, section_header, explanation
from . import settings


def get_per_base_scores(seqs, reads, circular, threads):
    section_header('Per-base quality scores')
    explanation('Trycycler now aligns all reads to each sequence and uses the alignments to '
                'create per-base quality scores for the entire sequence.')
    per_base_scores = {}
    for seq_name, seq in seqs.items():
        log(f'Aligning reads to sequence {seq_name}:')
        per_base_scores[seq_name] = \
            get_one_seq_per_base_scores(seq, reads, circular, threads)
        log()


def get_one_seq_per_base_scores(seq, reads, circular, threads):
    seq_len = len(seq)

    # For circular sequences, alignments are done to a doubled version of the sequence, to allow
    # for alignments around the junction point. We can then toss out any alignments that reside
    # entirely in the second half.
    if circular:
        ref_seq = seq + seq
        alignments = align_reads_to_seq(reads, ref_seq, threads)
        alignments = [a for a in alignments if a.ref_start < seq_len]
    else:
        ref_seq = seq
        alignments = align_reads_to_seq(reads, ref_seq, threads)

    log(f'  {len(alignments):,} alignments')

    per_base_scores = [0] * len(ref_seq)

    for i, a in enumerate(alignments):
        log(f'\r  calculating alignment scores: {i+1} / {len(alignments)}', end='')
        alignment_scores = get_alignment_scores(a)
        for j, s in enumerate(alignment_scores):
            ref_pos = a.ref_start + j
            per_base_scores[ref_pos] += s
    log()

    # If the sequence was doubled, we now have to undouble it by collapsing the two halves of the
    # per-base scores together.
    if circular:
        non_doubled_per_base_scores = [0] * seq_len
        for i in range(seq_len):
            score_1 = per_base_scores[i]
            score_2 = per_base_scores[i + seq_len]
            non_doubled_per_base_scores[i] = max(score_1, score_2)
        per_base_scores = non_doubled_per_base_scores

    total_score = sum(per_base_scores)
    log(f'  total score = {total_score:,}')
    return per_base_scores


def get_alignment_scores(a):
    # The expanded CIGARs has just the four characters (=, X, I and D) repeated (i.e. no numbers).
    expanded_cigar = a.get_expanded_cigar()

    # The simplified expanded CIGAR has only two characters: = for match and X for everything else.
    simplified_expanded_cigar = expanded_cigar.replace('I', 'X').replace('D', 'X')

    # The pass/fail string has two characters: P for good regions and F for bad regions.
    pass_fail = get_pass_fail_string(simplified_expanded_cigar)

    # We now make the score for each position of the expanded CIGAR. The score increases with
    # matches and resets to zero at fail regions. This is done in both forward and reverse
    # directions.
    forward_scores = get_cigar_scores_forward(expanded_cigar, pass_fail)
    reverse_scores = get_cigar_scores_reverse(expanded_cigar, pass_fail)

    # To make the final scores, we combine the forward and reverse scores, taking the minimum of
    # each. We also drop any insertion positions, so the scores match up with the corresponding
    # range of the reference sequence.
    final_scores = []
    assert len(forward_scores) == len(reverse_scores)
    for i, f in enumerate(forward_scores):
        r = reverse_scores[i]
        if expanded_cigar[i] != 'I':
            final_scores.append(min(f, r))
    assert len(final_scores) == a.ref_end - a.ref_start
    return final_scores


def get_pass_fail_string(simplified_expanded_cigar):
    pass_fail = ['P'] * len(simplified_expanded_cigar)
    i, j = 0, settings.SPLIT_ALIGNMENT_BAD_WINDOW
    while j <= len(simplified_expanded_cigar):
        cigar_window = simplified_expanded_cigar[i:j]
        if cigar_window.count('X') >= settings.SPLIT_ALIGNMENT_BAD_THRESHOLD:
            for k in range(i, j):
                pass_fail[k] = 'F'
        i += 1
        j += 1
    return ''.join(pass_fail)


def get_cigar_scores_forward(expanded_cigar, pass_fail):
    scores = [0] * len(expanded_cigar)
    score = 0
    for i in range(len(expanded_cigar)):
        if pass_fail[i] == 'F':
            score = 0
        elif expanded_cigar[i] == '=':
            score += 1
        scores[i] = score
    return scores


def get_cigar_scores_reverse(expanded_cigar, pass_fail):
    scores = [0] * len(expanded_cigar)
    score = 0
    for i in range(len(expanded_cigar) - 1, -1, -1):  # loop through indices backwards
        if pass_fail[i] == 'F':
            score = 0
        elif expanded_cigar[i] == '=':
            score += 1
        scores[i] = score
    return scores
