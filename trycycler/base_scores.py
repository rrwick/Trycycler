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

import matplotlib
import matplotlib.pyplot as plt
import multiprocessing
import re

from .alignment import align_reads_to_seq
from .log import log, section_header, explanation
from .misc import means_of_slices


def get_per_base_scores(seqs, cluster_dir, circular, threads, plot_qual,
                        plot_max=None, consensus=False):
    if consensus:
        section_header('Consensus quality scores')
        explanation('Trycycler now aligns all reads to the consensus sequence to create per-base '
                    'quality scores like it did for the input contigs, for comparison.')
    else:
        section_header('Per-base quality scores')
        explanation('Trycycler now aligns all reads to each sequence and uses the alignments to '
                    'create per-base quality scores for the entire sequence.')

    reads = cluster_dir / '4_reads.fastq'

    per_base_scores = {}
    for seq_name, seq in seqs.items():
        log(f'Aligning reads to sequence {seq_name}:')
        scores, total = get_one_seq_per_base_scores(seq, reads, circular, threads)
        per_base_scores[seq_name] = scores
        if plot_qual:
            averaging_window = max(1, min(100, len(scores) // 2000))
            plot_max = plot_per_base_scores(seq_name, scores, cluster_dir, plot_max, total,
                                            averaging_window=averaging_window)
        log()

    return per_base_scores, plot_max


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

    with multiprocessing.Pool(threads) as pool:
        i = 0
        for alignment_scores, a in pool.imap_unordered(get_alignment_scores, alignments):
            i += 1
            log(f'\r  calculating alignment scores: {i:,} / {len(alignments):,}', end='')
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
            if score_1 > score_2:
                non_doubled_per_base_scores[i] = score_1
            else:
                non_doubled_per_base_scores[i] = score_2
        per_base_scores = non_doubled_per_base_scores

    total_score = sum(per_base_scores)
    log(f'  total score = {total_score:,}')
    return per_base_scores, total_score


def get_alignment_scores(a):
    expanded_cigar = get_expanded_cigar(a)

    # We now make the score for each position of the expanded CIGAR. The score increases with
    # matches and resets to zero at fail regions. This is done in both forward and reverse
    # directions.
    forward_scores = get_cigar_scores(expanded_cigar, forward=True)
    reverse_scores = get_cigar_scores(expanded_cigar, forward=False)

    # To make the final scores, we combine the forward and reverse scores, taking the minimum of
    # each. We also drop any insertion positions, so the scores match up with the corresponding
    # range of the reference sequence.
    return combine_forward_and_reverse_scores(a, forward_scores, reverse_scores, expanded_cigar), a


def get_expanded_cigar(a):
    """
    An expanded CIGARs has just the four characters (=, X, I and D) repeated (i.e. no numbers).
    Here I store it as integers (0: =, 1: X, 2: I, 3: D) and return it as a Numpy array
    """
    cigar_parts = re.findall(r'\d+[IDX=]', a.cigar)
    cigar_parts = [(int(c[:-1]), c[-1]) for c in cigar_parts]
    expanded_cigar_size = sum(p[0] for p in cigar_parts)
    expanded_cigar = [0] * expanded_cigar_size

    i = 0
    for num, letter in cigar_parts:
        if letter == '=':
            v = 0
        elif letter == 'X':
            v = 1
        elif letter == 'I':
            v = 2
        elif letter == 'D':
            v = 3
        else:
            assert False
        for _ in range(num):
            expanded_cigar[i] = v
            i += 1
    assert i == expanded_cigar_size

    return expanded_cigar


def get_cigar_scores(expanded_cigar, forward=True):
    scores = [0] * len(expanded_cigar)
    score = 0
    if forward:
        iterator = range(len(expanded_cigar))
    else:
        iterator = range(len(expanded_cigar) - 1, -1, -1)

    for i in iterator:

        # Matches increase the score.
        if expanded_cigar[i] == 0:
            score += 1

        # Mismatches decrease the score.
        elif expanded_cigar[i] == 1:
            score -= 1

        # Insertions/deletions decrease the score more.
        else:
            score -= 2

        if score < 0:
            score = 0
        scores[i] = score
    return scores


def combine_forward_and_reverse_scores(a, forward_scores, reverse_scores, expanded_cigar):
    """
    The combined scores are the minimums values at each position, excluding insertion positions.
    """
    final_scores = [0] * (a.ref_end - a.ref_start)
    assert len(forward_scores) == len(reverse_scores)
    j = 0  # index in final_scores
    for i, f in enumerate(forward_scores):
        r = reverse_scores[i]
        if expanded_cigar[i] != 2:   # if not an insertion
            if f < r:
                final_scores[j] = f
            else:
                final_scores[j] = r
            j += 1
    assert len(final_scores) == j
    return final_scores


class MyAxes(matplotlib.axes.Axes):
    name = 'MyAxes'

    def drag_pan(self, button, _, x, y):
        matplotlib.axes.Axes.drag_pan(self, button, 'x', x, y)  # pretend key=='x'


matplotlib.projections.register_projection(MyAxes)


def plot_per_base_scores(seq_name, per_base_scores, out_dir, plot_max, total,
                         averaging_window=100):
    if plot_max is None:
        plot_max = max(per_base_scores) * 1.2
    positions = list(range(len(per_base_scores)))

    score_means = list(means_of_slices(per_base_scores, averaging_window))
    position_means = list(means_of_slices(positions, averaging_window))

    fig, ax1 = plt.subplots(1, 1, figsize=(12, 4), subplot_kw={'projection': 'MyAxes'}, dpi=300)
    ax1.plot(position_means, score_means, '-', color='#8F0505', linewidth=1)

    plt.xlabel('contig position')
    plt.ylabel('quality score')
    plt.title(f'{seq_name}, total = {total:,}')
    ax1.set_xlim([0, len(per_base_scores)])
    ax1.set_ylim([0, plot_max])

    plot_dir = out_dir / 'plots'
    plot_dir.mkdir(exist_ok=True)
    if seq_name == 'consensus':
        plot_filename = plot_dir / 'consensus.png'
    else:
        plot_filename = plot_dir / ('seq_' + seq_name + '.png')
    plt.savefig(str(plot_filename))

    return plot_max


def add_indels_to_per_base_scores(msa_seqs, per_base_scores):
    """
    This function takes the per-base scores as input and outputs a similar set of per-base scores,
    but this time with the MSA indels included. All output per-base score lists will therefore be
    the same length.
    For example:
      seq:    ACGACTAGCTACACG  ->  ACGACT--AGCTAC-ACG
      scores: 356289238951365  ->  356289222389511365
    """
    msa_per_base_scores = {}
    for seq_name, msa_seq in msa_seqs.items():

        # First we make the scores for the sequence using None as a placeholder for the indels.
        i = 0
        scores = []
        for base in msa_seq:
            if base == '-':
                scores.append(None)
            else:
                scores.append(per_base_scores[seq_name][i])
                i += 1

        # Sanity checks.
        assert i == len(per_base_scores[seq_name])
        assert len(scores) == len(msa_seq)

        # Then we replace the Nones with the lowest of the neighbouring scores.
        forward_scores, reverse_scores, final_scores = [], [], []
        current_score = None
        for i, score in enumerate(scores):
            if score is not None:
                current_score = score
            forward_scores.append(current_score)
        current_score = None
        for i in range(len(scores) - 1, -1, -1):
            score = scores[i]
            if score is not None:
                current_score = score
            reverse_scores.append(current_score)
        reverse_scores = reverse_scores[::-1]
        assert len(forward_scores) == len(reverse_scores)
        for i in range(len(scores)):
            forward_score, reverse_score = forward_scores[i], reverse_scores[i]
            assert forward_score is not None or reverse_score is not None
            if forward_score is None:
                final_scores.append(reverse_score)
            elif reverse_score is None:
                final_scores.append(forward_score)
            else:  # neither are None
                if msa_seq[i] == '-':
                    if forward_score < reverse_score:
                        final_scores.append(forward_score)
                    else:
                        final_scores.append(reverse_score)
                else:
                    assert forward_score == reverse_score
                    final_scores.append(forward_score)
        assert len(forward_scores) == len(final_scores)

        msa_per_base_scores[seq_name] = final_scores

    return msa_per_base_scores
