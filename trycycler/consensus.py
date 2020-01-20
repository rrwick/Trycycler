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


from .log import log, section_header, explanation


def get_consensus_seq(seqs, per_base_scores, pairwise_alignments):
    section_header('Consensus sequence')
    explanation('Trycycler now builds the consensus sequence by switching between contigs '
                'to stay on whichever has the highest local score.')

    seq_names = list(seqs.keys())
    other_seq_names = get_other_seq_names(seq_names)

    current_seq_name = choose_starting_sequence(seqs, per_base_scores)
    current_seq = seqs[current_seq_name]
    current_pos = 0

    consensus_seq = []
    counts = {n: 0 for n in seq_names}
    log('Consensus sequence composition:')
    while True:
        current_score = per_base_scores[current_seq_name][current_pos]
        best_other_score, best_other_seq_name, best_other_pos = 0, None, None
        for other_seq_name in other_seq_names[current_seq_name]:
            pairwise = pairwise_alignments[(current_seq_name, other_seq_name)]
            other_pos = pairwise[current_pos]
            if other_pos is not None:
                other_score = per_base_scores[other_seq_name][other_pos]
                if other_score > best_other_score:
                    best_other_score = other_score
                    best_other_seq_name = other_seq_name
                    best_other_pos = other_pos
        if best_other_score > current_score:
            current_seq_name = best_other_seq_name
            current_seq = seqs[current_seq_name]
            current_pos = best_other_pos
            log_proportion(counts)
        consensus_seq.append(current_seq[current_pos])
        counts[current_seq_name] += 1

        current_pos += 1
        if current_pos >= len(current_seq):
            break

    log_proportion(counts)
    log('\n')
    return ''.join(consensus_seq)


def choose_starting_sequence(seqs, per_base_scores):
    best_seq_name, best_score = None, 0
    for seq_name in seqs.keys():
        starting_score = per_base_scores[seq_name][0]
        if best_seq_name is None or starting_score > best_score:
            best_seq_name = seq_name
            best_score = starting_score
    return best_seq_name


def get_other_seq_names(seq_names):
    """
    Builds a dictionary where each seq name gives a list of the other seq names.
    """
    other_seq_names = {}
    for seq_name in seq_names:
        other_seq_names[seq_name] = [n for n in seq_names if n != seq_name]
    return other_seq_names


def log_proportion(counts):
    total = sum(counts.values())
    proportions = []
    for seq_name, count in counts.items():
        proportion = 100.0 * count / total
        proportions.append(f'{seq_name}: {proportion:.2f}%')
    log('\r  ' + ', '.join(proportions), end='    ')
