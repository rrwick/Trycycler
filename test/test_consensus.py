"""
This module contains some tests for Trycycler. To run them, execute `python3 -m pytest` from the
root Trycycler directory.

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

import trycycler.consensus
import trycycler.pairwise


def test_choose_starting_sequence_1():
    seqs = {'A': 'ACTACGACTA', 'B': 'ACTACGACTA', 'C': 'ACTACGACTA'}
    per_base_scores = {'A': [10] + [5] * 9, 'B': [5] + [5] * 9, 'C': [1] + [5] * 9}
    best_seq_name = trycycler.consensus.choose_starting_sequence(seqs, per_base_scores)
    assert best_seq_name == 'A'


def test_choose_starting_sequence_2():
    seqs = {'A': 'ACTACGACTA', 'B': 'ACTACGACTA', 'C': 'ACTACGACTA'}
    per_base_scores = {'A': [5] + [5] * 9, 'B': [10] + [5] * 9, 'C': [1] + [5] * 9}
    best_seq_name = trycycler.consensus.choose_starting_sequence(seqs, per_base_scores)
    assert best_seq_name == 'B'


def test_choose_starting_sequence_3():
    seqs = {'A': 'ACTACGACTA', 'B': 'ACTACGACTA', 'C': 'ACTACGACTA'}
    per_base_scores = {'A': [1] + [5] * 9, 'B': [5] + [5] * 9, 'C': [10] + [5] * 9}
    best_seq_name = trycycler.consensus.choose_starting_sequence(seqs, per_base_scores)
    assert best_seq_name == 'C'


def test_get_other_seq_names_1():
    seq_names = ['A', 'B', 'C']
    other_seq_names = trycycler.consensus.get_other_seq_names(seq_names)
    assert other_seq_names == {'A': ['B', 'C'], 'B': ['A', 'C'], 'C': ['A', 'B']}


def test_get_other_seq_names_2():
    seq_names = ['A', 'B', 'C', 'D']
    other_seq_names = trycycler.consensus.get_other_seq_names(seq_names)
    assert other_seq_names == {'A': ['B', 'C', 'D'], 'B': ['A', 'C', 'D'],
                               'C': ['A', 'B', 'D'], 'D': ['A', 'B', 'C']}


def test_get_other_seq_names_3():
    seq_names = ['A', 'B', 'C', 'D', 'E']
    other_seq_names = trycycler.consensus.get_other_seq_names(seq_names)
    assert other_seq_names == {'A': ['B', 'C', 'D', 'E'], 'B': ['A', 'C', 'D', 'E'],
                               'C': ['A', 'B', 'D', 'E'], 'D': ['A', 'B', 'C', 'E'],
                               'E': ['A', 'B', 'C', 'D']}


# These consensus tests use simple alignments (no indels) between three sequences.
#
# A seq: AAAA
# B seq: GGGG
# C seq: TTTT

def test_get_consensus_seq_1():
    seqs = {'A': 'AAAA', 'B': 'GGGG', 'C': 'TTTT'}
    per_base_scores = {'A': [10, 10, 10, 10], 'B': [5, 5, 5, 5], 'C': [1, 1, 1, 1]}
    pairwise_alignments = {('A', 'B'): [0, 1, 2, 3], ('B', 'A'): [0, 1, 2, 3],
                           ('A', 'C'): [0, 1, 2, 3], ('C', 'A'): [0, 1, 2, 3],
                           ('B', 'C'): [0, 1, 2, 3], ('C', 'B'): [0, 1, 2, 3]}
    consensus = trycycler.consensus.get_consensus_seq(seqs, per_base_scores, pairwise_alignments)
    assert consensus == 'AAAA'


def test_get_consensus_seq_2():
    seqs = {'A': 'AAAA', 'B': 'GGGG', 'C': 'TTTT'}
    per_base_scores = {'A': [10, 5, 1, 10], 'B': [5, 10, 5, 5], 'C': [1, 1, 10, 1]}
    pairwise_alignments = {('A', 'B'): [0, 1, 2, 3], ('B', 'A'): [0, 1, 2, 3],
                           ('A', 'C'): [0, 1, 2, 3], ('C', 'A'): [0, 1, 2, 3],
                           ('B', 'C'): [0, 1, 2, 3], ('C', 'B'): [0, 1, 2, 3]}
    consensus = trycycler.consensus.get_consensus_seq(seqs, per_base_scores, pairwise_alignments)
    assert consensus == 'AGTA'


# These consensus tests an alignment with indels between two sequences.
#
# A seq:  AC-TGCTC
# A pos:  01 23456
#
# B seq:  ACCTAC-C
# B pos:  012345 6

def test_get_consensus_seq_3():
    seqs = {'A': 'ACTGCTC', 'B': 'ACCTACC'}
    per_base_scores = {'A': [9, 9, 9, 9, 9, 9, 9], 'B': [1, 1, 1, 1, 1, 1, 1]}
    pairwise_alignments = {('A', 'B'): [0, 1, 3, 4, 5, None, 6],
                           ('B', 'A'): [0, 1, None, 2, 3, 4, 6]}
    consensus = trycycler.consensus.get_consensus_seq(seqs, per_base_scores, pairwise_alignments)
    assert consensus == 'ACTGCTC'  # all A


def test_get_consensus_seq_4():
    seqs = {'A': 'ACTGCTC', 'B': 'ACCTACC'}
    per_base_scores = {'A': [1, 1, 1, 1, 1, 1, 1], 'B': [9, 9, 9, 9, 9, 9, 9]}
    pairwise_alignments = {('A', 'B'): [0, 1, 3, 4, 5, None, 6],
                           ('B', 'A'): [0, 1, None, 2, 3, 4, 6]}
    consensus = trycycler.consensus.get_consensus_seq(seqs, per_base_scores, pairwise_alignments)
    assert consensus == 'ACCTACC'  # all B


def test_get_consensus_seq_5():
    seqs = {'A': 'ACTGCTC', 'B': 'ACCTACC'}
    per_base_scores = {'A': [9, 9, 9, 9, 1, 1, 1], 'B': [1, 1, 1, 1, 9, 9, 9]}
    pairwise_alignments = {('A', 'B'): [0, 1, 3, 4, 5, None, 6],
                           ('B', 'A'): [0, 1, None, 2, 3, 4, 6]}
    consensus = trycycler.consensus.get_consensus_seq(seqs, per_base_scores, pairwise_alignments)
    assert consensus == 'ACTGCC'  # starts with A, ends with B


def test_get_consensus_seq_6():
    seqs = {'A': 'ACTGCTC', 'B': 'ACCTACC'}
    per_base_scores = {'A': [9, 9, 9, 1, 1, 1, 1], 'B': [1, 1, 1, 9, 9, 9, 9]}
    pairwise_alignments = {('A', 'B'): [0, 1, 3, 4, 5, None, 6],
                           ('B', 'A'): [0, 1, None, 2, 3, 4, 6]}
    consensus = trycycler.consensus.get_consensus_seq(seqs, per_base_scores, pairwise_alignments)
    assert consensus == 'ACTACC'  # starts with A, ends with B


def test_get_consensus_seq_7():
    seqs = {'A': 'ACTGCTC', 'B': 'ACCTACC'}
    per_base_scores = {'A': [1, 1, 1, 9, 9, 9, 9], 'B': [9, 9, 9, 9, 1, 1, 1]}
    pairwise_alignments = {('A', 'B'): [0, 1, 3, 4, 5, None, 6],
                           ('B', 'A'): [0, 1, None, 2, 3, 4, 6]}
    consensus = trycycler.consensus.get_consensus_seq(seqs, per_base_scores, pairwise_alignments)
    assert consensus == 'ACCTGCTC'  # starts with B, ends with A


def test_get_consensus_seq_8():
    seqs = {'A': 'ACTGCTC', 'B': 'ACCTACC'}
    per_base_scores = {'A': [1, 1, 1, 1, 9, 9, 9], 'B': [9, 9, 9, 9, 9, 1, 1]}
    pairwise_alignments = {('A', 'B'): [0, 1, 3, 4, 5, None, 6],
                           ('B', 'A'): [0, 1, None, 2, 3, 4, 6]}
    consensus = trycycler.consensus.get_consensus_seq(seqs, per_base_scores, pairwise_alignments)
    assert consensus == 'ACCTACTC'  # starts with B, ends with A


def test_get_consensus_seq_9():
    """
    In this test, each of the three input sequences has a bad region which the consensus should
    avoid.
    """
    seqs = {'A': 'GTGTACCCCGCCACCTCGCCCGTGGCTGACCCTCCTACATAGCCCACGTTCTCTAAGGGAAGTGTGAATG',
            'B': 'GTGTACCCCGCCACCACGCCCGTGGCTGACCCTCGTACATAGCCCACGTTCTCTAAGGGAAGTGTGAATG',
            'C': 'GTGTACCCCGCCACCACGCCCGTGGCTGACCCTCCTACATAGCCCACGTTCTCTCAGGGAAGTGTGAATG'}
    per_base_scores = {'A': [9]*10 + [1]*10 + [9]*50,
                       'B': [9]*30 + [1]*10 + [9]*30,
                       'C': [9]*50 + [1]*10 + [9]*10}
    pairwise_alignments = trycycler.pairwise.get_pairwise_alignments(seqs)
    consensus = trycycler.consensus.get_consensus_seq(seqs, per_base_scores, pairwise_alignments)
    assert consensus == 'GTGTACCCCGCCACCACGCCCGTGGCTGACCCTCCTACATAGCCCACGTTCTCTAAGGGAAGTGTGAATG'

