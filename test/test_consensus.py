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


def test_get_consensus_seq_1():
    """
    Sequence A is always better and becomes the consensus.
    """
    msa_seqs = {'A': 'ACG-ACTGT',
                'B': 'ACGCAC-GT'}
    per_base_scores = {'A': '222222222',
                       'B': '111111111'}
    consensus_seq = trycycler.consensus.get_consensus_seq(msa_seqs, per_base_scores)
    assert consensus_seq == 'ACGACTGT'


def test_get_consensus_seq_2():
    """
    Sequence B is always better and becomes the consensus.
    """
    msa_seqs = {'A': 'ACG-ACTGT',
                'B': 'ACGCAC-GT'}
    per_base_scores = {'A': '111111111',
                       'B': '222222222'}
    consensus_seq = trycycler.consensus.get_consensus_seq(msa_seqs, per_base_scores)
    assert consensus_seq == 'ACGCACGT'


def test_get_consensus_seq_3():
    """
    Consensus is a blend of A and B (preferring deletions).
    """
    msa_seqs = {'A': 'ACG-ACTGT',
                'B': 'ACGCAC-GT'}
    per_base_scores = {'A': '222221111',
                       'B': '111112222'}
    consensus_seq = trycycler.consensus.get_consensus_seq(msa_seqs, per_base_scores)
    assert consensus_seq == 'ACGACGT'


def test_get_consensus_seq_4():
    """
    Consensus is a blend of A and B (preferring insertions).
    """
    msa_seqs = {'A': 'ACG-ACTGT',
                'B': 'ACGCAC-GT'}
    per_base_scores = {'A': '111112222',
                       'B': '222221111'}
    consensus_seq = trycycler.consensus.get_consensus_seq(msa_seqs, per_base_scores)
    assert consensus_seq == 'ACGCACTGT'
