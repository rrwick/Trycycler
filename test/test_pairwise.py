"""
This module contains some tests for Trycycler. To run them, execute `pytest` from the root
Trycycler directory.

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

import pytest

import trycycler.pairwise


def test_get_pairwise_alignments_1():
    # All seqs are the same.
    seqs = {'A': 'ATTAGCCGCTCGCACCACCTTGAAGATCGGCAACACATGCGCTCCTGGAT',
            'B': 'ATTAGCCGCTCGCACCACCTTGAAGATCGGCAACACATGCGCTCCTGGAT',
            'C': 'ATTAGCCGCTCGCACCACCTTGAAGATCGGCAACACATGCGCTCCTGGAT',
            'D': 'ATTAGCCGCTCGCACCACCTTGAAGATCGGCAACACATGCGCTCCTGGAT'}
    pairwise_cigars, percent_identities, max_indels = \
        trycycler.pairwise.get_pairwise_alignments(seqs)

    assert len(pairwise_cigars) == 6
    assert len(percent_identities) == 12
    assert len(max_indels) == 12

    for pair in pairwise_cigars:
        assert pairwise_cigars[pair] == '50='
        assert percent_identities[pair] == 100.0
        assert max_indels[pair] == 0


def test_get_pairwise_alignments_2():
    # Seq B has a single substitution.
    seqs = {'A': 'ATTAGCCGCTCGCACCACCTTGAAGATCGGCAACACATGCGCTCCTGGAT',
            'B': 'ATTAGCCGCTCGCACCACCTTGACGATCGGCAACACATGCGCTCCTGGAT',
            'C': 'ATTAGCCGCTCGCACCACCTTGAAGATCGGCAACACATGCGCTCCTGGAT'}
    pairwise_cigars, percent_identities, max_indels = \
        trycycler.pairwise.get_pairwise_alignments(seqs)

    assert len(pairwise_cigars) == 3
    assert len(percent_identities) == 6
    assert len(max_indels) == 6

    assert pairwise_cigars[('A', 'B')] == '23=1X26='
    assert pairwise_cigars[('A', 'C')] == '50='
    assert pairwise_cigars[('B', 'C')] == '23=1X26='

    assert percent_identities[('A', 'B')] == 98.0
    assert percent_identities[('A', 'C')] == 100.0
    assert percent_identities[('B', 'C')] == 98.0

    assert max_indels[('A', 'B')] == 0
    assert max_indels[('A', 'C')] == 0
    assert max_indels[('B', 'C')] == 0


def test_get_pairwise_alignments_3():
    # Seq B has one insertion and one deletion
    seqs = {'A': 'ATTAGCCGCTCGCACCACCTTGAAGATCGGCAACACATGCGCTCCTGGAT',
            'B': 'ATTAGCCGCTCGTCACCACCTTGAAGATCGGCAACACAGCGCTCCTGGAT',
            'C': 'ATTAGCCGCTCGCACCACCTTGAAGATCGGCAACACATGCGCTCCTGGAT'}
    pairwise_cigars, percent_identities, max_indels = \
        trycycler.pairwise.get_pairwise_alignments(seqs)

    assert len(pairwise_cigars) == 3
    assert len(percent_identities) == 6
    assert len(max_indels) == 6

    assert pairwise_cigars[('A', 'B')] == '12=1D25=1I12='
    assert pairwise_cigars[('A', 'C')] == '50='
    assert pairwise_cigars[('B', 'C')] == '12=1I25=1D12='

    assert percent_identities[('A', 'B')] == pytest.approx(96.078431372549)  # 49 / 51
    assert percent_identities[('A', 'C')] == 100.0
    assert percent_identities[('B', 'C')] == pytest.approx(96.078431372549)  # 49 / 51

    assert max_indels[('A', 'B')] == 1
    assert max_indels[('A', 'C')] == 0
    assert max_indels[('B', 'C')] == 1
