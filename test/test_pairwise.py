"""
This module contains some tests for Trycycler. To run them, execute `python3 -m pytest` from the
root Trycycler directory.

Copyright 2019 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Minipolish

This file is part of Minipolish. Minipolish is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Minipolish is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Minipolish.
If not, see <http://www.gnu.org/licenses/>.
"""

import trycycler.pairwise


def test_get_pairwise_coordinates_1():
    """
    A: ACTGACTA--CTA
    B: ACCG--TACGCTA
       ==X=II==DD===
    """
    cigar = '2=1X1=2I2=2D3='
    seq_lengths = {'A': 11, 'B': 11}
    a_to_b, b_to_a = trycycler.pairwise.get_pairwise_coordinates('A', 'B', cigar, seq_lengths)

    assert a_to_b == [0, 1, 2, 3, None, None, 4, 5, 8, 9, 10]
    assert b_to_a == [0, 1, 2, 3, 6, 7, None, None, 8, 9, 10]


def test_get_pairwise_alignments_1():
    """
    A: TGGTGTTACACTGCGGGCGAACCGTTCTGATACGTTCTTTTCATTGGTAC
    B: TGGTGTTACACTGCGGGCGAACC-TTCTGATACGTTCTTTTCATTGGTAC
       =======================I==========================
       23=1I26=
    """
    a = 'TGGTGTTACACTGCGGGCGAACCGTTCTGATACGTTCTTTTCATTGGTAC'
    b = 'TGGTGTTACACTGCGGGCGAACCTTCTGATACGTTCTTTTCATTGGTAC'
    seqs = {'A': a, 'B': b}
    pairwise_coordinates = trycycler.pairwise.get_pairwise_alignments(seqs)
    assert pairwise_coordinates[('A', 'B')] == list(range(0, 23)) + [None] + list(range(23, 49))
    assert pairwise_coordinates[('B', 'A')] == list(range(0, 23)) + list(range(24, 50))
