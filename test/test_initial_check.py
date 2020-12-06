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

import trycycler.initial_check


def test_get_length_thresholds_1():
    min_threshold, max_threshold = trycycler.initial_check.get_length_thresholds(1.5)
    assert min_threshold == pytest.approx(0.6666666666667)
    assert max_threshold == pytest.approx(1.5)


def test_get_length_thresholds_2():
    min_threshold, max_threshold = trycycler.initial_check.get_length_thresholds(2.0)
    assert min_threshold == pytest.approx(0.5)
    assert max_threshold == pytest.approx(2.0)


def test_get_length_thresholds_3():
    min_threshold, max_threshold = trycycler.initial_check.get_length_thresholds(1.1)
    assert min_threshold == pytest.approx(0.909090909090909)
    assert max_threshold == pytest.approx(1.1)


def test_check_length_ratios_1():
    length_matrix = {('A', 'B'): 1.1, ('B', 'A'): 0.909090909090909,
                     ('A', 'C'): 1.2, ('C', 'A'): 0.833333333333333,
                     ('B', 'C'): 1.090909090909091, ('C', 'B'): 0.916666666666667}
    trycycler.initial_check.check_length_ratios(length_matrix, 1.25)


def test_check_length_ratios_2():
    length_matrix = {('A', 'B'): 1.1, ('B', 'A'): 0.909090909090909,
                     ('A', 'C'): 1.2, ('C', 'A'): 0.833333333333333,
                     ('B', 'C'): 1.090909090909091, ('C', 'B'): 0.916666666666667}
    with pytest.raises(SystemExit):
        trycycler.initial_check.check_length_ratios(length_matrix, 1.15)


def test_check_mash_distances_1():
    distance_matrix = {('A', 'B'): 0.01, ('B', 'A'): 0.01,
                       ('A', 'C'): 0.02, ('C', 'A'): 0.02,
                       ('B', 'C'): 0.03, ('C', 'B'): 0.03}
    trycycler.initial_check.check_mash_distances(distance_matrix, 0.05)


def test_check_mash_distances_2():
    distance_matrix = {('A', 'B'): 0.01, ('B', 'A'): 0.01,
                       ('A', 'C'): 0.02, ('C', 'A'): 0.02,
                       ('B', 'C'): 0.03, ('C', 'B'): 0.03}
    with pytest.raises(SystemExit):
        trycycler.initial_check.check_mash_distances(distance_matrix, 0.005)


def test_get_length_ratio_matrix():
    seq_names = ['A', 'B', 'C']
    seqs = {'A': 'ACGACTCAGACTACGACTAC', 'B': 'ACGACTCAGA', 'C': 'ACGACTCAGACTACGA'}
    matrix = trycycler.initial_check.get_length_ratio_matrix(seq_names, seqs, 1.5)
    assert matrix[('A', 'B')] == pytest.approx(2.0)
    assert matrix[('B', 'A')] == pytest.approx(0.5)
    assert matrix[('A', 'C')] == pytest.approx(1.25)
    assert matrix[('C', 'A')] == pytest.approx(0.8)
    assert matrix[('B', 'C')] == pytest.approx(0.625)
    assert matrix[('C', 'B')] == pytest.approx(1.6)
