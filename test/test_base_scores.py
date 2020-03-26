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

import trycycler.base_scores


def test_add_indels_to_per_base_scores_1():
    """
    seq:    ACGACTAGCTACACG  ->  ACGACT--AGCTAC-ACG
    scores: 356289238951365  ->  356289222389511365
    """
    msa_seqs = {'seq': 'ACGACT--AGCTAC-ACG'}
    per_base_scores = {'seq': [3, 5, 6, 2, 8, 9, 2, 3, 8, 9, 5, 1, 3, 6, 5]}
    result = trycycler.base_scores.add_indels_to_per_base_scores(msa_seqs, per_base_scores)
    assert result['seq'] == [3, 5, 6, 2, 8, 9, 2, 2, 2, 3, 8, 9, 5, 1, 1, 3, 6, 5]


def test_add_indels_to_per_base_scores_2():
    """
    seq:    CGT  ->  ---C---G---T---
    scores: 341  ->  333333341111111
    """
    msa_seqs = {'seq': '---C---G---T---'}
    per_base_scores = {'seq': [3, 4, 1]}
    result = trycycler.base_scores.add_indels_to_per_base_scores(msa_seqs, per_base_scores)
    assert result['seq'] == [3, 3, 3, 3, 3, 3, 3, 4, 1, 1, 1, 1, 1, 1, 1]


def test_add_indels_to_per_base_scores_3():
    """
    seq:    AC  ->  -----A-C-----
    scores: 12  ->  1111111222222
    """
    msa_seqs = {'seq': '-----A-C-----'}
    per_base_scores = {'seq': [1, 2]}
    result = trycycler.base_scores.add_indels_to_per_base_scores(msa_seqs, per_base_scores)
    assert result['seq'] == [1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2]


def test_add_indels_to_per_base_scores_4():
    """
    seq:    ACGTCG  ->  ACGTCG
    scores: 123343  ->  123343
    """
    msa_seqs = {'seq': 'ACGTCG'}
    per_base_scores = {'seq': [1, 2, 3, 3, 4, 3]}
    result = trycycler.base_scores.add_indels_to_per_base_scores(msa_seqs, per_base_scores)
    assert result['seq'] == [1, 2, 3, 3, 4, 3]
