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


def test_partition_1():
    msa_names = ['1', '2', '3']
    msa_seqs = {'1': 'AAAAACAAAAAA',
                '2': 'AAAAAGAAAAAA',
                '3': 'AAAAATAAAAAA'}
    msa_length = 12
    chunks = trycycler.consensus.partition_msa(msa_seqs, msa_names, msa_length, 1)

    assert len(chunks) == 3

    assert chunks[0].type == 'same'
    assert chunks[0].seq == ['A', 'A', 'A', 'A', 'A']

    assert chunks[1].type == 'different'
    assert chunks[1].seqs == {'1': ['C'], '2': ['G'], '3': ['T']}

    assert chunks[2].type == 'same'
    assert chunks[2].seq == ['A', 'A', 'A', 'A', 'A', 'A']


def test_partition_2():
    msa_names = ['1', '2', '3']
    msa_seqs = {'1': 'AAAAACAACAAAA',
                '2': 'AAAAAGAAGAAAA',
                '3': 'AAAAATAATAAAA'}
    msa_length = 13
    chunks = trycycler.consensus.partition_msa(msa_seqs, msa_names, msa_length, 1)

    assert len(chunks) == 5

    assert chunks[0].type == 'same'
    assert chunks[0].seq == ['A', 'A', 'A', 'A', 'A']

    assert chunks[1].type == 'different'
    assert chunks[1].seqs == {'1': ['C'], '2': ['G'], '3': ['T']}

    assert chunks[2].type == 'same'
    assert chunks[2].seq == ['A', 'A']

    assert chunks[3].type == 'different'
    assert chunks[3].seqs == {'1': ['C'], '2': ['G'], '3': ['T']}

    assert chunks[4].type == 'same'
    assert chunks[4].seq == ['A', 'A', 'A', 'A']


def test_partition_3():
    msa_names = ['1', '2', '3']
    msa_seqs = {'1': 'AAAAACAACAAAA',
                '2': 'AAAAAGAAGAAAA',
                '3': 'AAAAATAATAAAA'}
    msa_length = 13
    chunks = trycycler.consensus.partition_msa(msa_seqs, msa_names, msa_length, 2)

    assert len(chunks) == 3

    assert chunks[0].type == 'same'
    assert chunks[0].seq == ['A', 'A', 'A', 'A', 'A']

    assert chunks[1].type == 'different'
    assert chunks[1].seqs == {'1': ['C', 'A', 'A', 'C'],
                              '2': ['G', 'A', 'A', 'G'],
                              '3': ['T', 'A', 'A', 'T']}

    assert chunks[2].type == 'same'
    assert chunks[2].seq == ['A', 'A', 'A', 'A']


def test_get_best_seq_1():
    """
    Tests get_best_seq() on a 'same' chunk: should just return the sequence.
    """
    c = trycycler.consensus.Chunk()
    c.add_bases({'1': 'A', '2': 'A', '3': 'A'})
    c.add_bases({'1': 'A', '2': 'A', '3': 'A'})
    c.add_bases({'1': 'A', '2': 'A', '3': 'A'})
    assert c.get_best_seq() == 'AAA'


def test_get_best_seq_2():
    """
    Tests get_best_seq() on simple 'different' chunk with a clear winner:
    AAA, AAA and CCC
    """
    c = trycycler.consensus.Chunk()
    c.add_bases({'1': 'A', '2': 'A', '3': 'C'})
    c.add_bases({'1': 'A', '2': 'A', '3': 'C'})
    c.add_bases({'1': 'A', '2': 'A', '3': 'C'})
    assert c.get_best_seq() == 'AAA'


def test_get_best_seq_3():
    """
    Tests get_best_seq() on 'different' chunk with a tie that can be broken by Hamming distance:
    AAA, AAA, CCC, CCC, CCT
    """
    c = trycycler.consensus.Chunk()
    c.add_bases({'1': 'A', '2': 'A', '3': 'C', '4': 'C', '5': 'C'})
    c.add_bases({'1': 'A', '2': 'A', '3': 'C', '4': 'C', '5': 'C'})
    c.add_bases({'1': 'A', '2': 'A', '3': 'C', '4': 'C', '5': 'T'})
    assert c.get_best_seq() == 'CCC'


def test_get_best_seq_4():
    """
    Tests get_best_seq() on 'different' chunk with a tie that cannot be broken by Hamming distance:
    AAA, AAA, CCC, CCC, TTT
    """
    c = trycycler.consensus.Chunk()
    c.add_bases({'1': 'A', '2': 'A', '3': 'C', '4': 'C', '5': 'T'})
    c.add_bases({'1': 'A', '2': 'A', '3': 'C', '4': 'C', '5': 'T'})
    c.add_bases({'1': 'A', '2': 'A', '3': 'C', '4': 'C', '5': 'T'})
    assert c.get_best_seq() == 'AAA'
