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

import collections

import trycycler.dotplot


def test_welcome_message(capsys):
    trycycler.dotplot.welcome_message()
    _, err = capsys.readouterr()
    assert 'Trycycler dotplot' in err


def test_finished_message(capsys):
    trycycler.dotplot.finished_message()
    _, err = capsys.readouterr()
    assert 'Finished' in err


def test_get_all_kmer_positions():
    seq = 'ACGATCGAC'
    forward_kmers, reverse_kmers = trycycler.dotplot.get_all_kmer_positions(3, seq)

    assert len(forward_kmers) == 6
    assert len(reverse_kmers) == 6

    assert forward_kmers['ACG'] == [0]
    assert forward_kmers['CGA'] == [1, 5]
    assert forward_kmers['GAT'] == [2]
    assert forward_kmers['ATC'] == [3]
    assert forward_kmers['TCG'] == [4]
    assert forward_kmers['GAC'] == [6]

    assert reverse_kmers['GTC'] == [6]
    assert reverse_kmers['TCG'] == [5, 1]
    assert reverse_kmers['CGA'] == [4]
    assert reverse_kmers['GAT'] == [3]
    assert reverse_kmers['ATC'] == [2]
    assert reverse_kmers['CGT'] == [0]


def test_create_dotplots_1():
    Args = collections.namedtuple('Args', ['kmer', 'res'])
    args = Args(kmer=8, res=250)

    seqs = {'A': 'CCCAGCCGATCGCAGCTTGGACGGATTGTCGCAGGGTCGTTGCTTTGCCCTGCGACAATCCGTCAGATGGACTTA',
            'B': 'CCCAGCGGATCGCAGCTTGGACGGATTGTCGCAGGATCGTTGCTTTGCCCTGCGACAATCCGTCAGATGGACTTA'}
    seq_names = ['A', 'B']
    image = trycycler.dotplot.create_dotplots(seq_names, seqs, args)
    assert image.size == (250, 250)


def test_create_dotplots_2():
    Args = collections.namedtuple('Args', ['kmer', 'res'])
    args = Args(kmer=8, res=250)

    seqs = {'This_is_the_first_sequence_and_it_has_a_very_long_name':
            'CCCAGCCGATCGCAGCTTGGACGGATTGTCGCAGGGTCGTTGCTTTGCCCTGCGACAATCCGTCAGATGGACTTA',
            'This_is_the_second_sequence_and_it_also_has_a_very_long_name':
            'CCCAGCGGATCGCAGCTTGGACGGATTGTCGCAGGATCGTTGCTTTGCCCTGCGACAATCCGTCAGATGGACTTA'}
    seq_names = ['This_is_the_first_sequence_and_it_has_a_very_long_name',
                 'This_is_the_second_sequence_and_it_also_has_a_very_long_name']
    image = trycycler.dotplot.create_dotplots(seq_names, seqs, args)
    assert image.size == (250, 250)
