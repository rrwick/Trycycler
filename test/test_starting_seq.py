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

import trycycler.starting_seq
import trycycler.misc


def test_look_for_known_starting_seq():
    seqs = {'test': trycycler.misc.load_fasta('test/test_starting_seq/test.fasta')[0][1]}
    starting_genes = trycycler.misc.load_fasta('trycycler/data/starting_genes.fasta')
    starting_seq = trycycler.starting_seq.look_for_known_starting_seq(seqs, 1)
    assert starting_seq == starting_genes[0][1]


def test_flip_seqs_as_necessary():
    seqs = {'test': trycycler.misc.load_fasta('test/test_starting_seq/test.fasta')[0][1]}
    starting_seq = trycycler.misc.load_fasta('trycycler/data/starting_genes.fasta')[0][1]
    strand_fixed_seqs = trycycler.starting_seq.flip_seqs_as_necessary(seqs, starting_seq)
    assert strand_fixed_seqs['test'] == seqs['test']

