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

import trycycler.cluster


def test_input_assemblies_1():
    # Tests the function on two FASTAs without problems.
    total_lengths = trycycler.cluster.check_input_assemblies(['test/test_cluster/contigs_1.fasta',
                                                              'test/test_cluster/contigs_2.fasta'])
    assert total_lengths == {'A': 48, 'B': 24}


def test_input_assemblies_2():
    # Ensures that the function doesn't accept duplicate contig names.
    with pytest.raises(SystemExit) as e:
        trycycler.cluster.check_input_assemblies(['test/test_cluster/contigs_1.fasta',
                                                  'test/test_cluster/contigs_2.fasta',
                                                  'test/test_cluster/duplicate_contig_names.fasta'])
    assert 'duplicate' in str(e.value)


def test_input_assemblies_3():
    # Ensures that the function doesn't accept only one assembly.
    with pytest.raises(SystemExit) as e:
        trycycler.cluster.check_input_assemblies(['test/test_cluster/contigs_1.fasta'])
    assert 'two or more' in str(e.value)


def test_input_assemblies_4():
    # Ensures that the function doesn't accept non-FASTA files.
    with pytest.raises(SystemExit) as e:
        trycycler.cluster.check_input_assemblies(['test/test_cluster/contigs_1.fasta',
                                                  'test/test_cluster/not_a_fasta_file'])
    assert 'not in FASTA format' in str(e.value)


def test_input_assemblies_5():
    # Ensures that the function doesn't accept files which don't exist.
    with pytest.raises(SystemExit) as e:
        trycycler.cluster.check_input_assemblies(['test/test_cluster/contigs_1.fasta',
                                                  'test/test_cluster/not_a_file'])
    assert 'could not find' in str(e.value)


# TODO: make some tests with hashes in contig names: might cause problems with depths
