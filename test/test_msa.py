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
import pathlib
import pytest
import tempfile

import trycycler.msa
from trycycler.misc import load_fasta


def test_welcome_message(capsys):
    trycycler.msa.welcome_message()
    _, err = capsys.readouterr()
    assert 'Trycycler MSA' in err


def test_bad_input_file():
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = pathlib.Path(temp_dir)
        with open(temp_dir / '2_all_seqs.fasta', 'wt') as fasta:
            fasta.write(f'?bad_format\nACGATCGACATCAGT\n')
        with pytest.raises(SystemExit) as e:
            trycycler.msa.check_input_sequences(temp_dir)
        assert 'not in FASTA format' in str(e.value)


def test_only_one_input():
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = pathlib.Path(temp_dir)
        with open(temp_dir / '2_all_seqs.fasta', 'wt') as fasta:
            fasta.write(f'>only_seq\nACGATCGACATCAGT\n')
        with pytest.raises(SystemExit) as e:
            trycycler.msa.check_input_sequences(temp_dir)
        assert 'must have two or more sequences' in str(e.value)


def test_cluster_directory_not_exist():
    cluster_dir = pathlib.Path('not_a_real_dir')
    with pytest.raises(SystemExit) as e:
        trycycler.msa.check_cluster_directory(cluster_dir)
    assert 'does not exist' in str(e.value)


def test_cluster_directory_is_a_file():
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_file = pathlib.Path(temp_dir) / 'temp_file'
        temp_file.touch()
        with pytest.raises(SystemExit) as e:
            trycycler.msa.check_cluster_directory(temp_file)
        assert 'already exists as a file' in str(e.value)


def test_incomplete_muscle_1():
    # Tests missing files.
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = pathlib.Path(temp_dir)
        with pytest.raises(SystemExit) as e:
            trycycler.msa.check_muscle_results(temp_dir, 2)
        assert 'MUSCLE failed to complete' in str(e.value)


def test_incomplete_muscle_2():
    # Tests empty files.
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = pathlib.Path(temp_dir)
        file_0 = temp_dir / f'{0:012d}_msa.fasta'
        open(file_0, 'a').close()
        file_1 = temp_dir / f'{1:012d}_msa.fasta'
        open(file_1, 'a').close()
        with pytest.raises(SystemExit) as e:
            trycycler.msa.check_muscle_results(temp_dir, 2)
        assert 'MUSCLE failed to complete' in str(e.value)


def create_args(temp_dir, kmer, step, lookahead, threads):
    Args = collections.namedtuple('Args', ['cluster_dir', 'kmer', 'step', 'lookahead', 'threads'])
    return Args(cluster_dir=temp_dir, kmer=kmer, step=step, lookahead=lookahead, threads=threads)


def save_input_sequences(temp_dir, seqs):
    with open(temp_dir / '2_all_seqs.fasta', 'wt') as fasta:
        for name, seq in seqs:
            fasta.write(f'>{name}\n{seq}\n')


def test_msa_1():
    # Using default partitioning parameters, one thread.
    seqs = [('A', 'ATGTAAAGGTTCCGGGGCACTTAGCAGCTCCACAAATCCATTCCAACCTATA'),
            ('B', 'ATGTAAAGGTTGGGGCACTTAGCACTCCACACATCCATTGCCAACCTATA')]
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = pathlib.Path(temp_dir)
        save_input_sequences(temp_dir, seqs)
        trycycler.msa.msa(create_args(temp_dir, kmer=32, step=1000, lookahead=10000, threads=1))
        aligned_seqs = load_fasta(temp_dir / '3_msa.fasta')

    assert aligned_seqs == [('A', 'ATGTAAAGGTTCCGGGGCACTTAGCAGCTCCACAAATCCATT-CCAACCTATA'),
                            ('B', 'ATGTAAAGGTT--GGGGCACTTAGCA-CTCCACACATCCATTGCCAACCTATA')]


def test_msa_2():
    # Same sequences, but with multiple threads.
    seqs = [('A', 'ATGTAAAGGTTCCGGGGCACTTAGCAGCTCCACAAATCCATTCCAACCTATA'),
            ('B', 'ATGTAAAGGTTGGGGCACTTAGCACTCCACACATCCATTGCCAACCTATA')]
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = pathlib.Path(temp_dir)
        save_input_sequences(temp_dir, seqs)
        trycycler.msa.msa(create_args(temp_dir, kmer=32, step=1000, lookahead=10000, threads=8))
        aligned_seqs = load_fasta(temp_dir / '3_msa.fasta')

    assert aligned_seqs == [('A', 'ATGTAAAGGTTCCGGGGCACTTAGCAGCTCCACAAATCCATT-CCAACCTATA'),
                            ('B', 'ATGTAAAGGTT--GGGGCACTTAGCA-CTCCACACATCCATTGCCAACCTATA')]


def test_msa_3():
    # Same sequences, but with small partitioning parameters.
    seqs = [('A', 'ATGTAAAGGTTCCGGGGCACTTAGCAGCTCCACAAATCCATTCCAACCTATA'),
            ('B', 'ATGTAAAGGTTGGGGCACTTAGCACTCCACACATCCATTGCCAACCTATA')]
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = pathlib.Path(temp_dir)
        save_input_sequences(temp_dir, seqs)
        trycycler.msa.msa(create_args(temp_dir, kmer=8, step=12, lookahead=24, threads=8))
        aligned_seqs = load_fasta(temp_dir / '3_msa.fasta')

    print(aligned_seqs)

    assert aligned_seqs == [('A', 'ATGTAAAGGTTCCGGGGCACTTAGCAGCTCCACAAATCCATT-CCAACCTATA'),
                            ('B', 'ATGTAAAGGTT--GGGGCACTTAGCA-CTCCACACATCCATTGCCAACCTATA')]


def test_msa_4():
    # Same sequences, but with lower case bases (should be made uppercase in alignment).
    seqs = [('A', 'ATGtAAAGGTTcCGGGGCACttAGCaGCTCCACAaAtCcATTCCAACcTaTA'),
            ('B', 'ATGTaAaGGtTgGgGCAcTTAGCACTCCaCACATcCAttGCCaACCTATA')]
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = pathlib.Path(temp_dir)
        save_input_sequences(temp_dir, seqs)
        trycycler.msa.msa(create_args(temp_dir, kmer=32, step=1000, lookahead=10000, threads=8))
        aligned_seqs = load_fasta(temp_dir / '3_msa.fasta')

    assert aligned_seqs == [('A', 'ATGTAAAGGTTCCGGGGCACTTAGCAGCTCCACAAATCCATT-CCAACCTATA'),
                            ('B', 'ATGTAAAGGTT--GGGGCACTTAGCA-CTCCACACATCCATTGCCAACCTATA')]
