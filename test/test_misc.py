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

import gzip
import pathlib
import pytest
import sys
import tempfile
import unittest.mock

import trycycler.misc


def test_get_compression_type_1():
    assert trycycler.misc.get_compression_type('test/test_misc/test.txt') == 'plain'


def test_get_compression_type_2():
    assert trycycler.misc.get_compression_type('test/test_misc/test.gz') == 'gz'


def test_get_compression_type_3():
    with pytest.raises(SystemExit) as e:
        trycycler.misc.get_compression_type('test/test_misc/test.bz2')
    assert 'cannot use bzip2' in str(e.value)


def test_get_compression_type_4():
    with pytest.raises(SystemExit) as e:
        trycycler.misc.get_compression_type('test/test_misc/test.zip')
    assert 'cannot use zip' in str(e.value)


def test_get_open_func_1():
    assert trycycler.misc.get_open_func('test/test_misc/test.txt') == open


def test_get_open_func_2():
    assert trycycler.misc.get_open_func('test/test_misc/test.gz') == gzip.open


def test_get_sequence_file_type_1():
    assert trycycler.misc.get_sequence_file_type('test/test_misc/test.fasta') == 'FASTA'


def test_get_sequence_file_type_2():
    assert trycycler.misc.get_sequence_file_type('test/test_misc/test.fastq') == 'FASTQ'


def test_get_sequence_file_type_3():
    assert trycycler.misc.get_sequence_file_type('test/test_misc/test.txt') == 'neither'


def test_get_sequence_file_type_4():
    assert trycycler.misc.get_sequence_file_type('test/test_misc/empty') == 'neither'


def test_get_sequence_file_type_5():
    assert trycycler.misc.get_sequence_file_type('test/test_misc/test.fasta.gz') == 'FASTA'


def test_get_sequence_file_type_6():
    assert trycycler.misc.get_sequence_file_type('test/test_misc/test.fastq.gz') == 'FASTQ'


def test_get_sequence_file_type_7():
    assert trycycler.misc.get_sequence_file_type('test/test_misc/test.gz') == 'neither'


def test_get_sequence_file_type_8():
    assert trycycler.misc.get_sequence_file_type('test/test_misc/not_unicode') == 'neither'


def test_iterate_fastq_1():
    seqs = list(trycycler.misc.iterate_fastq('test/test_misc/test.fastq'))
    assert len(seqs) == 2
    assert seqs[0][0] == 'A'
    assert seqs[0][1] == '@A info'
    assert seqs[0][2].startswith('TTGCCTGTAGTCGGGACC')
    assert seqs[0][3].startswith('##$#%#%&++3*&&&-.7')
    assert seqs[1][0] == 'B'
    assert seqs[1][1] == '@B stuff'
    assert seqs[1][2].startswith('ATTCTCAGAATGGCGTAG')
    assert seqs[1][3].startswith(':;@@AHD98/.5C*-CEC')


def test_iterate_fastq_2():
    # Tests a FASTQ with extra line breaks.
    seqs = list(trycycler.misc.iterate_fastq('test/test_misc/bad_1.fastq'))
    assert len(seqs) == 2
    assert seqs[0][0] == 'A'
    assert seqs[0][1] == '@A info'
    assert seqs[0][2].startswith('TTGCCTGTAGTCGGGACC')
    assert seqs[0][3].startswith('##$#%#%&++3*&&&-.7')
    assert seqs[1][0] == 'B'
    assert seqs[1][1] == '@B stuff'
    assert seqs[1][2].startswith('ATTCTCAGAATGGCGTAG')
    assert seqs[1][3].startswith(':;@@AHD98/.5C*-CEC')


def test_iterate_fastq_3():
    # Tests a FASTQ with an extra line of text.
    seqs = list(trycycler.misc.iterate_fastq('test/test_misc/bad_2.fastq'))
    assert len(seqs) == 2
    assert seqs[0][0] == 'A'
    assert seqs[0][1] == '@A info'
    assert seqs[0][2].startswith('TTGCCTGTAGTCGGGACC')
    assert seqs[0][3].startswith('##$#%#%&++3*&&&-.7')
    assert seqs[1][0] == 'B'
    assert seqs[1][1] == '@B stuff'
    assert seqs[1][2].startswith('ATTCTCAGAATGGCGTAG')
    assert seqs[1][3].startswith(':;@@AHD98/.5C*-CEC')


def test_iterate_fastq_4():
    with pytest.raises(SystemExit) as e:
        _ = list(trycycler.misc.iterate_fastq('test/test_misc/test.fasta'))
    assert 'not FASTQ format' in str(e.value)


def test_load_fastq_as_dict():
    seqs = trycycler.misc.load_fastq_as_dict('test/test_misc/test.fastq')
    assert len(seqs) == 2
    assert seqs['A'][0] == '@A info'
    assert seqs['A'][1].startswith('TTGCCTGTAGTCGGGACC')
    assert seqs['A'][2].startswith('##$#%#%&++3*&&&-.7')
    assert seqs['B'][0] == '@B stuff'
    assert seqs['B'][1].startswith('ATTCTCAGAATGGCGTAG')
    assert seqs['B'][2].startswith(':;@@AHD98/.5C*-CEC')


def test_get_fastq_stats():
    read_count, total_size, n50 = trycycler.misc.get_fastq_stats('test/test_misc/test.fastq')
    assert read_count == 2
    assert total_size == 200
    assert n50 == 100


def test_get_n50_1():
    assert trycycler.misc.get_n50([1, 2, 3, 4, 1000]) == 1000


def test_get_n50_2():
    assert trycycler.misc.get_n50([12, 23455, 15, 12433, 15343, 9, 10]) == 15343


def test_get_n50_3():
    assert trycycler.misc.get_n50([]) == 0


def test_load_fasta_1():
    seqs = trycycler.misc.load_fasta('test/test_misc/test.fasta')
    assert len(seqs) == 2
    assert seqs[0][0] == 'A'
    assert seqs[0][1].startswith('TTGCCTGTAGTCGGGACC')
    assert seqs[1][0] == 'B'
    assert seqs[1][1].startswith('ATTCTCAGAATGGCGTAG')


def test_load_fasta_2():
    seqs = trycycler.misc.load_fasta('test/test_misc/test.fasta', include_full_header=True)
    assert len(seqs) == 2
    assert seqs[0][0] == 'A'
    assert seqs[0][1] == 'A info'
    assert seqs[0][2].startswith('TTGCCTGTAGTCGGGACC')
    assert seqs[1][0] == 'B'
    assert seqs[1][1] == 'B stuff'
    assert seqs[1][2].startswith('ATTCTCAGAATGGCGTAG')


def test_load_fasta_3():
    seqs = trycycler.misc.load_fasta('test/test_misc/test.fasta.gz')
    assert len(seqs) == 2
    assert seqs[0][0] == 'A'
    assert seqs[0][1].startswith('TTGCCTGTAGTCGGGACC')
    assert seqs[1][0] == 'B'
    assert seqs[1][1].startswith('ATTCTCAGAATGGCGTAG')


def test_load_fasta_4():
    seqs = trycycler.misc.load_fasta('test/test_misc/bad_1.fasta')
    assert len(seqs) == 2
    assert seqs[0][0] == 'A'
    assert seqs[0][1].startswith('TTGCCTGTAGTCGGGACC')
    assert seqs[1][0] == 'B'
    assert seqs[1][1].startswith('ATTCTCAGAATGGCGTAG')


def test_load_fasta_5():
    seqs = trycycler.misc.load_fasta('test/test_misc/lowercase.fasta')
    assert len(seqs) == 2
    assert seqs[0][0] == 'A'
    assert seqs[0][1].startswith('TTGCCTGTAGTCGGGACC')
    assert seqs[1][0] == 'B'
    assert seqs[1][1].startswith('ATTCTCAGAATGGCGTAG')


def test_get_default_thread_count():
    assert 1 <= trycycler.misc.get_default_thread_count() <= 16


def test_write_seq_to_fasta_1():
    with tempfile.TemporaryDirectory() as temp_dir:
        filename = pathlib.Path(temp_dir) / 'temp.fasta'
        trycycler.misc.write_seq_to_fasta('CAGAATGGCGT', 'name', filename)
        seqs = trycycler.misc.load_fasta(filename)
        assert len(seqs) == 1
        assert seqs[0][0] == 'name'
        assert seqs[0][1] == 'CAGAATGGCGT'


def test_write_seq_to_fasta_2():
    # Same test, but with lowercase bases in input (should be made uppercase in saved file).
    with tempfile.TemporaryDirectory() as temp_dir:
        filename = pathlib.Path(temp_dir) / 'temp.fasta'
        trycycler.misc.write_seq_to_fasta('CAgaaTgGcgt', 'name', filename)
        seqs = trycycler.misc.load_fasta(filename)
        assert len(seqs) == 1
        assert seqs[0][0] == 'name'
        assert seqs[0][1] == 'CAGAATGGCGT'


def test_reverse_complement_1():
    assert trycycler.misc.reverse_complement('GGGGaaaaaaaatttatatat') == 'atatataaattttttttCCCC'


def test_reverse_complement_2():
    assert trycycler.misc.reverse_complement('atatataaattttttttCCCC') == 'GGGGaaaaaaaatttatatat'


def test_reverse_complement_3():
    assert trycycler.misc.reverse_complement('ACGT123') == 'NNNACGT'


def test_remove_duplicates_1():
    assert trycycler.misc.remove_duplicates([1, 4, 3, 4, 2]) == [1, 4, 3, 2]


def test_remove_duplicates_2():
    assert trycycler.misc.remove_duplicates(['a', 'a', 'a', 'b', 'a']) == ['a', 'b']


def test_check_python_version_1():
    with unittest.mock.patch.object(sys, 'version_info') as v_info:
        v_info.major = 3
        v_info.minor = 6
        trycycler.misc.check_python_version()


def test_check_python_version_2():
    with unittest.mock.patch.object(sys, 'version_info') as v_info:
        v_info.major = 3
        v_info.minor = 8
        trycycler.misc.check_python_version()


def test_check_python_version_3():
    with pytest.raises(SystemExit) as e:
        with unittest.mock.patch.object(sys, 'version_info') as v_info:
            v_info.major = 3
            v_info.minor = 5
            trycycler.misc.check_python_version()
    assert 'requires Python 3.6 or later' in str(e.value)


def test_check_python_version_4():
    with pytest.raises(SystemExit) as e:
        with unittest.mock.patch.object(sys, 'version_info') as v_info:
            v_info.major = 2
            v_info.minor = 7
            trycycler.misc.check_python_version()
    assert 'requires Python 3.6 or later' in str(e.value)


def test_check_output_directory_1():
    with pytest.raises(SystemExit) as e:
        trycycler.misc.check_output_directory(pathlib.Path('test/test_misc/test.fasta'))
    assert 'already exists as a file' in str(e.value)


def test_check_output_directory_2():
    with tempfile.TemporaryDirectory() as temp_dir:
        out_dir = pathlib.Path(temp_dir) / 'output'
        trycycler.misc.check_output_directory(out_dir)
        assert out_dir.is_dir()
        trycycler.misc.check_output_directory(out_dir)
        assert out_dir.is_dir()
        temp_file = out_dir / 'temp'
        open(temp_file, 'a').close()
        trycycler.misc.check_output_directory(out_dir)
        assert out_dir.is_dir()


def test_count_substrings_1():
    assert trycycler.misc.count_substrings('000123000123', '123') == 2


def test_count_substrings_2():
    assert trycycler.misc.count_substrings('000123000123', 'abc') == 0


def test_range_overlap_1():
    assert trycycler.misc.range_overlap(0, 10, 5, 20)


def test_range_overlap_2():
    assert trycycler.misc.range_overlap(0, 10, 9, 20)


def test_range_overlap_3():
    assert not trycycler.misc.range_overlap(0, 10, 10, 20)


def test_range_overlap_4():
    assert not trycycler.misc.range_overlap(0, 10, 11, 20)


def test_check_input_reads_1():
    read_count, total_size = trycycler.misc.check_input_reads('test/test_misc/test.fastq.gz')
    assert read_count == 2
    assert total_size == 200


def test_check_input_reads_2():
    file_size = trycycler.misc.check_input_reads('test/test_misc/test.fastq.gz',
                                                 file_size_only=True)
    assert file_size > 100


def test_check_input_reads_3():
    with pytest.raises(SystemExit) as e:
        trycycler.misc.check_input_reads('test/test_misc/test.fasta')
    assert 'not in FASTQ format' in str(e.value)


def test_get_ascii_art():
    assert "| || '__|| | | | / __|| | | |" in trycycler.misc.get_ascii_art()


def test_count_lines_1():
    assert trycycler.misc.count_lines('test/test_misc/test.fasta') == 4


def test_count_lines_2():
    assert trycycler.misc.count_lines('test/test_misc/test.fastq.gz') == 8
