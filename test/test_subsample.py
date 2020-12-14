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

import pathlib
import pytest
import tempfile

import trycycler.subsample
import trycycler.misc


def test_welcome_message(capsys):
    trycycler.subsample.welcome_message()
    _, err = capsys.readouterr()
    assert 'Trycycler subsample' in err


def test_finished_message(capsys):
    trycycler.subsample.finished_message()
    _, err = capsys.readouterr()
    assert 'Finished' in err


def test_interpret_genome_size_01():
    assert trycycler.subsample.interpret_genome_size('1000') == 1000


def test_interpret_genome_size_02():
    assert trycycler.subsample.interpret_genome_size('283746283746') == 283746283746


def test_interpret_genome_size_03():
    assert trycycler.subsample.interpret_genome_size('5M') == 5000000


def test_interpret_genome_size_04():
    assert trycycler.subsample.interpret_genome_size('6.76m') == 6760000


def test_interpret_genome_size_05():
    assert trycycler.subsample.interpret_genome_size('0.23M') == 230000


def test_interpret_genome_size_06():
    assert trycycler.subsample.interpret_genome_size('35k') == 35000


def test_interpret_genome_size_07():
    assert trycycler.subsample.interpret_genome_size('4567.457K') == 4567457


def test_interpret_genome_size_08():
    assert trycycler.subsample.interpret_genome_size('2G') == 2000000000


def test_interpret_genome_size_09():
    assert trycycler.subsample.interpret_genome_size('1.45735g') == 1457350000


def test_interpret_genome_size_10():
    assert trycycler.subsample.interpret_genome_size('12323423.1') == 12323423


def test_interpret_genome_size_11():
    assert trycycler.subsample.interpret_genome_size('12323423.8') == 12323424


def test_interpret_genome_size_12():
    with pytest.raises(SystemExit) as e:
        trycycler.subsample.interpret_genome_size('abcdefg')
    assert 'cannot interpret genome size' in str(e.value)


def test_interpret_genome_size_13():
    with pytest.raises(SystemExit) as e:
        trycycler.subsample.interpret_genome_size('1.2.3')
    assert 'cannot interpret genome size' in str(e.value)


def test_interpret_genome_size_14():
    with pytest.raises(SystemExit) as e:
        trycycler.subsample.interpret_genome_size('123q')
    assert 'cannot interpret genome size' in str(e.value)


def test_interpret_genome_size_15():
    with pytest.raises(SystemExit) as e:
        trycycler.subsample.interpret_genome_size('123q')
    assert 'cannot interpret genome size' in str(e.value)


def test_interpret_genome_size_16():
    with pytest.raises(SystemExit) as e:
        trycycler.subsample.interpret_genome_size('23mk')
    assert 'cannot interpret genome size' in str(e.value)


def test_interpret_genome_size_17():
    with pytest.raises(SystemExit) as e:
        trycycler.subsample.interpret_genome_size('2o3m')
    assert 'cannot interpret genome size' in str(e.value)


def test_calculate_subsets_01():
    assert trycycler.subsample.calculate_subsets(1000, 25000000, 1000000, 25) == 1000


def test_calculate_subsets_02():
    assert trycycler.subsample.calculate_subsets(1000, 100000000, 1000000, 25) == 500


def test_calculate_subsets_03():
    assert trycycler.subsample.calculate_subsets(4000, 400000000, 1000000, 25) == 750


def test_calculate_subsets_04():
    with pytest.raises(SystemExit):
        trycycler.subsample.calculate_subsets(1000, 1000000, 1000000, 25)


def test_shuffle_reads():
    read_order = trycycler.subsample.shuffle_reads('test/test_subsample/reads_1.fastq')
    assert len(read_order) == 20
    assert sorted(read_order) == list(range(20))


def test_get_gfa_stats():
    contig_count, total_size, n50 = \
        trycycler.subsample.get_gfa_stats('test/test_subsample/graph.gfa', 0)
    assert contig_count == 20
    assert total_size == 2441
    assert n50 == 148


def test_save_subsets_01():
    output_count = 0
    read_names = set()
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = pathlib.Path(temp_dir)
        trycycler.subsample.save_subsets('test/test_subsample/reads_1.fastq', 12, 10, temp_dir)

        for f in temp_dir.glob('*.fastq'):
            output_count += 1
            assert len(list(trycycler.misc.iterate_fastq(f))) == 10
            for name, _, _, _ in trycycler.misc.iterate_fastq(f):
                read_names.add(name)

    assert output_count == 12
    assert len(read_names) == 20
