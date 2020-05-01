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

import pytest

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
    c.set_best_seq_as_most_common()
    assert c.best_seq == 'AAA'
    assert not c.had_tie


def test_get_best_seq_2():
    """
    Tests get_best_seq() on simple 'different' chunk with a clear winner:
    AAA, AAA and CCC
    """
    c = trycycler.consensus.Chunk()
    c.add_bases({'1': 'A', '2': 'A', '3': 'C'})
    c.add_bases({'1': 'A', '2': 'A', '3': 'C'})
    c.add_bases({'1': 'A', '2': 'A', '3': 'C'})
    c.set_best_seq_as_most_common()
    assert c.best_seq == 'AAA'
    assert not c.had_tie


def test_get_best_seq_3():
    """
    Tests get_best_seq() on 'different' chunk with a tie that can be broken by Hamming distance:
    AAA, AAA, CCC, CCC, CCT
    """
    c = trycycler.consensus.Chunk()
    c.add_bases({'1': 'A', '2': 'A', '3': 'C', '4': 'C', '5': 'C'})
    c.add_bases({'1': 'A', '2': 'A', '3': 'C', '4': 'C', '5': 'C'})
    c.add_bases({'1': 'A', '2': 'A', '3': 'C', '4': 'C', '5': 'T'})
    c.set_best_seq_as_most_common()
    assert c.best_seq == 'CCC'
    assert not c.had_tie


def test_get_best_seq_4():
    """
    Tests get_best_seq() on 'different' chunk with a tie that cannot be broken by Hamming distance:
    AAA, AAA, CCC, CCC, TTT
    """
    c = trycycler.consensus.Chunk()
    c.add_bases({'1': 'A', '2': 'A', '3': 'C', '4': 'C', '5': 'T'})
    c.add_bases({'1': 'A', '2': 'A', '3': 'C', '4': 'C', '5': 'T'})
    c.add_bases({'1': 'A', '2': 'A', '3': 'C', '4': 'C', '5': 'T'})
    c.set_best_seq_as_most_common()
    assert c.best_seq == 'AAA'
    assert c.had_tie


def test_get_best_seq_5():
    """
    Tests get_best_seq() on an empty chunk.
    """
    c = trycycler.consensus.Chunk()
    c.set_best_seq_as_most_common()
    assert c.best_seq == ''


def test_make_ungapped_pos_to_gapped_pos_dict_1():
    with_gaps = 'AAA'
    without_gaps = 'AAA'
    ungapped_to_gapped = trycycler.consensus.make_ungapped_pos_to_gapped_pos_dict(with_gaps,
                                                                                  without_gaps)
    assert ungapped_to_gapped == {0: 0, 1: 1, 2: 2, 3: 3}


def test_make_ungapped_pos_to_gapped_pos_dict_2():
    with_gaps = 'A-A'
    without_gaps = 'AA'
    ungapped_to_gapped = trycycler.consensus.make_ungapped_pos_to_gapped_pos_dict(with_gaps,
                                                                                  without_gaps)
    assert ungapped_to_gapped == {0: 0, 1: 2, 2: 3}


def test_make_ungapped_pos_to_gapped_pos_dict_3():
    with_gaps = '-CGA-GAC--A-'
    without_gaps = 'CGAGACA'
    ungapped_to_gapped = trycycler.consensus.make_ungapped_pos_to_gapped_pos_dict(with_gaps,
                                                                                  without_gaps)
    assert ungapped_to_gapped == {0: 1, 1: 2, 2: 3, 3: 5, 4: 6, 5: 7, 6: 10, 7: 12}


def make_chunks_for_test_sequence():
    msa_names = ['1', '2', '3']
    msa_seqs = {'1': 'ACGACGAGCAC-GCAGACGACACGCGAACTAGCGCA-CATCGC',
                '2': 'A-GACGAGCACTGCA-ACGACACGCGAACTAGCGCAGCATCGC',
                '3': 'ACGACGAGCACTGCAGACGACACGCG-ACTAGCGCAGCATCGC'}
    msa_length = 43
    chunks = trycycler.consensus.partition_msa(msa_seqs, msa_names, msa_length, 0)
    assert len(chunks) == 11
    for chunk in chunks:
        chunk.set_best_seq_as_most_common()
    assert [c.best_seq for c in chunks] == ['A', 'C', 'GACGAGCAC', 'T', 'GCA', 'G', 'ACGACACGCG',
                                            'A', 'ACTAGCGCA', 'G', 'CATCGC']
    return chunks


def test_make_initial_consensus():
    chunks = make_chunks_for_test_sequence()
    with_gaps, without_gaps = trycycler.consensus.make_initial_consensus(chunks)
    assert with_gaps == 'ACGACGAGCACTGCAGACGACACGCGAACTAGCGCAGCATCGC'
    assert without_gaps == 'ACGACGAGCACTGCAGACGACACGCGAACTAGCGCAGCATCGC'


def test_build_test_sequence_linear_1():
    chunks = make_chunks_for_test_sequence()
    test_seq = trycycler.consensus.build_test_sequence(1, chunks, 'C', False, 1)
    assert test_seq == 'ACG'


def test_build_test_sequence_linear_2():
    chunks = make_chunks_for_test_sequence()
    test_seq = trycycler.consensus.build_test_sequence(1, chunks, '-', False, 1)
    assert test_seq == 'AG'


def test_build_test_sequence_linear_3():
    chunks = make_chunks_for_test_sequence()
    test_seq = trycycler.consensus.build_test_sequence(1, chunks, 'C', False, 2)
    assert test_seq == 'ACGA'


def test_build_test_sequence_linear_4():
    chunks = make_chunks_for_test_sequence()
    test_seq = trycycler.consensus.build_test_sequence(1, chunks, 'C', False, 3)
    assert test_seq == 'ACGAC'


def test_build_test_sequence_linear_5():
    # Using a very large margin on a linear sequence should return the entire sequence.
    chunks = make_chunks_for_test_sequence()
    test_seq = trycycler.consensus.build_test_sequence(1, chunks, 'C', False, 500)
    assert test_seq == 'ACGACGAGCACTGCAGACGACACGCGAACTAGCGCAGCATCGC'


def test_build_test_sequence_linear_6():
    chunks = make_chunks_for_test_sequence()
    test_seq = trycycler.consensus.build_test_sequence(9, chunks, 'G', False, 0)
    assert test_seq == 'G'


def test_build_test_sequence_linear_7():
    chunks = make_chunks_for_test_sequence()
    test_seq = trycycler.consensus.build_test_sequence(9, chunks, 'G', False, 1)
    assert test_seq == 'AGC'


def test_build_test_sequence_linear_8():
    chunks = make_chunks_for_test_sequence()
    test_seq = trycycler.consensus.build_test_sequence(9, chunks, 'G', False, 5)
    assert test_seq == 'GCGCAGCATCG'


def test_build_test_sequence_linear_9():
    chunks = make_chunks_for_test_sequence()
    test_seq = trycycler.consensus.build_test_sequence(9, chunks, 'G', False, 13)
    assert test_seq == 'GCGAACTAGCGCAGCATCGC'


def test_build_test_sequence_circular_1():
    chunks = make_chunks_for_test_sequence()
    test_seq = trycycler.consensus.build_test_sequence(1, chunks, 'C', True, 1)
    assert test_seq == 'ACG'


def test_build_test_sequence_circular_2():
    chunks = make_chunks_for_test_sequence()
    test_seq = trycycler.consensus.build_test_sequence(1, chunks, 'C', True, 3)
    assert test_seq == 'GCACGAC'


def test_build_test_sequence_circular_3():
    chunks = make_chunks_for_test_sequence()
    test_seq = trycycler.consensus.build_test_sequence(1, chunks, 'C', True, 12)
    assert test_seq == 'CGCAGCATCGCACGACGAGCACTGC'


def test_build_test_sequence_circular_4():
    chunks = make_chunks_for_test_sequence()
    test_seq = trycycler.consensus.build_test_sequence(1, chunks, 'C', True, 500)
    assert test_seq == 'GCGAACTAGCGCAGCATCGCACGACGAGCACTGCAGACGACAC'


def test_build_test_sequence_circular_5():
    chunks = make_chunks_for_test_sequence()
    test_seq = trycycler.consensus.build_test_sequence(9, chunks, '-', True, 2)
    assert test_seq == 'CACA'


def test_build_test_sequence_circular_6():
    chunks = make_chunks_for_test_sequence()
    test_seq = trycycler.consensus.build_test_sequence(9, chunks, '-', True, 1000)
    assert test_seq == 'GACGACACGCGAACTAGCGCACATCGCACGACGAGCACTGCA'


def test_hamming_distance_1():
    assert trycycler.consensus.hamming_distance('AAAA', 'AAAA') == 0


def test_hamming_distance_2():
    assert trycycler.consensus.hamming_distance('AAAA', 'AAAC') == 1


def test_hamming_distance_3():
    assert trycycler.consensus.hamming_distance('AAAA', 'CACC') == 3


def test_hamming_distance_4():
    # Strings of different length don't work.
    with pytest.raises(AssertionError):
        trycycler.consensus.hamming_distance('AAAA', 'AAA')


def make_different_chunk():
    c = trycycler.consensus.Chunk()
    c.add_bases({'1': 'A', '2': 'C'})
    c.add_bases({'1': 'A', '2': '-'})
    c.add_bases({'1': 'A', '2': '-'})
    c.add_bases({'1': 'A', '2': '-'})
    c.add_bases({'1': 'A', '2': 'C'})
    assert c.type == 'different'
    return c


def make_same_chunk():
    c = trycycler.consensus.Chunk()
    c.add_bases({'1': 'A'})
    c.add_bases({'1': 'A'})
    c.add_bases({'1': 'A'})
    c.add_bases({'1': 'A'})
    c.add_bases({'1': 'A'})
    assert c.type == 'same'
    return c


def make_none_chunk():
    c = trycycler.consensus.Chunk()
    assert c.type is None
    return c


def make_bogus_chunk():
    c = trycycler.consensus.Chunk()
    c.type = 'not_a_type'
    return c


def test_has_long_indel_1():
    c = make_different_chunk()
    assert c.has_long_indel(1)
    assert c.has_long_indel(2)
    assert c.has_long_indel(3)
    assert not c.has_long_indel(4)
    assert not c.has_long_indel(5)


def test_has_long_indel_2():
    c = make_same_chunk()
    assert not c.has_long_indel(1)
    assert not c.has_long_indel(5)


def test_has_long_indel_3():
    c = make_none_chunk()
    assert not c.has_long_indel(1)
    assert not c.has_long_indel(5)


def test_get_length_1():
    c = make_none_chunk()
    assert c.get_length() == 0


def test_get_length_2():
    c = make_same_chunk()
    assert c.get_length() == 5


def test_get_length_3():
    c = make_different_chunk()
    assert c.get_length() == 5


def test_get_length_4():
    c = make_bogus_chunk()
    with pytest.raises(AssertionError):
        c.get_length()
