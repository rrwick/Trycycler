#!/usr/bin/env python3
"""
This script trims circular contigs in a Canu assembly based on the 'suggestCircular' and 'trim'
values in the FASTA header. It takes one argument (a Canu assembly FASTA filename) and it outputs
(to stdout) the same assembly where circular contigs are trimmed to have no overlap.

Copyright 2022 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Trycycler

This file is part of Trycycler. Trycycler is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Trycycler is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Trycycler.
If not, see <https://www.gnu.org/licenses/>.
"""

import argparse
import gzip
import re
import sys


def get_arguments():
    parser = argparse.ArgumentParser(description='Canu circular contig trimmer')
    parser.add_argument('input', type=str,
                        help='Filename of Canu assembly in FASTA format')
    args = parser.parse_args()
    return args


def main():
    args = get_arguments()
    assembly = load_fasta(args.input)
    for header, seq in assembly:
        header, seq = trim_seq(header, seq)
        print(f'>{header}')
        print(f'{seq}')


def trim_seq(header, seq):
    if 'suggestCircular=yes' in header:
        result = re.search(r'trim=(\d+)-(\d+)', header)
        if result is not None:
            trim_string = result.group(0)
            start = int(result.group(1))
            end = int(result.group(2))
            seq = seq[start:end]
            header = header.replace(trim_string, f'trim=0-{len(seq)}')
            header = re.sub(r'len=\d+', f'len={len(seq)}', header)
    return header, seq


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    https://stackoverflow.com/questions/13044562
    """
    magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
                  'bz2': (b'\x42', b'\x5a', b'\x68'),
                  'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
    max_len = max(len(x) for x in magic_dict)
    unknown_file = open(str(filename), 'rb')
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = 'plain'
    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = file_type
    if compression_type == 'bz2':
        sys.exit('\nError: cannot use bzip2 format - use gzip instead')
    if compression_type == 'zip':
        sys.exit('\nError: cannot use zip format - use gzip instead')
    return compression_type


def get_open_func(filename):
    if get_compression_type(filename) == 'gz':
        return gzip.open
    else:  # plain text
        return open


def load_fasta(fasta_filename):
    fasta_seqs = []
    with get_open_func(fasta_filename)(fasta_filename, 'rt') as fasta_file:
        name = ''
        sequence = []
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':  # Header line = start of new contig
                if name:
                    fasta_seqs.append((name, ''.join(sequence)))
                    sequence = []
                name = line[1:]
            else:
                sequence.append(line.upper())
        if name:
            fasta_seqs.append((name, ''.join(sequence)))
    return fasta_seqs


if __name__ == '__main__':
    main()


# Unit tests for Pytest
# =====================

def test_trim_seq_1():
    header = '>tig00000001 len=60 reads=50 class=contig suggestRepeat=no suggestBubble=no ' \
             'suggestCircular=no trim=0-60'
    seq = 'AGTAGCCAAACTATTTAATGCTAGAGATGCTGCATATCAAAAAATAATCAAACAATTATC'
    new_header, new_seq = trim_seq(header, seq)
    assert new_header == header
    assert new_seq == seq


def test_trim_seq_2():
    header = '>tig00000001 len=60 reads=50 class=contig suggestRepeat=no suggestBubble=no ' \
             'suggestCircular=yes trim=0-50'
    seq = 'AGTAGCCAAACTATTTAATGCTAGAGATGCTGCATATCAAAAAATAATCAAACAATTATC'
    new_header, new_seq = trim_seq(header, seq)
    assert new_header == '>tig00000001 len=50 reads=50 class=contig suggestRepeat=no ' \
                         'suggestBubble=no suggestCircular=yes trim=0-50'
    assert new_seq == 'AGTAGCCAAACTATTTAATGCTAGAGATGCTGCATATCAAAAAATAATCA'


def test_trim_seq_3():
    header = '>tig00000001 len=60 reads=50 class=contig suggestRepeat=no suggestBubble=no ' \
             'suggestCircular=yes trim=10-60'
    seq = 'AGTAGCCAAACTATTTAATGCTAGAGATGCTGCATATCAAAAAATAATCAAACAATTATC'
    new_header, new_seq = trim_seq(header, seq)
    assert new_header == '>tig00000001 len=50 reads=50 class=contig suggestRepeat=no ' \
                         'suggestBubble=no suggestCircular=yes trim=0-50'
    assert new_seq == 'CTATTTAATGCTAGAGATGCTGCATATCAAAAAATAATCAAACAATTATC'


def test_trim_seq_4():
    header = '>tig00000001 len=60 reads=50 class=contig suggestRepeat=no suggestBubble=no ' \
             'suggestCircular=yes trim=10-50'
    seq = 'AGTAGCCAAACTATTTAATGCTAGAGATGCTGCATATCAAAAAATAATCAAACAATTATC'
    new_header, new_seq = trim_seq(header, seq)
    assert new_header == '>tig00000001 len=40 reads=50 class=contig suggestRepeat=no ' \
                         'suggestBubble=no suggestCircular=yes trim=0-40'
    assert new_seq == 'CTATTTAATGCTAGAGATGCTGCATATCAAAAAATAATCA'
