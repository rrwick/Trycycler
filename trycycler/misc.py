#!/usr/bin/env python3
"""
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

import gzip
import itertools
import multiprocessing
import os
import sys

from .log import log


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
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


def get_sequence_file_type(filename):
    """
    Determines whether a file is FASTA or FASTQ.
    """
    if not os.path.isfile(filename):
        sys.exit(f'\nError: could not find {filename}')
    if get_compression_type(filename) == 'gz':
        open_func = gzip.open
    else:  # plain text
        open_func = open
    with open_func(filename, 'rt') as seq_file:
        try:
            first_char = seq_file.read(1)
        except UnicodeDecodeError:
            first_char = ''
    if first_char == '>':
        return 'FASTA'
    elif first_char == '@':
        return 'FASTQ'
    else:
        return 'neither'


def iterate_fastq(filename):
    if get_sequence_file_type(filename) != 'FASTQ':
        sys.exit('\nError: {} is not FASTQ format'.format(filename))
    with get_open_func(filename)(filename, 'rt') as fastq:
        for line in fastq:
            line = line.strip()
            if len(line) == 0:
                continue
            if not line.startswith('@'):
                continue
            name = line[1:].split()[0]
            header = line
            sequence = next(fastq).strip()
            _ = next(fastq)
            qualities = next(fastq).strip()
            yield name, header, sequence, qualities


def get_fastq_stats(filename):
    seq_lengths = [len(s) for _, _, s, _ in iterate_fastq(filename)]
    read_count = len(seq_lengths)
    total_size = sum(seq_lengths)
    n50 = get_n50(seq_lengths)
    return read_count, total_size, n50


def get_n50(seq_lengths):
    seq_lengths = sorted(seq_lengths, reverse=True)
    total_bases = sum(seq_lengths)
    target_bases = total_bases * 0.5
    bases_so_far = 0
    for sequence_length in seq_lengths:
        bases_so_far += sequence_length
        if bases_so_far >= target_bases:
            return sequence_length
    return 0


def count_reads(filename):
    count = 0
    for _ in iterate_fastq(filename):
        count += 1
    return count


def load_fasta(fasta_filename, include_full_header=False):
    if get_compression_type(fasta_filename) == 'gz':
        open_func = gzip.open
    else:  # plain text
        open_func = open
    fasta_seqs = []
    with open_func(fasta_filename, 'rt') as fasta_file:
        name = ''
        sequence = []
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':  # Header line = start of new contig
                if name:
                    if include_full_header:
                        fasta_seqs.append((name.split()[0], name, ''.join(sequence)))
                    else:
                        fasta_seqs.append((name.split()[0], ''.join(sequence)))
                    sequence = []
                name = line[1:]
            else:
                sequence.append(line)
        if name:
            if include_full_header:
                fasta_seqs.append((name.split()[0], name, ''.join(sequence)))
            else:
                fasta_seqs.append((name.split()[0], ''.join(sequence)))
    return fasta_seqs


def get_default_thread_count():
    return min(multiprocessing.cpu_count(), 16)


def write_seq_to_fasta(seq, name, filename):
    with open(filename, 'wt') as f:
        f.write(f'>{name}\n')
        f.write(f'{seq}\n')


REV_COMP_DICT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                 'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
                 'D': 'H', 'H': 'D', 'N': 'N', 'r': 'y', 'y': 'r', 's': 's', 'w': 'w', 'k': 'm',
                 'm': 'k', 'b': 'v', 'v': 'b', 'd': 'h', 'h': 'd', 'n': 'n', '.': '.', '-': '-',
                 '?': '?'}


def complement_base(base):
    try:
        return REV_COMP_DICT[base]
    except KeyError:
        return 'N'


def reverse_complement(seq):
    return ''.join([complement_base(x) for x in seq][::-1])


def remove_duplicates(lst):
    """
    https://stackoverflow.com/questions/480214/how-do-you-remove-duplicates-from-a-list-whilst-preserving-order
    """
    seen = set()
    seen_add = seen.add
    return [x for x in lst if not (x in seen or seen_add(x))]


def means_of_slices(iterable, slice_size):
    """
    https://stackoverflow.com/questions/53373362/the-average-value-of-a-list-in-chunks-of-100-items
    """
    iterator = iter(iterable)
    while True:
        s = list(itertools.islice(iterator, slice_size))
        if s:
            yield sum(s) / len(s)
        else:
            return


def check_python_version():
    if sys.version_info.major < 3 or sys.version_info.minor < 6:
        sys.exit('\nError: Trycycler requires Python 3.6 or later')


def check_output_directory(directory):
    if directory.is_file():
        sys.exit(f'\nError: output directory ({directory}) already exists as a file')
    if directory.is_dir():
        log(f'Output directory ({directory}) already exists - files may be overwritten.')
    else:
        log(f'Creating output directory: {directory}')
        directory.mkdir(parents=True)
    log()


def count_substrings(s, substring):
    """
    https://stackoverflow.com/questions/11476713
    """
    string_size = len(s)
    substring_size = len(substring)
    count = 0
    for i in range(0,string_size-substring_size+1):
        if s[i:i+substring_size] == substring:
            count+=1
    return count


def range_overlap(x1, x2, y1, y2):
    """
    Returns true if the range (x1, x2) overlaps with the range (y1, y2).
    """
    return x1 < y2 and y1 < x2
