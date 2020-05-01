"""
This module contains a class for describing read-to-reference alignments (as made by minimap2) and
related functions.

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
import os
import pathlib
import sys
import subprocess
import tempfile

from .misc import write_seq_to_fasta


class Alignment(object):

    def __init__(self, paf_line):
        line_parts = paf_line.strip().split('\t')
        if len(line_parts) < 11:
            sys.exit('\nError: alignment file does not seem to be in PAF format')

        self.query_name = line_parts[0]
        self.query_length = int(line_parts[1])
        self.query_start = int(line_parts[2])
        self.query_end = int(line_parts[3])
        self.strand = line_parts[4]

        self.ref_name = line_parts[5]
        self.ref_length = int(line_parts[6])
        self.ref_start = int(line_parts[7])
        self.ref_end = int(line_parts[8])

        self.matching_bases = int(line_parts[9])
        self.num_bases = int(line_parts[10])
        self.percent_identity = 100.0 * self.matching_bases / self.num_bases

        self.query_cov = 100.0 * (self.query_end - self.query_start) / self.query_length

        self.cigar, self.alignment_score = None, None
        for part in line_parts:
            if part.startswith('cg:Z:'):
                self.cigar = part[5:]
            if part.startswith('AS:i:'):
                self.alignment_score = int(part[5:])

    def __repr__(self):
        return self.query_name + ':' + str(self.query_start) + '-' + str(self.query_end) + \
               '(' + self.strand + '), ' + \
               self.ref_name + ':' + str(self.ref_start) + '-' + str(self.ref_end) + \
               ' (' + ('%.3f' % self.percent_identity) + '%)'


def align_a_to_b(seq_a, seq_b, preset='asm20'):
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = pathlib.Path(temp_dir)
        temp_a = temp_dir / 'a.fasta'
        temp_b = temp_dir / 'b.fasta'
        write_seq_to_fasta(seq_a, 'A', temp_a)
        write_seq_to_fasta(seq_b, 'B', temp_b)
        with open(os.devnull, 'w') as dev_null:
            out = subprocess.check_output(['minimap2', '--eqx', '-c', '-x', preset,
                                           str(temp_b), str(temp_a)], stderr=dev_null)
    out = out.decode()
    alignment_lines = out.splitlines()
    alignments = [Alignment(x) for x in alignment_lines]
    return alignments


def align_reads_to_seq(reads, seq, threads, include_cigar=True):
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = pathlib.Path(temp_dir)
        temp_fasta = temp_dir / 'seq.fasta'
        write_seq_to_fasta(seq, 'seq', temp_fasta)
        minimap2_command = ['minimap2']
        if include_cigar:
            minimap2_command += ['--eqx', '-c']
        minimap2_command += ['-x', 'map-ont', '-t', str(threads), str(temp_fasta), str(reads)]
        with open(os.devnull, 'w') as dev_null:
            out = subprocess.check_output(minimap2_command, stderr=dev_null)
    out = out.decode()
    alignment_lines = out.splitlines()
    alignments = [Alignment(x) for x in alignment_lines]
    return alignments


def get_best_alignment_per_read(alignments):
    alignments_per_read = collections.defaultdict(list)
    for a in alignments:
        alignments_per_read[a.query_name].append(a)
    best_alignments = []
    for read_name, alignments in alignments_per_read.items():
        alignments = sorted(alignments, key=lambda x: x.alignment_score, reverse=True)
        best_alignments.append(alignments[0])
    return best_alignments
