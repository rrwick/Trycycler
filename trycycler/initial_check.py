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

import os
import pathlib
import subprocess
import sys
import tempfile

from .log import log, section_header, explanation, dim, red
from .misc import write_seq_to_fasta, reverse_complement
from . import settings


def initial_sanity_check(seqs):
    section_header('Checking closeness of contigs')
    explanation('Before proceeding, Trycycler ensures that the input contigs appear sufficiently '
                'close to each other to make a consensus. If not, the program will quit and the '
                'user must fix the input contigs (make them more similar to each other) or exclude '
                'some before trying again.')

    seq_names = list(seqs.keys())

    log('Relative sequence lengths:')
    length_matrix = get_length_ratio_matrix(seq_names, seqs)
    check_length_ratios(length_matrix)

    log('Mash distances:')
    mash_matrix = get_mash_dist_matrix(seq_names, seqs)
    check_mash_distances(mash_matrix)

    log('Contigs have passed the initial check - they seem sufficiently close to reconcile.')
    log()


def get_length_ratio_matrix(seq_names, seqs):
    min_threshold, max_threshold = get_length_thresholds()
    length_matrix = {}
    for a in seq_names:
        log(f'  {a}: ', end='')
        for b in seq_names:
            ratio = len(seqs[a]) / len(seqs[b])
            ratio_str = f'{ratio:.3f}'
            if a == b:
                log(dim(ratio_str), end='')
            elif ratio < min_threshold or ratio > max_threshold:
                log(red(ratio_str), end='')
            else:
                log(ratio_str, end='')
            if b != seq_names[-1]:  # if not the last one in the row
                log('  ', end='')
            length_matrix[(a, b)] = ratio
        log()
    log()
    return length_matrix


def get_mash_dist_matrix(seq_names, seqs):
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = pathlib.Path(temp_dir)
        pos_sketches, neg_sketches = make_mash_sketches(seq_names, seqs, temp_dir)
        mash_matrix = {}
        for a in seq_names:
            log(f'  {a}: ', end='')
            for b in seq_names:
                if a == b:
                    distance = 0.0
                else:
                    distance = min(get_mash_dist(pos_sketches[a], pos_sketches[b]),
                                   get_mash_dist(pos_sketches[a], neg_sketches[b]))
                    mash_matrix[(a, b)] = distance
                if a == b:
                    log(dim(f'{distance:.3f}'), end='')
                elif distance > settings.MASH_DISTANCE_THRESHOLD:
                    log(red(f'{distance:.3f}'), end='')
                else:
                    log(f'{distance:.3f}', end='')
                if b != seq_names[-1]:  # if not the last one in the row
                    log('  ', end='')
                mash_matrix[(b, a)] = distance
            log()
    log()
    return mash_matrix


def check_length_ratios(length_matrix):
    """
    This function looks at the whole length ratio matrix and quits the program if any value is too
    out of bounds.
    """
    min_pair = min(length_matrix, key=length_matrix.get)
    max_pair = max(length_matrix, key=length_matrix.get)
    min_ratio = length_matrix[min_pair]
    max_ratio = length_matrix[max_pair]
    min_threshold, max_threshold = get_length_thresholds()
    if min_ratio < min_threshold or max_ratio > max_threshold:
        sys.exit('Error: there is too much length difference between contigs')


def get_length_thresholds():
    min_threshold = settings.LENGTH_DIFFERENCE_THRESHOLD
    max_threshold = 1.0 / min_threshold
    assert min_threshold < 1.0
    assert max_threshold > 1.0
    return min_threshold, max_threshold


def check_mash_distances(mash_matrix):
    """
    This function looks at the whole Mash distance matrix and quits the program if any value is too
    out of bounds.
    """
    max_pair = max(mash_matrix, key=mash_matrix.get)
    max_dist = mash_matrix[max_pair]
    if max_dist > settings.MASH_DISTANCE_THRESHOLD:
        sys.exit(f'Error: there is too much Mash distance between contigs')


def make_mash_sketches(seq_names, seqs, temp_dir):
    pos_sketches, neg_sketches = {}, {}
    for seq_name in seq_names:
        seq_pos = seqs[seq_name]
        seq_neg = reverse_complement(seq_pos)
        fasta_pos = temp_dir / (seq_name + '_pos.fasta')
        fasta_neg = temp_dir / (seq_name + '_neg.fasta')
        write_seq_to_fasta(seq_pos, seq_name, fasta_pos)
        write_seq_to_fasta(seq_neg, seq_name, fasta_neg)
        sketch_pos = temp_dir / (seq_name + '_pos.msh')
        sketch_neg = temp_dir / (seq_name + '_neg.msh')
        with open(os.devnull, 'w') as dev_null:
            _ = subprocess.check_output(['mash', 'sketch', '-n', '-o', str(sketch_pos),
                                         str(fasta_pos)], stderr=dev_null)
        with open(os.devnull, 'w') as dev_null:
            _ = subprocess.check_output(['mash', 'sketch', '-n', '-o', str(sketch_neg),
                                         str(fasta_neg)], stderr=dev_null)
        pos_sketches[seq_name] = sketch_pos
        neg_sketches[seq_name] = sketch_neg
    return pos_sketches, neg_sketches


def get_mash_dist(sketch_a, sketch_b):
    with open(os.devnull, 'w') as dev_null:
        out = subprocess.check_output(['mash', 'dist', str(sketch_a), str(sketch_b)],
                                      stderr=dev_null)
    out = out.decode()
    parts = out.split('\t')
    return float(parts[2])
