"""
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

import os
import pathlib
import subprocess
import tempfile
import itertools
from concurrent.futures import ProcessPoolExecutor, as_completed

from .log import log, dim, red
from .misc import write_seq_to_fasta, reverse_complement
from . import settings

def compute_pairwise_distance(args):
    a, b, pos_sketches, neg_sketches = args
    if a == b:
        return (a, b, 0.0)
    d = min(get_mash_dist(pos_sketches[a], pos_sketches[b]),
            get_mash_dist(pos_sketches[a], neg_sketches[b]))
    return (a, b, d)

def get_mash_dist_matrix(seq_names, seqs, distance_threshold, threads, indent=True):
    max_seq_name_len = max(len(x) for x in seq_names)
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = pathlib.Path(temp_dir)
        pos_sketches, neg_sketches = make_mash_sketches(seq_names, seqs, temp_dir)
        mash_matrix = {}
        
        # Generate (a, b) combinations
        pairs = list(itertools.combinations_with_replacement(seq_names, 2))
        # Prepare input tuples for multiprocessing
        work_items = [(a, b, pos_sketches, neg_sketches) for a, b in pairs]

        # Run parallel distance computation
        with ProcessPoolExecutor(max_workers=threads) as executor:
            futures = [executor.submit(compute_pairwise_distance, args) for args in work_items]
            for future in as_completed(futures):
                a, b, dist = future.result()
                mash_matrix[(a, b)] = dist
                mash_matrix[(b, a)] = dist

        # keep output
        for a in seq_names:
            if indent:
                log('  ', end='')
            log(a, end=':')
            log(' ' * (max_seq_name_len - len(a)), end=' ')
            for b in seq_names:
                dist = mash_matrix[(a, b)]
                if a == b:
                    log(dim(f'{dist:.3f}'), end='')
                elif dist > distance_threshold:
                    log(red(f'{dist:.3f}'), end='')
                else:
                    log(f'{dist:.3f}', end='')
                if b != seq_names[-1]:
                    log('  ', end='')
            log()
        log()

    return mash_matrix


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
    return min(float(parts[2]), settings.MAX_MASH_DISTANCE)
