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

import collections
import pathlib
import string
import subprocess
import sys
import tempfile

from scipy.cluster.hierarchy import linkage, fcluster
import scipy.spatial.distance as ssd
import numpy as np

from .alignment import align_reads_to_fasta, get_best_alignment_per_read
from .log import log, section_header, explanation, red
from .mash import get_mash_dist_matrix
from .misc import get_sequence_file_type, load_fasta, check_input_reads
from .software import check_mash, check_r, check_ape, check_phangorn


def cluster(args):
    welcome_message()
    assembly_lengths = check_inputs_and_requirements(args)
    seqs, seq_names, fasta_names = load_assembly_sequences(args.assemblies)
    depths, depth_filter = get_contig_depths(args.assemblies, seqs, seq_names, fasta_names,
                                             args.reads, args.threads, assembly_lengths,
                                             args.min_contig_depth)
    seq_names = filter_contigs(args.assemblies, fasta_names, seq_names, seqs, args.min_contig_len,
                               args.min_contig_depth, depth_filter)
    matrix = distance_matrix(seqs, seq_names, args.distance)
    cluster_numbers = complete_linkage(seqs, seq_names, depths, matrix, args.distance, args.out_dir)
    build_tree(seq_names, seqs, depths, matrix, args.out_dir, cluster_numbers)
    finished_message()


def welcome_message():
    section_header('Starting Trycycler clustering')
    explanation('Trycycler cluster is a tool for clustering the contigs from multiple different '
                'assemblies (e.g. from different assemblers) into highly-similar groups.')


def finished_message():
    section_header('Finished!')
    explanation('Now you must decide which clusters are good (i.e. contain well-assembled contigs '
                'for replicons in the genome) and which are bad (i.e. contain incomplete or '
                'spurious contigs). You can then delete the directories corresponding to the bad '
                'clusters and proceed to the next step in the pipeline: trycycler reconcile.')


def check_inputs_and_requirements(args):
    assembly_lengths = check_input_assemblies(args.assemblies)
    check_output_directory(args.out_dir)
    check_input_reads(args.reads)
    check_required_software()
    return assembly_lengths


def check_input_assemblies(filenames):
    if len(filenames) < 2:
        sys.exit('Error: two or more input assemblies are required')
    if len(filenames) > 26:
        sys.exit('Error: you cannot give more than 26 input assemblies')
    log_lines = []
    total_lengths = {}
    for i, f in enumerate(filenames):
        contig_type = get_sequence_file_type(f)
        if contig_type != 'FASTA':
            sys.exit(f'\nError: input assembly file ({f}) is not in FASTA format')
        seqs = load_fasta(f)
        contig_names = set()
        for contig_name, _ in seqs:
            if contig_name in contig_names:
                sys.exit(f'\nError: duplicate contig name: {contig_name}')
            contig_names.add(contig_name)
        contig_count = len(seqs)
        total_length = sum(len(s[1]) for s in seqs)
        letter = string.ascii_uppercase[i]
        total_lengths[letter] = total_length
        log_lines.append((letter, f, contig_count, total_length))

    longest_filename = max(len(f) for _, f, _, _ in log_lines)
    longest_length = max(len(f'{l:,}') for _, _, _, l in log_lines)
    longest_contigs = max(len(f'{c:,}') for _, _, c, _ in log_lines)
    log(f'Input assemblies:')
    for letter, f, contig_count, total_length in log_lines:
        noun = 'contig' if contig_count == 1 else 'contigs'
        padded_f = f.ljust(longest_filename)
        padded_length = f'{total_length:,}'.rjust(longest_length)
        padded_contigs = f'{contig_count:,}'.rjust(longest_contigs)
        log(f'  {letter}: {padded_f} ({padded_length} bp, {padded_contigs} {noun})')
    log()
    return total_lengths


def check_output_directory(directory):
    if directory.is_file():
        sys.exit(f'Error: output directory ({directory}) already exists as a file')
    if directory.is_dir():
        if len(list(directory.iterdir())) > 0:
            sys.exit(f'Error: output directory ({directory}) already exists and is not empty')
        else:
            log(f'Output directory already exists: {directory}')
    else:
        log(f'Creating output directory: {directory}')
        directory.mkdir(parents=True)
    log()


def check_required_software():
    log('Checking required software:')
    check_mash()
    check_r()
    check_ape()
    check_phangorn()
    log()


def load_assembly_sequences(filenames):
    assembly_seqs, fasta_names = {}, {}
    for i, f in enumerate(filenames):
        letter = string.ascii_uppercase[i]
        fasta_names[letter] = f
        seqs = load_fasta(f)
        for name, seq in seqs:
            if len(seq) < 32:  # very short sequences won't work in Mash
                continue
            name = name.replace('#', '_')  # hashes in names can cause downstream problems
            full_name = f'{letter}_{name}'
            assembly_seqs[full_name] = seq
    seq_names = list(assembly_seqs.keys())
    return assembly_seqs, seq_names, fasta_names


def get_contig_depths(assembly_filenames, seqs, seq_names, fasta_names, read_filename, threads,
                      assembly_lengths, min_contig_depth):
    section_header('Getting contig depths')
    explanation('Trycycler now aligns the reads to each of the assemblies to assign a read depth '
                'value to each of the contigs. Contigs displayed in red have a low read depth and '
                'will be filtered out.')

    name_to_letter = {v: k for k, v in fasta_names.items()}
    depths = {n: 0.0 for n in seq_names}
    depth_filter = {}
    for assembly_filename in assembly_filenames:
        letter = name_to_letter[assembly_filename]
        log(f'{letter} ({assembly_filename}):')
        alignments = align_reads_to_fasta(read_filename, assembly_filename, threads)
        alignments = get_best_alignment_per_read(alignments)
        log(f'  {len(alignments):,} alignments', end='')
        total_alignment_length = sum((a.ref_end - a.ref_start) for a in alignments)
        mean_depth = total_alignment_length / assembly_lengths[letter]
        threshold_depth = mean_depth * min_contig_depth
        log(f', mean depth = {mean_depth:.1f}x')
        for a in alignments:
            seq_name = f'{letter}_{a.ref_name}'.replace('#', '_')
            depth_contribution = (a.ref_end - a.ref_start) / a.ref_length
            assert depth_contribution > 0.0
            depths[seq_name] += depth_contribution
        log_lines = []
        for seq_name, depth in depths.items():
            if seq_name.startswith(letter + '_'):
                seq_length = len(seqs[seq_name])
                log_lines.append((seq_name, seq_length, depth))
                depth_filter[seq_name] = (depth >= threshold_depth)

        longest_seq_name = max(len(n) + 1 for n, _, _ in log_lines)
        longest_length = max(len(f'{l:,}') for _, l, _, in log_lines)
        longest_depth = max(len(f'{d:.1f}') for _, _, d in log_lines)
        for seq_name, seq_length, depth in log_lines:
            padded_seq_name = f'{seq_name}:'.ljust(longest_seq_name)
            padded_length = f'{seq_length:,}'.rjust(longest_length)
            padded_depth = f'{depth:.1f}'.rjust(longest_depth)
            output_line = f'  {padded_seq_name} {padded_length} bp, {padded_depth}x'
            if depth_filter[seq_name]:
                log(output_line)
            else:
                log(red(output_line))
        log()

    return depths, depth_filter


def filter_contigs(assembly_filenames, fasta_names, seq_names, seqs, min_contig_len,
                   min_contig_depth, depth_filter):
    section_header('Filtering contigs')
    explanation(f'Contigs are now filtered out if they are too short (<{min_contig_len:,} bp) or '
                f'too low depth (<{min_contig_depth:.1f} times the mean depth).')

    final_seq_names = []
    name_to_letter = {v: k for k, v in fasta_names.items()}
    for assembly_filename in assembly_filenames:
        pass_names, fail_length_names, fail_depth_names = [], [], []
        letter = name_to_letter[assembly_filename]
        log(f'{letter} ({assembly_filename}):')
        for seq_name in seq_names:
            if seq_name.startswith(letter + '_'):
                if len(seqs[seq_name]) < min_contig_len:
                    fail_length_names.append(seq_name)
                elif not depth_filter[seq_name]:
                    fail_depth_names.append(seq_name)
                else:
                    pass_names.append(seq_name)
        if not fail_depth_names and not fail_length_names:
            log('  all contigs passed filtering')
        else:
            if pass_names:
                log('  passed filtering:         ' + ', '.join(pass_names))
            if fail_length_names:
                log('  removed for short length: ' + ', '.join(fail_length_names))
            if fail_depth_names:
                log('  removed for low depth:    ' + ', '.join(fail_depth_names))
        final_seq_names += pass_names
        log()

    return final_seq_names


def distance_matrix(seqs, seq_names, distance):
    section_header('Building distance matrix')
    explanation('Mash is used to build a distance matrix of all contigs in the assemblies.')
    mash_matrix = get_mash_dist_matrix(seq_names, seqs, distance, indent=False)
    return mash_matrix


def save_matrix_to_phylip(seq_names, seqs, depths, matrix, out_dir, cluster_numbers):
    phylip = out_dir / 'contigs.phylip'
    log(f'saving distance matrix: {phylip}')
    with open(phylip, 'wt') as f:
        f.write(str(len(seq_names)))
        f.write('\n')
        for a in seq_names:
            seq_len = len(seqs[a])
            depth = depths[a]
            f.write(f'cluster_{cluster_numbers[a]}_{a}_{seq_len}_bp_{depth:.1f}x')
            for b in seq_names:
                f.write('\t')
                f.write(str(matrix[(a, b)]))
            f.write('\n')
    return phylip


def build_tree(seq_names, seqs, depths, matrix, out_dir, cluster_numbers):
    section_header('Building FastME tree')
    explanation('R (ape and phangorn) are used to build a FastME tree of the relationships '
                'between the contigs.')
    phylip = save_matrix_to_phylip(seq_names, seqs, depths, matrix, out_dir, cluster_numbers)
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = pathlib.Path(temp_dir)
        tree_script, newick = create_tree_script(temp_dir, phylip)
        log(f'saving tree: {newick}')
        subprocess.check_output(['Rscript', tree_script])
    log()


def create_tree_script(temp_dir, phylip):
    phylip = str(phylip.absolute())
    newick = phylip.replace('.phylip', '.newick')
    tree_script = temp_dir / 'tree.R'
    with open(tree_script, 'wt') as f:
        f.write('#!/usr/bin/Rscript\n')
        f.write('library(ape)\n')
        f.write('library(phangorn)\n')
        f.write(f'distances <- readDist("{phylip}")\n')
        f.write('tree <- fastme.bal(distances)\n')
        f.write('tree$edge.length <- pmax(tree$edge.length, 0.0)\n')
        f.write('tree <- midpoint(tree)\n')
        f.write(f'write.tree(tree, "{newick}")\n')
    return str(tree_script), pathlib.Path(newick).relative_to(pathlib.Path.cwd())


def complete_linkage(seqs, seq_names, depths, distances, threshold, out_dir):
    section_header('Clustering')
    explanation('The contigs are now split into clusters using a complete-linkage hierarchical '
                'approach.')

    # Build the distance matrix as a numpy array.
    matrix = []
    for x in seq_names:
        row = [distances[(x, y)] for y in seq_names]
        matrix.append(row)
    matrix = np.array(matrix)
    np.set_printoptions(precision=3, edgeitems=10, linewidth=1000)

    dist_array = ssd.squareform(matrix)
    z = linkage(dist_array, 'complete')
    result = fcluster(z, threshold, criterion='distance')

    clusters = collections.defaultdict(list)
    for name, cluster_name in zip(seq_names, result):
        clusters[cluster_name].append((name, seqs[name]))

    # Sort clusters based on the sum of their sequence lengths. This makes cluster with longer
    # sequences (e.g. chromosomes) and more representative sequences come earlier in the list.
    cluster_names = sorted(clusters.keys(), reverse=True,
                           key=lambda c: sum(len(s[1]) for s in clusters[c]))

    cluster_dir = out_dir / f'cluster_{0:03d}' / '1_contigs'
    longest_fasta = max(len(str(cluster_dir / f'{name}.fasta')) + 1 for name in seq_names)
    longest_length = max(len(f'{len(s):,}') for s in seqs.values())
    longest_depth = max(len(f'{d:.1f}') for d in depths.values())

    cluster_numbers = {}
    for i, cluster_name in enumerate(cluster_names):
        cluster_num = i + 1
        cluster_dir = out_dir / f'cluster_{cluster_num:03d}' / '1_contigs'
        cluster_dir.mkdir(parents=True)
        log(f'{cluster_dir}:')
        log_lines = []
        for name, seq in clusters[cluster_name]:
            cluster_numbers[name] = cluster_num
            seq_fasta = cluster_dir / f'{name}.fasta'
            seq_length = len(seq)
            seq_depth = depths[name]
            with open(seq_fasta, 'wt') as f:
                f.write(f'>{name}\n')
                f.write(f'{seq}\n')
            log_lines.append((seq_fasta, seq_length, seq_depth))

        for seq_fasta, seq_length, seq_depth in log_lines:
            padded_fasta = f'{seq_fasta}:'.ljust(longest_fasta)
            padded_length = f'{seq_length:,}'.rjust(longest_length)
            padded_depth = f'{seq_depth:.1f}'.rjust(longest_depth)
            log(f'  {padded_fasta} {padded_length} bp, {padded_depth}x')
        log()

    return cluster_numbers
