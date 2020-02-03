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

import pathlib
import string
import subprocess
import sys
import tempfile

from .log import log, section_header, explanation
from .mash import get_mash_dist_matrix
from .misc import get_sequence_file_type, load_fasta, get_fastq_stats


def cluster(args):
    welcome_message()
    check_inputs_and_requirements(args)
    seqs, fasta_names = load_assembly_sequences(args.assemblies)
    phylip = distance_matrix(seqs, args.distance, args.out_dir)
    build_tree(phylip, args.out_dir)


def welcome_message():
    section_header('Starting Trycycler clustering')
    explanation('Trycycler cluster is a tool for clustering the contigs from multiple different '
                'assemblies (e.g. from different assemblers) into highly-similar groups.')


def check_inputs_and_requirements(args):
    check_input_assemblies(args.assemblies)
    check_output_directory(args.out_dir)
    check_required_software()


def check_input_assemblies(filenames):
    if len(filenames) < 2:
        sys.exit('Error: two or more input assemblies are required')
    log(f'Input assemblies:')
    for i, f in enumerate(filenames):
        contig_type = get_sequence_file_type(f)
        if contig_type != 'FASTA':
            sys.exit(f'Error: input assembly file ({f}) is not in FASTA format')
        seqs = load_fasta(f)
        contig_names = set()
        for contig_name, _ in seqs:
            if contig_name in contig_names:
                sys.exit(f'Error: duplicate contig name: {contig_name}')
            contig_names.add(contig_name)
        contig_count = len(seqs)
        total_length = sum(len(s[1]) for s in seqs)
        noun = 'contig' if contig_count == 1 else 'contigs'
        letter = string.ascii_uppercase[i]
        log(f'  {letter}: {f} ({contig_count} {noun}, {total_length:,} bp)')
    log()


def check_output_directory(directory):
    if directory.is_file():
        sys.exit(f'Error: output directory ({directory}) already exists as a file')
    if directory.is_dir():
        log(f'Output directory ({directory}) already exists - files may be overwritten.')
    else:
        log(f'Creating output directory: {directory}')
        directory.mkdir(parents=True)
    log()


def check_required_software():
    pass
    # TODO
    # TODO
    # TODO
    # TODO
    # TODO


def load_assembly_sequences(filenames):
    assembly_seqs, fasta_names = {}, {}
    for i, f in enumerate(filenames):
        letter = string.ascii_uppercase[i]
        fasta_names[letter] = f
        seqs = load_fasta(f)
        for name, seq in seqs:
            full_name = f'{letter}_{name}'
            assembly_seqs[full_name] = seq
    return assembly_seqs, fasta_names


def distance_matrix(seqs, distance, out_dir):
    section_header('Distance matrix')
    explanation('Mash is used to build a distance matrix of all contigs in the assemblies.')
    seq_names = list(seqs.keys())
    mash_matrix = get_mash_dist_matrix(seq_names, seqs, distance)
    return save_matrix_to_phylip(seq_names, mash_matrix, out_dir)


def save_matrix_to_phylip(seq_names, matrix, out_dir):
    phylip = out_dir / 'mash.phylip'
    with open(phylip, 'wt') as f:
        f.write(str(len(seq_names)))
        f.write('\n')
        for a in seq_names:
            f.write(a)
            for b in seq_names:
                f.write('\t')
                f.write(str(matrix[(a, b)]))
            f.write('\n')
    return phylip


def build_tree(phylip, out_dir):
    section_header('Build FastME tree')
    explanation('R (ape and phangorn) are used to build a FastME tree of the relationships '
                'between the contigs.')
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = pathlib.Path(temp_dir)
        tree_script, newick = create_tree_script(temp_dir, phylip)
        log(f'saving tree: {newick}')
        subprocess.check_output(['Rscript', tree_script])


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
