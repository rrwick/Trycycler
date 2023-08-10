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

import math
import os
import pathlib
import random
import subprocess
import sys
import tempfile

from .log import log, section_header, explanation
from .misc import check_input_reads, check_output_directory, count_lines, get_n50, iterate_fastq
from .software import check_minimap2, check_miniasm


def subsample(args):
    random.seed(0)
    welcome_message()
    read_count, read_bases = check_inputs_and_requirements(args)
    if args.genome_size == 'auto':
        genome_size = determine_approx_genome_size(args)
    else:
        genome_size = interpret_genome_size(args.genome_size)
    reads_per_subset = calculate_subsets(read_count, read_bases, genome_size, args.min_read_depth)
    save_subsets(args.reads, args.count, reads_per_subset, args.out_dir)
    finished_message()


def welcome_message():
    section_header('Starting Trycycler read subsampling')
    explanation('Trycycler subsample is a tool for subsampling a long-read set into subsets that '
                'are maximally independent from each other.')


def finished_message():
    section_header('Finished!')
    explanation('You must now assemble each of the subsampled read sets to produce a set of '
                'assemblies you can input into Trycycler cluster.')


def check_inputs_and_requirements(args):
    check_output_directory(args.out_dir)
    read_count, read_bases = check_input_reads(args.reads)
    check_required_software()
    return read_count, read_bases


def check_required_software():
    log('Checking required software:')
    check_minimap2()
    check_miniasm()
    log()


def determine_approx_genome_size(args):
    section_header('Estimating genome size')
    explanation('Since you did not provide a genome size, Trycycler will perform a quick miniasm '
                'assembly of the reads to get a genome size estimate.')

    with tempfile.TemporaryDirectory() as temp_dir, open(os.devnull, 'w') as dev_null:
        temp_dir = pathlib.Path(temp_dir)

        log('Doing an all-vs-all read alignment:')
        paf_filename = temp_dir / 'alignments.paf'
        minimap2_command = ['minimap2', '-x', 'ava-ont', '-t', str(args.threads),
                            str(args.reads), str(args.reads)]
        with open(paf_filename, 'w') as out:
            subprocess.run(minimap2_command, stdout=out, stderr=dev_null)
        alignment_count = count_lines(paf_filename)
        log(f'  {alignment_count:,} alignments')
        log()

        log('Running a miniasm assembly:')
        gfa_filename = temp_dir / 'assembly.gfa'
        miniasm_command = ['miniasm', '-f', str(args.reads), str(paf_filename)]
        with open(gfa_filename, 'w') as out:
            subprocess.run(miniasm_command, stdout=out, stderr=dev_null)
        contig_count, total_size, n50 = get_gfa_stats(gfa_filename, 10000)
        plural = '' if contig_count == 1 else 's'
        log(f'  {contig_count:,} contig{plural} ({total_size:,} bp)')
        log(f'  N50 = {n50:,} bp')
        log()

    if total_size == 0:
        sys.exit('Error: genome size estimate has failed, use --genome_size to manually set a '
                 'genome size')
    log(f'Estimated genome size: {total_size:,} bp')
    log()
    return total_size


def get_gfa_stats(gfa_filename, min_contig_length):
    contig_count, total_size = 0, 0
    seq_lengths = []
    with open(gfa_filename, 'rt') as gfa:
        for line in gfa:
            parts = line.strip().split('\t')
            if parts[0] == 'S':
                contig_count += 1
                seq_length = len(parts[2])
                # Miniasm has a tendency to duplicate small plasmids, so we ignore short contigs.
                if seq_length >= min_contig_length:
                    total_size += seq_length
                seq_lengths.append(seq_length)
    n50 = get_n50(seq_lengths)
    return contig_count, total_size, n50


def interpret_genome_size(genome_size_str: str):
    genome_size = interpret_genome_size_2(genome_size_str)
    if genome_size < 1:
        sys.exit('Error: genome size must be a positive value')
    return genome_size


def interpret_genome_size_2(genome_size_str: str):
    """
    Converts a genome size string (e.g. 5.5M or 3000k) to an integer.
    """
    try:
        return int(round(float(genome_size_str)))
    except ValueError:
        pass
    if genome_size_str.endswith('k') or genome_size_str.endswith('K'):
        try:
            return int(round(1000 * float(genome_size_str[:-1])))
        except ValueError:
            sys.exit('Error: cannot interpret genome size')
    if genome_size_str.endswith('m') or genome_size_str.endswith('M'):
        try:
            return int(round(1000000 * float(genome_size_str[:-1])))
        except ValueError:
            sys.exit('Error: cannot interpret genome size')
    if genome_size_str.endswith('g') or genome_size_str.endswith('G'):
        try:
            return int(round(1000000000 * float(genome_size_str[:-1])))
        except ValueError:
            sys.exit('Error: cannot interpret genome size')
    sys.exit('Error: cannot interpret genome size')


def calculate_subsets(read_count, read_bases, genome_size, min_depth):
    section_header('Calculating subset size')
    explanation('Trycycler will now calculate the number of reads to put in each subset.')

    total_depth = read_bases / genome_size
    mean_read_length = int(round(read_bases / read_count))
    log(f'Total read depth: {total_depth:.1f}x')
    log(f'Mean read length: {mean_read_length:,} bp')
    log()

    if total_depth < min_depth:
        sys.exit('Error: input reads are too shallow to subset')

    log('Calculating subset sizes:')
    log(f'  subset_depth = {min_depth} * log_2(4 * total_depth / {min_depth}) / 2')
    subset_depth = min_depth * math.log(4 * total_depth / min_depth) / (2 * math.log(2))
    log(f'               = {subset_depth:.1f}x')

    subset_ratio = subset_depth / total_depth
    reads_per_subset = int(round(subset_ratio * read_count))
    log(f'  reads per subset: {reads_per_subset:,}')
    log()

    return reads_per_subset


def save_subsets(reads, count, reads_per_subset, out_dir):
    section_header('Subsetting reads')
    explanation('This step shuffles the reads and saves them into subset files.')
    read_order = shuffle_reads(reads)
    read_count = len(read_order)

    for i in range(count):
        log(f'subset {i}:')

        start_1 = int(round(i * read_count / count))
        end_1 = start_1 + reads_per_subset
        if end_1 > read_count:
            start_2 = 0
            end_2 = end_1 - read_count
            end_1 = read_count
            log(f'  reads {start_1+1}-{end_1} and {start_2+1}-{end_2}')
        else:
            start_2, end_2 = None, None
            log(f'  reads {start_1+1}-{end_1}')
        subset = set()
        for j in range(start_1, end_1):
            subset.add(read_order[j])
        if start_2 is not None:
            for j in range(start_2, end_2):
                subset.add(read_order[j])
        assert len(subset) == reads_per_subset

        filename = out_dir / f'sample_{i+1:02d}.fastq'
        log(f'  {filename}')
        total_size = 0
        with open(filename, 'wt') as out_file:
            for j, read in enumerate(iterate_fastq(reads)):
                if j in subset:
                    _, header, seq, qual = read
                    out_file.write(f'{header}\n{seq}\n+\n{qual}\n')
                    total_size += len(seq)
        log(f'  {total_size:,} bp')
        log()


def shuffle_reads(reads):
    log('Shuffling reads... ', end='')
    read_order = []
    for i, _ in enumerate(iterate_fastq(reads)):
        read_order.append(i)
    random.shuffle(read_order)
    log('done\n')
    return read_order
