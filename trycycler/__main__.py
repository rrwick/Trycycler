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

import argparse
import pathlib
import sys

from .align import align
from .cluster import cluster
from .consensus import consensus
from .help_formatter import MyParser, MyHelpFormatter
from .log import bold
from .misc import get_default_thread_count, check_python_version
from .partition import partition
from .version import __version__


def main():
    check_python_version()
    args = parse_args(sys.argv[1:])
    if args.subparser_name == 'cluster':
        cluster(args)

    elif args.subparser_name == 'align':
        align(args)

    elif args.subparser_name == 'partition':
        partition(args)

    elif args.subparser_name == 'consensus':
        consensus(args)


def parse_args(args):
    parser = MyParser(description=bold('Trycycler: a consensus long-read assembly tool'),
                      formatter_class=MyHelpFormatter, add_help=False)

    subparsers = parser.add_subparsers(title='Commands', dest='subparser_name')
    cluster_subparser(subparsers)
    align_subparser(subparsers)
    partition_subparser(subparsers)
    consensus_subparser(subparsers)

    longest_choice_name = max(len(c) for c in subparsers.choices)
    subparsers.help = 'R|'
    for choice, choice_parser in subparsers.choices.items():
        padding = ' ' * (longest_choice_name - len(choice))
        subparsers.help += choice + ': ' + padding
        d = choice_parser.description
        subparsers.help += d[0].lower() + d[1:]  # don't capitalise the first letter
        subparsers.help += '\n'

    help_args = parser.add_argument_group('Help')
    help_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                           help='Show this help message and exit')
    help_args.add_argument('--version', action='version', version='Trycycler v' + __version__,
                           help="Show program's version number and exit")

    # If no arguments were used, print the base-level help which lists possible commands.
    if len(args) == 0:
        parser.print_help(file=sys.stderr)
        sys.exit(1)

    return parser.parse_args(args)


def cluster_subparser(subparsers):
    group = subparsers.add_parser('cluster', description='cluster assembly contigs by similarity ',
                                  formatter_class=MyHelpFormatter, add_help=False)

    required_args = group.add_argument_group('Required arguments')
    required_args.add_argument('-a', '--assemblies', type=str, required=True, nargs='+',
                               help='Input assemblies whose contigs will be clustered (multiple '
                                    'FASTA files)')
    required_args.add_argument('-o', '--out_dir', type=pathlib.Path, required=True,
                               help='Output directory')

    setting_args = group.add_argument_group('Settings')
    setting_args.add_argument('-d', '--distance', type=float, default=0.01,
                              help='Mash distance complete-linkage clustering threshold')
    setting_args.add_argument('-t', '--threads', type=int, default=get_default_thread_count(),
                              help='Number of threads to use for alignment')

    other_args = group.add_argument_group('Other')
    other_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                            help='Show this help message and exit')
    other_args.add_argument('--version', action='version', version='Trycycler v' + __version__,
                            help="Show program's version number and exit")


def align_subparser(subparsers):
    group = subparsers.add_parser('align', description='Check and align contig sequences',
                                  formatter_class=MyHelpFormatter, add_help=False)

    required_args = group.add_argument_group('Required arguments')
    required_args.add_argument('-c', '--cluster_dir', type=pathlib.Path, required=True,
                               help='Cluster directory (should contain a 1_contigs subdirectory)')
    required_args.add_argument('-r', '--reads', type=str, required=True,
                               help='Long reads (FASTQ format) used to generate the assemblies')

    setting_args = group.add_argument_group('Settings')
    setting_args.add_argument('--linear', action='store_true',
                              help='The input contigs are not circular (default: assume the input '
                                   'contigs are circular)')
    setting_args.add_argument('--max_mash_dist', type=float, default=0.02,
                              help='Maximum allowed pairwise Mash distance')
    setting_args.add_argument('--max_length_diff', type=float, default=1.1,
                              help='Maximum allowed pairwise relative length difference')
    setting_args.add_argument('--min_identity', type=float, default=99.0,
                              help='Minimum allowed pairwise percent identity')
    setting_args.add_argument('-t', '--threads', type=int, default=get_default_thread_count(),
                              help='Number of threads to use for alignment')

    other_args = group.add_argument_group('Other')
    other_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                            help='Show this help message and exit')
    other_args.add_argument('--version', action='version', version='Trycycler v' + __version__,
                            help="Show program's version number and exit")


def partition_subparser(subparsers):
    group = subparsers.add_parser('partition', description='partition reads by cluster',
                                  formatter_class=MyHelpFormatter, add_help=False)

    required_args = group.add_argument_group('Required arguments')
    required_args.add_argument('-c', '--cluster_dirs', type=pathlib.Path, required=True, nargs='+',
                               help='Cluster directories (each should contain 2_all_seqs.fasta '
                                    'and 3_pairwise_alignments files)')
    required_args.add_argument('-r', '--reads', type=str, required=True,
                               help='Long reads (FASTQ format) used to generate the assemblies')

    setting_args = group.add_argument_group('Settings')
    setting_args.add_argument('--min_aligned_len', type=int, default=1000,
                              help='Do not consider reads with less than this many bases aligned')
    setting_args.add_argument('--min_read_cov', type=int, default=90.0,
                              help='Do not consider reads with less than this percentages of '
                                   'their length covered by alignments')
    setting_args.add_argument('-t', '--threads', type=int, default=get_default_thread_count(),
                              help='Number of threads to use for alignment')

    other_args = group.add_argument_group('Other')
    other_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                            help='Show this help message and exit')
    other_args.add_argument('--version', action='version', version='Trycycler v' + __version__,
                            help="Show program's version number and exit")


def consensus_subparser(subparsers):
    group = subparsers.add_parser('consensus', description='Derive a consensus contig sequence',
                                  formatter_class=MyHelpFormatter, add_help=False)

    required_args = group.add_argument_group('Required arguments')
    required_args.add_argument('-c', '--cluster_dir', type=pathlib.Path, required=True,
                               help='Cluster directory (should contain 2_all_seqs.fasta, '
                                    '3_pairwise_alignments and 4_reads.fastq files)')

    setting_args = group.add_argument_group('Settings')
    setting_args.add_argument('--linear', action='store_true',
                              help='The input contigs are not circular (default: assume the input '
                                   'contigs are circular)')
    setting_args.add_argument('--plot_qual', action='store_true',
                              help='Save plots of the per-base quality across the contigs')
    setting_args.add_argument('-t', '--threads', type=int, default=get_default_thread_count(),
                              help='Number of threads to use for alignment')

    other_args = group.add_argument_group('Other')
    other_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                            help='Show this help message and exit')
    other_args.add_argument('--version', action='version', version='Trycycler v' + __version__,
                            help="Show program's version number and exit")


# TODO: parameter checking
#       --max_length_diff should be >1 and <=2
#       --max_mash_dist should be >0 and <=1


if __name__ == '__main__':
    main()
