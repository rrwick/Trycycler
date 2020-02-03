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

import random
import string
import sys

from .base_scores import get_per_base_scores
from .circularisation import circularise
from .initial_check import initial_sanity_check
from .log import log, section_header, explanation
from .misc import get_sequence_file_type, load_fasta, get_fastq_stats
from .pairwise import get_pairwise_alignments
from .starting_seq import get_starting_seq, rotate_to_starting_seq
from . import settings


def cluster(args):
    welcome_message()


def welcome_message():
    section_header('Starting Trycycler clustering')
    explanation('Trycycler cluster is a tool for clustering the contigs from multiple different '
                'assemblies (e.g. from different assemblers) into highly-similar groups.')
