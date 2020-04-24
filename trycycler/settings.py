"""
This module contains some hard-coded settings used in various parts of Trycycler. These are
probably too low-level to expose to the user (via command-line arguments) but developers may want
to tweak them.

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

MAX_INPUT_CONTIGS = 20

MAX_MASH_DISTANCE = 0.25

START_END_SIZE = 500
START_END_COV_THRESHOLD = 95.0
START_END_IDENTITY_THRESHOLD = 95.0
START_END_GAP_SIZE = 5000
START_END_OVERLAP_SIZE = 200000

CIRCULARISATION_CHOICE_ALIGNMENT_SIZE = 1000
CIRCULARISATION_TRIM_SIZE = 500
CIRCULARISATION_MAX_TRIM_COUNT = 5

RANDOM_STARTING_SEQ_LEN = 250
RANDOM_STARTING_SEQ_TRIAL_COUNT = 1000
RANDOM_STARTING_SEQ_MIN_IDENTITY = 97.5

KNOWN_STARTING_SEQ_MIN_IDENTITY = 95.0
KNOWN_STARTING_SEQ_MIN_COVERAGE = 95.0

CHUNK_COMBINE_SIZE = 100
