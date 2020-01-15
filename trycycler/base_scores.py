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


from .log import log, section_header, explanation


def get_per_base_scores(seqs, reads, circular):
    section_header('Per-base quality scores')
    explanation('')

    per_base_scores = {}
    for name, seq in seqs:
        per_base_scores[name] = get_one_seq_per_base_scores(seq, reads, circular)


def get_one_seq_per_base_scores(seq, reads, circular):
    pass
    # TODO: if circular, save a doubled version of the seq to a temp file. If not circular, save
    #       the normal seq instead.

    # TODO: align the reads to the doubled (or not) sequence

    # TODO: clean up redundant alignments (residing entirely in second duplicated half)

    # TODO: GET THE PER-BASE SCORES! (somehow)
