"""
This module contains some tests for Trycycler. To run them, execute `pytest` from the root
Trycycler directory.

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

import trycycler.software


def test_parse_muscle_version_1():
    output = 'MUSCLE v3.8.1551 by Robert C. Edgar'
    assert trycycler.software.parse_muscle_version(output) == '3.8.1551'


def test_parse_muscle_version_2():
    output = 'muscle 5.0.1430_osx64'
    assert trycycler.software.parse_muscle_version(output) == '5.0.1430'


def test_parse_muscle_version_3():
    output = 'Not the correct output'
    assert trycycler.software.parse_muscle_version(output) == '?'


def test_parse_r_version_1():
    output = 'R version 3.6.2 (2019-12-12) -- "Dark and Stormy Night"\n' \
             'Copyright (C) 2019 The R Foundation for Statistical Computing\n' \
             'Platform: x86_64-apple-darwin19.2.0 (64-bit)'
    assert trycycler.software.parse_r_version(output) == '3.6.2'


def test_parse_r_version_2():
    output = 'Not the correct output'
    assert trycycler.software.parse_r_version(output) == '?'


def test_parse_ape_version_1():
    output = '> packageVersion("ape")\n' \
             '[1] ‘5.3’\n' \
             '>'
    assert trycycler.software.parse_ape_version(output) == '5.3'


def test_parse_ape_version_2():
    output = 'Not the correct output'
    assert trycycler.software.parse_ape_version(output) == '?'


def test_parse_phangorn_version_1():
    output = '> packageVersion("phangorn")\n' \
             '[1] ‘2.5.5’\n' \
             '>'
    assert trycycler.software.parse_ape_version(output) == '2.5.5'


def test_parse_phangorn_version_2():
    output = 'Not the correct output'
    assert trycycler.software.parse_ape_version(output) == '?'
