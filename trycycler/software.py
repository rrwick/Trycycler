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

import subprocess
import sys

from .log import log


def check_minimap2():
    try:
        output = subprocess.check_output(['minimap2', '--version'], stderr=subprocess.STDOUT)
    except FileNotFoundError:
        sys.exit('\nError: unable to find minimap2 - make sure that minimap2 is installed and '
                 'available on the path, then try again.')
    except subprocess.CalledProcessError:
        sys.exit('\nError: unable to determine minimap2 version - make sure that minimap2 is '
                 'correctly installed, then try again.')
    output = output.decode().strip()
    log(f'  minimap2: v{output}')


def check_muscle():
    try:
        output = subprocess.check_output(['muscle', '-version'], stderr=subprocess.STDOUT)
    except FileNotFoundError:
        sys.exit('\nError: unable to find MUSCLE - make sure that MUSCLE is installed and '
                 'available on the path, then try again.')
    except subprocess.CalledProcessError:
        sys.exit('\nError: unable to determine MUSCLE version - make sure that MUSCLE is '
                 'correctly installed, then try again.')
    output = output.decode().strip()
    if 'MUSCLE ' in output:
        output = output.split('MUSCLE ')[1]
    if ' by ' in output:
        output = output.split(' by ')[0]
    output = output.strip()
    log(f'   MUSCLE: {output}')


def check_mash():
    try:
        output = subprocess.check_output(['mash', '--version'], stderr=subprocess.STDOUT)
    except FileNotFoundError:
        sys.exit('\nError: unable to find Mash - make sure that Mash is installed and '
                 'available on the path, then try again.')
    except subprocess.CalledProcessError:
        sys.exit('\nError: unable to determine Mash version - make sure that Mash is '
                 'correctly installed, then try again.')
    output = output.decode().strip()
    log(f'  Mash: v{output}')
