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

import datetime
import os
import textwrap
import sys


def log(message='', end='\n'):
    print(message, file=sys.stderr, flush=True, end=end)


def log_with_wrapping(message):
    terminal_width, _ = get_terminal_size_stderr()
    for line in textwrap.wrap(message, width=terminal_width - 1):
        log(line)


def section_header(text):
    log()
    time = get_timestamp()
    time_str = dim('(' + time + ')')
    header = bold_yellow_underline(text)
    print(header + ' ' + time_str, file=sys.stderr, flush=True)


END_FORMATTING = '\033[0m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'
RED = '\033[31m'
GREEN = '\033[32m'
MAGENTA = '\033[35m'
YELLOW = '\033[93m'
DIM = '\033[2m'


def bold(text):
    return BOLD + text + END_FORMATTING


def bold_yellow(text):
    return YELLOW + BOLD + text + END_FORMATTING


def bold_yellow_underline(text):
    return YELLOW + BOLD + UNDERLINE + text + END_FORMATTING


def dim(text):
    return DIM + text + END_FORMATTING


def red(text):
    return RED + text + END_FORMATTING


def bold_red(text):
    return RED + BOLD + text + END_FORMATTING


def explanation(text, indent_size=4):
    text = ' ' * indent_size + text
    terminal_width, _ = get_terminal_size_stderr()
    for line in textwrap.wrap(text, width=terminal_width - 1):
        log(dim(line))
    log()


def quit_with_error(text):
    log()
    log_with_wrapping(text)
    log()
    sys.exit()


def get_terminal_size_stderr(fallback=(80, 24)):
    """
    Unlike shutil.get_terminal_size, which looks at stdout, this looks at stderr.
    """
    try:
        size = os.get_terminal_size(sys.__stderr__.fileno())
    except (AttributeError, ValueError, OSError):
        size = os.terminal_size(fallback)
    return size


def get_timestamp():
    return '{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now())
