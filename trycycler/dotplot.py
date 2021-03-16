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
from PIL import Image, ImageDraw, ImageFont

from .log import log, section_header, explanation
from .misc import reverse_complement
from .reconcile import check_cluster_directory, check_input_contigs, load_contig_sequences


# Some hard-coded settings for the dotplot generation.
INITIAL_TOP_LEFT_GAP = 200
BORDER_GAP = 30
BETWEEN_SEQ_GAP = 30
OUTLINE_WIDTH = 4
TEXT_GAP = 5
BACKGROUND_COLOUR = (235, 235, 235)
SELF_VS_SELF_COLOUR = (220, 220, 220)
SELF_VS_OTHER_COLOUR = 'white'
TEXT_COLOUR = 'black'
FORWARD_STRAND_DOT_COLOUR = (0, 0, 127)
REVERSE_STRAND_DOT_COLOUR = (127, 0, 0)
MAX_FONT_SIZE = 50


def dotplot(args):
    welcome_message()
    check_inputs_and_requirements(args)
    seqs, fasta_names = load_contig_sequences(args.cluster_dir)
    seq_names = sorted(seqs.keys())
    image = create_dotplots(seq_names, seqs, args)

    dotplots_filename = args.cluster_dir / 'dotplots.png'
    image.save(dotplots_filename)
    log(f'Saving dotplots to: {dotplots_filename}')
    log()
    finished_message()


def welcome_message():
    section_header('Starting Trycycler dotplot')
    explanation('Trycycler dotplot is a tool for drawing all pairwise dotplots for a set of '
                'contigs in a cluster. This step is optional but may help with interpretation '
                'of difficult clusters when performing reconciliation.')


def finished_message():
    section_header('Finished!')
    explanation('You can now examine the dotplots.png file in '
                'Now you must decide which clusters are good (i.e. contain well-assembled contigs '
                'for replicons in the genome) and which are bad (i.e. contain incomplete or '
                'spurious contigs). You can then delete the directories corresponding to the bad '
                'clusters and proceed to the next step in the pipeline: trycycler reconcile.')


def check_inputs_and_requirements(args):
    check_cluster_directory(args.cluster_dir)
    check_input_contigs(args.cluster_dir)


def create_dotplots(seq_names, seqs, args):
    section_header('Creating dotplots')
    explanation('Each pairwise combination of sequences is now compared to generate a dotplot, '
                'all of which will be combined into a single image.')
    seq_names = sorted(seq_names)

    # We create an initial image to test the label sizes.
    start_positions, end_positions, bp_per_pixel = \
        get_positions(args, seq_names, seqs, args.kmer, INITIAL_TOP_LEFT_GAP, BORDER_GAP)
    image = Image.new('RGB', (args.dot_plot_res, args.dot_plot_res), BACKGROUND_COLOUR)
    min_font_size, max_text_height = \
        draw_labels(image, seq_names, start_positions, end_positions, font_size=None)

    # Now that we know the values for min_font_size and max_text_height, we start over, this time
    # limiting the font size to the minimum (so all text is the same size) and readjusting the
    # top-left gap (so it isn't bigger than necessary).
    new_top_left_gap = max_text_height + BORDER_GAP
    start_positions, end_positions, bp_per_pixel = \
        get_positions(args, seq_names, seqs, args.kmer, new_top_left_gap, BORDER_GAP)
    image = Image.new('RGB', (args.dot_plot_res, args.dot_plot_res), BACKGROUND_COLOUR)
    draw_sequence_boxes(image, seq_names, start_positions, end_positions)
    draw_labels(image, seq_names, start_positions, end_positions, font_size=min_font_size)

    for name_a in seq_names:
        seq_a = seqs[name_a]
        for name_b in seq_names:
            seq_b = seqs[name_b]
            log(f'  {name_a} vs {name_b}')
            draw_dots(image, name_a, name_b, seq_a, seq_b, start_positions, bp_per_pixel, args.kmer)

    # The boxes are drawn once more, this time with no fill. This is to overwrite any dots which
    # leaked into the outline, which would look messy.
    draw_sequence_boxes(image, seq_names, start_positions, end_positions, fill=False)

    log()
    return image


def get_positions(args, seq_names, seqs, kmer_size, top_left_gap, bottom_right_gap):
    """
    This function returns the image coordinates that start/end each sequence. Since the dot plot is
    symmetrical, there is only one start/end per sequence (used for both x and y coordinates).
    """
    all_gaps = top_left_gap + bottom_right_gap + BETWEEN_SEQ_GAP * (len(seq_names) - 1)
    pixels_for_sequence = args.dot_plot_res - all_gaps
    total_seq_length = sum(len(seqs[n]) - kmer_size for n in seq_names)
    bp_per_pixel = total_seq_length / pixels_for_sequence

    start_positions, end_positions = {}, {}
    current_pos = top_left_gap
    for name in seq_names:
        start_positions[name] = current_pos
        rect_size = int(round(len(seqs[name]) / bp_per_pixel))
        current_pos += rect_size
        end_positions[name] = current_pos
        current_pos += BETWEEN_SEQ_GAP

    return start_positions, end_positions, bp_per_pixel


def draw_sequence_boxes(image, seq_names, start_positions, end_positions, fill=True):
    """
    This function draws the box for each of the dot plots in the full image.
    """
    draw = ImageDraw.Draw(image)
    for name_a in seq_names:
        start_a = start_positions[name_a] - OUTLINE_WIDTH
        end_a = end_positions[name_a] + OUTLINE_WIDTH
        for name_b in seq_names:
            start_b = start_positions[name_b] - OUTLINE_WIDTH
            end_b = end_positions[name_b] + OUTLINE_WIDTH
            if fill:
                if name_a == name_b:
                    fill_colour = SELF_VS_SELF_COLOUR
                else:
                    fill_colour = SELF_VS_OTHER_COLOUR
                draw.rectangle([(start_a, start_b), (end_a, end_b)],
                               fill=fill_colour, outline='black', width=OUTLINE_WIDTH)
            else:
                draw.rectangle([(start_a, start_b), (end_a, end_b)],
                               outline='black', width=OUTLINE_WIDTH)


def draw_labels(image, seq_names, start_positions, end_positions, font_size):
    draw = ImageDraw.Draw(image)
    if font_size is None:
        font_size = MAX_FONT_SIZE

    min_pos = min(p for p in start_positions.values())
    font_sizes, text_heights = [], []
    for name in seq_names:
        font, text_width, text_height, font_size = \
            get_font(draw, name, font_size, start_positions[name], end_positions[name])
        font_sizes.append(font_size)
        text_heights.append(text_height)

        # Horizontal labels on the top side.
        pos = min_pos - text_height - OUTLINE_WIDTH - TEXT_GAP
        draw.text((start_positions[name], pos), name, font=font, fill=TEXT_COLOUR)

        # Vertical labels on the left side.
        # text = Image.new('L', (text_width, text_height))
        # draw2 = ImageDraw.Draw(text)
        # draw2.text((0, 0), label, font=font, fill='black')
        # rotated_text = text.rotate(90, expand=1)
        # image.paste(rotated_text, (pos, end_positions[name]), rotated_text)
        image_2 = Image.new('RGBA', (text_width, text_height), BACKGROUND_COLOUR)
        draw_2 = ImageDraw.Draw(image_2)
        draw_2.text((0, 0), text=name, font=font, fill=TEXT_COLOUR)
        image_2 = image_2.rotate(90, expand=1)
        sx, sy = image_2.size
        image.paste(image_2, (pos, end_positions[name] - sy, pos + sx, end_positions[name]),
                    image_2)
    return min(font_sizes), max(text_heights)


def draw_dots(image, name_a, name_b, seq_a, seq_b, start_positions, bp_per_pixel, kmer_size):
    pixels = image.load()
    # draw = ImageDraw.Draw(image)
    a_start_pos = start_positions[name_a]
    b_start_pos = start_positions[name_b]

    a_forward_kmers, a_reverse_kmers = get_all_kmer_positions(kmer_size, seq_a)

    for j in range(len(seq_b) - kmer_size + 1):
        j_pixel = int(round(j / bp_per_pixel)) + b_start_pos
        k = seq_b[j:j+kmer_size]
        if k in a_reverse_kmers:
            for i in a_reverse_kmers[k]:
                i_pixel = int(round(i / bp_per_pixel)) + a_start_pos
                pixels[i_pixel, j_pixel] = REVERSE_STRAND_DOT_COLOUR
                # pixels[i_pixel+1, j_pixel] = REVERSE_STRAND_DOT_COLOUR
                # pixels[i_pixel-1, j_pixel] = REVERSE_STRAND_DOT_COLOUR
                # pixels[i_pixel, j_pixel+1] = REVERSE_STRAND_DOT_COLOUR
                # pixels[i_pixel, j_pixel-1] = REVERSE_STRAND_DOT_COLOUR
                # draw.point((i_pixel, j_pixel), fill=REVERSE_STRAND_DOT_COLOUR)
        if k in a_forward_kmers:
            for i in a_forward_kmers[k]:
                i_pixel = int(round(i / bp_per_pixel)) + a_start_pos
                pixels[i_pixel, j_pixel] = FORWARD_STRAND_DOT_COLOUR
                # pixels[i_pixel+1, j_pixel] = FORWARD_STRAND_DOT_COLOUR
                # pixels[i_pixel-1, j_pixel] = FORWARD_STRAND_DOT_COLOUR
                # pixels[i_pixel, j_pixel+1] = FORWARD_STRAND_DOT_COLOUR
                # pixels[i_pixel, j_pixel-1] = FORWARD_STRAND_DOT_COLOUR
                # draw.point((i_pixel, j_pixel), fill=FORWARD_STRAND_DOT_COLOUR)


def get_all_kmer_positions(kmer_size, seq):
    forward_kmers, reverse_kmers = collections.defaultdict(list), collections.defaultdict(list)
    rev_comp_seq = reverse_complement(seq)
    seq_len = len(seq) - kmer_size + 1
    for i in range(seq_len):
        k = seq[i:i+kmer_size]
        forward_kmers[k].append(i)
        k = rev_comp_seq[i:i+kmer_size]
        reverse_kmers[k].append(seq_len - i)
    assert len(forward_kmers) < len(seq)
    assert len(reverse_kmers) < len(seq)
    return forward_kmers, reverse_kmers


def get_font(draw, label, font_size, start_position, end_position):
    font, is_default_font = load_font(font_size)

    # If we had to resort to the default font, then we can't size it.
    if is_default_font:
        text_width, text_height = draw.textsize(label, font=font)
        return font, text_width, text_height

    # If we have a size-able font, then we adjust the size down until the text fits in the
    # available space.
    available_width = end_position - start_position
    text_width, text_height = draw.textsize(label, font=font)
    while text_width > available_width:
        font_size -= 1
        font, _ = load_font(font_size)
        text_width, text_height = draw.textsize(label, font=font)
    return font, text_width, text_height, font_size


def load_font(font_size):
    try:
        return ImageFont.truetype('DejaVuSans.ttf', font_size), False
    except OSError:
        pass
    try:
        return ImageFont.truetype('OpenSans-Regular.ttf', font_size), False
    except OSError:
        pass
    try:
        return ImageFont.truetype('Arial.ttf', font_size), False
    except OSError:
        pass
    return ImageFont.load_default(), True
