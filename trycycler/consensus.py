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
import math
import pathlib
import sys
import tempfile

from .alignment import align_reads_to_seq, get_best_alignment_per_read
from .log import log, section_header, explanation
from .misc import get_sequence_file_type, load_fasta, check_input_reads, range_overlap, \
    load_fastq_as_dict
from .software import check_minimap2
from . import settings


def consensus(args):
    welcome_message()
    check_inputs_and_requirements(args)
    circular = not args.linear

    seqs, seq_names, seq_lengths = load_seqs(args.cluster_dir)
    msa_seqs, msa_names, msa_length = load_msa(args.cluster_dir)
    sanity_check_msa(seqs, seq_names, seq_lengths, msa_seqs, msa_names, msa_length)

    chunks = partition_msa(msa_seqs, msa_names, msa_length, settings.CHUNK_COMBINE_SIZE)
    save_chunks_to_gfa(chunks, args.cluster_dir / '5_chunked_sequence.gfa', len(msa_names),
                       circular)

    consensus_seq_with_gaps, consensus_seq_without_gaps = make_initial_consensus(chunks)
    save_seqs_to_fasta({args.cluster_dir.name + '_consensus': consensus_seq_without_gaps},
                       args.cluster_dir / '6_initial_consensus.fasta')

    index_reads(args.cluster_dir, chunks, consensus_seq_with_gaps, consensus_seq_without_gaps,
                circular, args.threads, args.min_read_cov, args.min_aligned_len)
    choose_best_chunk_options(chunks, args.cluster_dir, args.threads, args.verbose, circular)

    final_consensus_seq = ''.join([c.best_seq for c in chunks]).replace('-', '')
    save_seqs_to_fasta({args.cluster_dir.name + '_consensus': final_consensus_seq},
                       args.cluster_dir / '7_final_consensus.fasta')


def welcome_message():
    section_header('Starting Trycycler consensus')
    explanation('Trycycler consensus is the final stage of the Trycycler pipeline. It operates '
                'on one replicon (i.e. cluster) at a time. It takes the multiple sequence '
                'alignment of alternative contig sequences and combines them into a single '
                'consensus sequence. Where needed, it will use read alignments to help choose '
                'which variants to include/exclude from the consensus sequence. If all goes '
                'well, the final consensus will be free of any large-scale errors.')


def check_inputs_and_requirements(args):
    check_input_reads(args.cluster_dir / '4_reads.fastq')
    check_cluster_directory(args.cluster_dir)
    check_seqs(args.cluster_dir)
    check_required_software()


def partition_msa(msa_seqs, seq_names, msa_length, combine_size):
    section_header('Partitioning MSA')
    explanation('The multiple sequence alignment is now partitioned into chunks. Chunks where the '
                'input contig sequences are all in agreement are called "same" chunks, and those '
                'where the input contig sequences disagree are called "different" chunks. '
                'The consensus sequence will be made by choosing a best option for each of the '
                'different chunks.')

    chunks, chunk_count = [], 0
    current_chunk = Chunk()
    for msa_seq in msa_seqs.values():
        assert len(msa_seq) == msa_length

    for i in range(msa_length):
        bases = {n: msa_seqs[n][i] for n in seq_names}
        if not current_chunk.can_add_bases(bases):
            chunks.append(current_chunk)
            chunk_count += 1
            if chunk_count % 100 == 0:
                log(f'\rchunks: {chunk_count:,}', end='')
            current_chunk = Chunk()
        current_chunk.add_bases(bases)

    if not current_chunk.is_empty():
        chunks.append(current_chunk)
        chunk_count += 1

    same_count = len([c for c in chunks if c.type == 'same'])
    different_count = len([c for c in chunks if c.type == 'different'])
    log(f'\rchunks: {chunk_count:,} ({same_count:,} same, {different_count:,} different)', end='')
    log()

    sanity_check_chunks(chunks, msa_length)
    chunks = combine_chunks(chunks, combine_size)
    sanity_check_chunks(chunks, msa_length)
    log()

    return chunks


def make_initial_consensus(chunks):
    section_header('Initial consensus')
    explanation('Trycycler now makes an initial consensus sequence by choosing a sequence for '
                'each of the different chunks. The chosen sequence is the one with the lowest '
                'total Hamming distance to the other sequences. For example, a chunk with options '
                'of TT, TT, CC, CC and TA will give a consensus of TT. If the total Hamming '
                'distances fail to break a tie or if all sequences differ, the chunk will be '
                'flagged for read-based assessment.')
    total_length, different_needing_assessment, different_not_needing_assessment = 0, 0, 0
    for i, chunk in enumerate(chunks):
        chunk.prepare_chunk()
        if chunk.type == 'different':
            if chunk.needs_assessment:
                different_needing_assessment += 1
            else:
                different_not_needing_assessment += 1

        assert chunk.best_seq is not None
        total_length += len(chunk.best_seq.replace('-', ''))

    log(f'Consensus length: {total_length:,} bp')
    log('')
    log(f'Different chunks needing assessment:     {different_needing_assessment:,}')
    log(f'Different chunks not needing assessment: {different_not_needing_assessment:,}')
    log('')
    consensus_seq_with_gaps = ''.join([c.best_seq for c in chunks])
    consensus_seq_without_gaps = consensus_seq_with_gaps.replace('-', '')
    return consensus_seq_with_gaps, consensus_seq_without_gaps


def index_reads(cluster_dir, chunks, consensus_seq_with_gaps, consensus_seq_without_gaps,
                circular, threads, min_read_cov, min_aligned_len):
    section_header('Indexing reads')
    explanation('Trycycler now aligns all reads to the initial consensus to form an index of '
                'which reads span each of the chunks. This makes the following step faster, as '
                'only relevant reads will be used when conducting read-based assessment of chunks.')

    needs_assessment_count = len([c for c in chunks if c.needs_assessment])
    if needs_assessment_count == 0:
        log('No chunks need read-based assessment. Skipping this step.\n')
        return

    ungapped_to_gapped = make_ungapped_pos_to_gapped_pos_dict(consensus_seq_with_gaps,
                                                              consensus_seq_without_gaps)
    reads = cluster_dir / '4_reads.fastq'
    ungapped_len = len(consensus_seq_without_gaps)
    gapped_len = len(consensus_seq_with_gaps)

    log('Aligning reads to the initial consensus:')
    if circular:
        ref_seq = consensus_seq_without_gaps + consensus_seq_without_gaps
        alignments = align_reads_to_seq(reads, ref_seq, threads)
        alignments = [a for a in alignments if a.ref_start < ungapped_len]
    else:
        ref_seq = consensus_seq_without_gaps
        alignments = align_reads_to_seq(reads, ref_seq, threads)
    log(f'  {len(alignments):,} alignments')
    log()

    log('Filtering for best alignment per read:')
    alignments = get_best_alignment_per_read(alignments)
    alignments = [a for a in alignments if a.query_cov >= min_read_cov
                  and a.query_end - a.query_start >= min_aligned_len]
    log(f'  {len(alignments):,} alignments')
    log()

    chunk_start = 0
    completed = 0
    for chunk in chunks:
        chunk_end = chunk_start + chunk.get_length()
        if chunk.needs_assessment:
            assert chunk.type == 'different'
            for a in alignments:
                gapped_start = ungapped_to_gapped[a.ref_start]
                if a.ref_end <= ungapped_len:  # if we're not spanning the circular gap
                    gapped_end = ungapped_to_gapped[a.ref_end]
                    if range_overlap(chunk_start, chunk_end, gapped_start, gapped_end):
                        chunk.read_names.add(a.query_name)
                else:  # if we are spanning the circular gap
                    gapped_end = ungapped_to_gapped[a.ref_end - ungapped_len]
                    if range_overlap(chunk_start, chunk_end, gapped_start, gapped_len) or \
                            range_overlap(chunk_start, chunk_end, 0, gapped_end):
                        chunk.read_names.add(a.query_name)
            completed += 1
            log(f'\rGathering reads for chunks: {completed:,} / {needs_assessment_count:,}', end='')
        chunk_start = chunk_end
    assert chunk_start == len(consensus_seq_with_gaps)
    log('\n')


def make_ungapped_pos_to_gapped_pos_dict(consensus_seq_with_gaps, consensus_seq_without_gaps):
    ungapped_to_gapped = {}
    gapped_pos, ungapped_pos = 0, 0
    for base in consensus_seq_with_gaps:
        ungapped_to_gapped[ungapped_pos] = gapped_pos
        gapped_pos += 1
        if base != '-':
            ungapped_pos += 1
    ungapped_to_gapped[ungapped_pos] = gapped_pos
    assert gapped_pos == len(consensus_seq_with_gaps)
    assert ungapped_pos == len(consensus_seq_without_gaps)
    return ungapped_to_gapped


def choose_best_chunk_options(chunks, cluster_dir, threads, verbose, circular):
    section_header('Choosing best options with reads')
    explanation('For each of the chunks to be assessed, Trycycler now aligns the relevant reads '
                'to each alternative sequence. Whichever option gives the strongest read '
                'alignments (defined as the total alignment score for each of the read\'s best '
                'alignment) is chosen as the best. This should result in a consensus sequence '
                'which is more accurate than the initial consensus.')
    needs_assessment_count = len([c for c in chunks if c.needs_assessment])
    if needs_assessment_count == 0:
        log('No chunks need read-based assessment. Skipping this step.\n')
        return

    reads = load_fastq_as_dict(cluster_dir)

    new_best_seqs = {}
    completed, kept, changed = 0, 0, 0
    for i, chunk in enumerate(chunks):
        if not chunk.needs_assessment:
            continue
        assert chunk.type == 'different'

        best_seq, output_lines = choose_best_chunk_option(i, reads, chunks, threads, circular)

        new_best_seqs[i] = best_seq
        if chunk.best_seq == best_seq:
            kept += 1
        else:
            changed += 1

        if verbose:
            for line in output_lines:
                log(line)
            log()
        else:
            completed += 1
            log(f'\rProcessing chunks: {completed:,} / {needs_assessment_count:,}', end='')

    for i, chunk in enumerate(chunks):
        if not chunk.needs_assessment:
            continue
        assert chunk.type == 'different'
        chunk.best_seq = new_best_seqs[i]

    if not verbose:
        log('\n')
    log('Chunks where sequence is...')
    log(f'  the same as in the initial consensus: {kept:,}')
    log(f'  different to the initial consensus:   {changed:,}')
    log()


def choose_best_chunk_option(i, reads, chunks, threads, circular):
    chunk = chunks[i]
    assert chunk.type == 'different'
    assert len(chunk.seqs) > 1

    chunk_reads = [reads[r] for r in sorted(chunk.read_names)]

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = pathlib.Path(temp_dir)

        # Save the reads relevant to this chunk to a file.
        fastq_filename = temp_dir / 'chunk_reads.fastq'
        with open(fastq_filename, 'wt') as fastq_file:
            for header, seq, qual in chunk_reads:
                fastq_file.write(f'@{header}\n{seq}\n+\n{qual}\n')

        scores = {}
        all_seqs = sorted(''.join(s) for s in chunk.seqs.values())
        option_seqs = sorted(set(all_seqs))  # get rid of duplicates
        distances = get_hamming_totals(all_seqs, option_seqs)

        counts = collections.defaultdict(int)
        for s in all_seqs:
            counts[s] += 1

        # Get the alignment scores for each option.
        for option_seq in option_seqs:
            if option_seq in chunk.assessment_options:
                test_sequence = build_test_sequence(i, chunks, option_seq, circular,
                                                    settings.CHUNK_TEST_MARGIN)
                alignments = align_reads_to_seq(fastq_filename, test_sequence, threads,
                                                scores=(1, 1, 1, 1))
                alignments = get_best_alignment_per_read(alignments)
                score_sum = sum(a.alignment_score for a in alignments)
                scores[option_seq] = score_sum
            else:
                scores[option_seq] = 0

    best_seq = max(scores, key=scores.get)

    output_lines = [f'chunk {i + 1}']
    for seq in option_seqs:
        line = f'  {seq}: count = {counts[seq]}, Hamming total = {distances[seq]}'
        if scores[seq] > 0:
            line += f', score = {scores[seq]:,}'
        else:
            line += ', not assessed'
        if seq == chunk.best_seq:
            line += ', initial'
        if seq == best_seq:
            line += ', best'
        output_lines.append(line)

    return best_seq, output_lines


def build_test_sequence(i, chunks, option_seq, circular, chunk_test_margin):
    """
    Builds a piece of the consensus sequence with the to-be-tested option in the middle.
    """
    chunk_count = len(chunks)
    assert chunks[i].type == 'different'
    if chunk_test_margin == 0:
        return option_seq.replace('-', '')

    # Set some limits for circular sequences so the option seq is in the middle.
    total_length_with_gaps = 0
    for j, c in enumerate(chunks):
        if i == j:
            total_length_with_gaps += len(option_seq)
        else:  # for most chunks, we just use the initial consensus option
            total_length_with_gaps += len(c.best_seq)
    length_minus_option = total_length_with_gaps - len(option_seq)
    max_forward = math.ceil(length_minus_option / 2)
    max_backward = math.floor(length_minus_option / 2)

    # Build forward.
    j = i
    forward_length = 0
    forward_sequence = []
    while True:
        j += 1
        if j == chunk_count:
            if not circular:
                break
            else:
                j = 0
        new_seq = chunks[j].best_seq
        forward_sequence.append(new_seq)
        forward_length += len(new_seq)
        if forward_length >= chunk_test_margin:
            break
        if circular and forward_length >= max_forward:
            break
    forward_sequence = ''.join(forward_sequence)
    forward_sequence = forward_sequence[:chunk_test_margin]
    if circular:
        forward_sequence = forward_sequence[:max_forward]

    # Build backward.
    j = i
    backward_length = 0
    backward_sequence = []
    while True:
        j -= 1
        if j == -1:
            if not circular:
                break
            else:
                j = chunk_count - 1
        new_seq = chunks[j].best_seq
        backward_sequence.insert(0, new_seq)  # add to front
        backward_length += len(new_seq)
        if backward_length >= chunk_test_margin:
            break
    backward_sequence = ''.join(backward_sequence)
    backward_sequence = backward_sequence[-chunk_test_margin:]
    if circular:
        backward_sequence = backward_sequence[-max_backward:]

    return ''.join([backward_sequence, option_seq, forward_sequence]).replace('-', '')


def sanity_check_chunks(chunks, msa_length):
    """
    Makes sure that chunks alternate in type: same, different, same, different, etc.
    """
    total_length = 0
    for i, chunk in enumerate(chunks):
        assert chunk.type is not None
        if i > 0:
            prev_chunk = chunks[i-1]
            if chunk.type == 'same':
                assert prev_chunk.type == 'different'
            elif chunk.type == 'different':
                assert prev_chunk.type == 'same'
            else:
                assert False
        total_length += chunk.get_length()
    assert total_length == msa_length


def combine_chunks(chunks, combine_size):
    """
    This function combines chunks when there is a very small 'same' chunk between two 'different'
    chunks.
    """
    log('combining small chunks: ', end='')
    combined_chunks = []
    for i, chunk in enumerate(chunks):
        if i == 0 or chunk.type == 'different':
            combined_chunks.append(chunk)
        else:
            assert chunk.type == 'same'
            assert combined_chunks[-1].type == 'different'
            if chunk.get_length() <= combine_size:
                combined_chunks[-1].add_one_seq_to_seqs(chunk.seq)
            else:
                combined_chunks.append(chunk)

    # We are now in a position where two adjacent chunks might both be 'different' chunks, so
    # we need to merge them together.
    new_chunks = []
    for i, chunk in enumerate(combined_chunks):
        if i == 0 or chunk.type == 'same':
            new_chunks.append(chunk)
        else:
            assert chunk.type == 'different'
            if new_chunks[-1].type == 'different':
                new_chunks[-1].add_multiple_seqs_to_seqs(chunk.seqs)
            else:
                new_chunks.append(chunk)

    same_count = len([c for c in new_chunks if c.type == 'same'])
    different_count = len([c for c in new_chunks if c.type == 'different'])
    log(f'{len(new_chunks):,} ({same_count:,} same, {different_count:,} different)')
    return new_chunks


class Chunk(object):
    """
    This class holds a chunk of the MSA, which can either be a 'same' chunk (where all sequences
    agree) or a 'different' chunk (where there is at least one difference).
    """
    def __init__(self):
        self.type = None  # will be either 'same' or 'different'
        self.seq = None  # will hold the sequence for a 'same' chunk
        self.seqs = None  # will hold the multiple alternative sequences for a 'different' chunk
        self.read_names = set()  # will hold read names relevant for assessing this chunk
        self.best_seq = None  # will hold the chunk's best sequence as a string
        self.had_tie = False  # whether or not this chunk's consensus sequence involved a tie

        self.needs_assessment = None
        self.assessment_options = None

    def add_bases(self, bases):
        assert self.can_add_bases(bases)

        # If this is a new chunk, we'll set its type now.
        if self.type is None:
            base_count = len(set(bases.values()))
            if base_count == 1:
                self.type = 'same'
                self.seq = []
            else:
                self.type = 'different'
                self.seqs = {n: [] for n in bases.keys()}

        if self.type == 'same':
            base = list(bases.values())[0]
            self.seq.append(base)

        else:
            assert self.type == 'different'
            for name, base in bases.items():
                self.seqs[name].append(base)

    def can_add_bases(self, bases):
        """
        Tests to see whether the given bases are incompatible with the current chunk type.
        """
        if self.type is None:
            return True
        base_count = len(set(bases.values()))
        if self.type == 'same' and base_count == 1:
            return True
        if self.type == 'different' and base_count > 1:
            return True
        return False

    def is_empty(self):
        return self.type is None

    def get_length(self):
        if self.type is None:
            return 0
        elif self.type == 'same':
            return len(self.seq)
        elif self.type == 'different':
            lengths = set(len(seq) for seq in self.seqs.values())
            assert len(lengths) == 1  # all seqs should be the same length
            return list(lengths)[0]
        else:
            assert False

    def __str__(self):
        if self.type is None:
            return ''
        elif self.type == 'same':
            return ''.join(self.seq)
        elif self.type == 'different':
            seq_lines = []
            longest_name = max(len(name) for name in self.seqs.keys())
            for name, seq in self.seqs.items():
                seq = ''.join(seq)
                seq_lines.append(f'{name.rjust(longest_name)}: {seq}')
            return '\n'.join(seq_lines)

    def add_one_seq_to_seqs(self, additional_seq):
        assert self.type == 'different'
        new_seqs = {}
        for name, seq in self.seqs.items():
            new_seqs[name] = seq + additional_seq
        self.seqs = new_seqs

    def add_multiple_seqs_to_seqs(self, additional_seqs):
        assert self.type == 'different'
        new_seqs = {}
        for name, seq in self.seqs.items():
            new_seqs[name] = seq + additional_seqs[name]
        self.seqs = new_seqs

    def prepare_chunk(self):
        """
        This function is run after the chunk has had all of its sequences added. It:
          * sets an initial best sequence (though this may change later with read-based assessment)
          * sets whether or not read-based assessment is needed
          * sets which of the option sequences will be considered during read-bases assessment
        """
        assert self.type is not None
        if self.type == 'same':
            self.prepare_same_chunk()
        elif self.type == 'different':
            self.prepare_different_chunk()

    def prepare_same_chunk(self):
        assert self.type == 'same'
        self.best_seq = ''.join(self.seq)
        self.needs_assessment = False
        self.assessment_options = None

    def prepare_different_chunk(self):
        assert self.type == 'different'
        all_options = [''.join(seq) for seq in self.seqs.values()]
        unique_options = list(set(all_options))

        hamming_distances = get_hamming_totals(all_options, unique_options)
        hamming_distances = sorted(hamming_distances.items(), key=lambda x: x[1])
        best_hamming_distance = hamming_distances[0][1]
        best_seqs = [x[0] for x in hamming_distances if x[1] == best_hamming_distance]

        if len(best_seqs) == 1:  # a clear winner
            self.best_seq = best_seqs[0]
            self.needs_assessment = False
            self.assessment_options = None

        else:  # a tie
            self.needs_assessment = True
            self.assessment_options = best_seqs

            # Try to break the tie using counts, if possible.
            counts = {seq: count for seq, count in collections.Counter(all_options).items()
                      if seq in best_seqs}
            best_count = max(counts.values())
            best_seqs = [seq for seq in best_seqs if counts[seq] == best_count]
            self.best_seq = sorted(best_seqs)[0]  # lexicographically first

        # If there was only one instance of each option, then we will assess the chunk using all
        # sequences.
        if len(all_options) == len(unique_options):
            self.needs_assessment = True
            self.assessment_options = all_options


def get_hamming_totals(all_options, unique_options):
    hamming_distances = {x: 0 for x in unique_options}
    for x in unique_options:
        for y in all_options:
            hamming_distances[x] += hamming_distance(x, y)
    return hamming_distances


def hamming_distance(s1, s2):
    dist = 0
    assert len(s1) == len(s2)
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            dist += 1
    return dist


def check_seqs(cluster_dir):
    filename = cluster_dir / '2_all_seqs.fasta'
    log(f'Input contigs: {filename}')
    contig_type = get_sequence_file_type(filename)
    if contig_type != 'FASTA':
        sys.exit(f'\nError: input contig file ({filename}) is not in FASTA format')
    seqs = load_fasta(filename)
    if len(seqs) == 0:
        sys.exit(f'\nError: input contig file ({filename}) contains no sequences')
    contig_names = set()
    for contig_name, seq in seqs:
        if contig_name in contig_names:
            sys.exit(f'\nError: duplicate contig name: {contig_name}')
        contig_names.add(contig_name)
        log(f'  {contig_name}: {len(seq):,} bp')
    log()


def check_cluster_directory(directory):
    if directory.is_file():
        sys.exit(f'\nError: output directory ({directory}) already exists as a file')
    if not directory.is_dir():
        sys.exit(f'\nError: output directory ({directory}) does not exist')

    seq_file = directory / '2_all_seqs.fasta'
    if not seq_file.is_file():
        sys.exit(f'\nError: output directory ({directory}) does not contain2_all_seqs.fasta')

    pairwise_file = directory / '3_msa.fasta'
    if not pairwise_file.is_file():
        sys.exit(f'\nError: output directory ({directory}) does not contain 3_msa.fasta')

    reads_file = directory / '4_reads.fastq'
    if not reads_file.is_file():
        sys.exit(f'\nError: output directory ({directory}) does not contain 4_reads.fastq')


def check_required_software():
    log('Checking required software:')
    check_minimap2()
    log()


def load_seqs(cluster_dir):
    filename = cluster_dir / '2_all_seqs.fasta'
    seqs = dict(load_fasta(filename))
    seq_names = sorted(seqs.keys())
    seq_lengths = {name: len(seq) for name, seq in seqs.items()}
    return seqs, seq_names, seq_lengths


def load_msa(cluster_dir):
    filename = cluster_dir / '3_msa.fasta'
    seqs = dict(load_fasta(filename))
    seq_names = sorted(seqs.keys())
    seq_lengths = {name: len(seq) for name, seq in seqs.items()}

    seq_length_set = set(seq_lengths.values())
    assert len(seq_length_set) == 1
    msa_length = list(seq_length_set)[0]

    return seqs, seq_names, msa_length


def sanity_check_msa(seqs, seq_names, seq_lengths, msa_seqs, msa_names, msa_length):
    assert seq_names == msa_names
    for n in seq_names:
        assert seq_lengths[n] <= msa_length
        assert seqs[n] == msa_seqs[n].replace('-', '')


def save_seqs_to_fasta(seqs, filename, extra_newline=True):
    seq_word = 'sequence' if len(seqs) == 1 else 'sequences'
    log(f'Saving {seq_word} to file: {filename}')
    with open(filename, 'wt') as fasta:
        for name, seq in seqs.items():
            fasta.write(f'>{name}\n')
            fasta.write(f'{seq}\n')
    if extra_newline:
        log()


def save_chunks_to_gfa(chunks, filename, input_count, circular, extra_newline=True):
    chunk_word = 'sequence' if len(chunks) == 1 else 'sequences'
    log(f'Saving {chunk_word} to graph: {filename}')
    with open(filename, 'wt') as gfa:
        gfa.write('H\tVN:Z:1.0\tbn:Z:--linear --singlearr\n')  # header line with Bandage options
        link_lines = []
        prev_chunk_names = None
        first_chunk_names, last_chunk_names = [], []
        for i, chunk in enumerate(chunks):
            if chunk.type == 'same':
                assert chunk.seq is not None
                chunk_seq = ''.join(chunk.seq)
                chunk_name = str(i+1)
                if i == 0:
                    first_chunk_names.append(chunk_name)
                if i == len(chunks) - 1:
                    last_chunk_names.append(chunk_name)
                gfa.write(f'S\t{chunk_name}\t{chunk_seq}\tdp:f:{input_count}\n')
                if prev_chunk_names is not None:
                    assert len(prev_chunk_names) > 1  # same chunks are preceded by diff chunks
                    for prev_name in prev_chunk_names:
                        link_lines.append(f'L\t{prev_name}\t+\t{chunk_name}\t+\t0M\n')
                prev_chunk_names = [chunk_name]

            elif chunk.type == 'different':
                assert chunk.seqs is not None
                chunk_seq_counts = collections.defaultdict(int)
                chunk_names = []
                for s in chunk.seqs.values():
                    chunk_seq_counts[''.join(s)] += 1
                j = 1
                for chunk_seq, count in chunk_seq_counts.items():
                    chunk_name = f'{i+1}_{j}'
                    chunk_names.append(chunk_name)
                    if i == 0:
                        first_chunk_names.append(chunk_name)
                    if i == len(chunks) - 1:
                        last_chunk_names.append(chunk_name)
                    gfa.write(f'S\t{chunk_name}\t{chunk_seq}\tdp:f:{count}\n')
                    j += 1
                if prev_chunk_names is not None:
                    assert len(prev_chunk_names) == 1  # diff chunks are preceded by same chunks
                    prev_name = prev_chunk_names[0]
                    for chunk_name in chunk_names:
                        link_lines.append(f'L\t{prev_name}\t+\t{chunk_name}\t+\t0M\n')
                prev_chunk_names = chunk_names
            else:
                assert False
        if circular:
            for first_chunk_name in first_chunk_names:
                for last_chunk_name in last_chunk_names:
                    link_lines.append(f'L\t{last_chunk_name}\t+\t{first_chunk_name}\t+\t0M\n')
        for link_line in link_lines:
            gfa.write(link_line)

    if extra_newline:
        log()
