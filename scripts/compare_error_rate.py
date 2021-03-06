#!/usr/bin/env python
#
# ############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################
import pysam
from Bio import Seq, SeqIO
import numpy as np
import sys
sys.path.append("../")
from gtf2db import *
from src.junction_comparator import *
from src.common import *


logger = logging.getLogger('CSA')


class ErrorRateStat:
    def __init__(self, name=""):
        self.name = name
        self.mismatches = np.zeros(100)
        self.insertions = np.zeros(100)
        self.deletions = np.zeros(100)
        self.total = np.zeros(100)
        self.mismapped = 0
        self.unmapped = 0
        self.no_poly = 0

    def dump(self, outf):
        outf.write("##%s\n" % self.name)
        outf.write("#unmapped\tmismapped\tnopolya\n")
        outf.write("%d\t%d\t%d\t\n" % (self.unmapped, self.mismapped, self.no_poly))
        total_nucls = sum(self.total)
        mm_rate = sum(self.mismatches) / total_nucls
        d_rate = sum(self.deletions) / total_nucls
        i_rate = sum(self.insertions) / total_nucls
        outf.write("#total\tmismatch\tdeletion\tinsertion\n")
        outf.write("%.6f\t%.6f\t%.6f\t%.6f\n" % (mm_rate + i_rate + d_rate, mm_rate, d_rate, i_rate))
        outf.write("#position\ttotal\tX\tI\tD\n")
        for i in range(len(self.total)):
            outf.write("%d\t%d\t%d\t%d\t%d\n" %
                       (i, self.total[i], self.mismatches[i], self.insertions[i], self.deletions[i]))

    def add_unmapped(self):
        self.unmapped += 1

    def add_mismapped(self):
        self.mismapped += 1

    def add_no_polya(self):
        self.no_poly += 1

    def get_index(self, position, read_len, reverse):
        if reverse:
            return math.floor(float(read_len-position-1) * 100 / read_len)
        return math.floor(float(position) * 100 / read_len)

    def add_match(self, position, read_len, reverse):
        index = self.get_index(position, read_len, reverse)
        self.total[index] += 1

    def add_mismatch(self, position, read_len, reverse):
        index = self.get_index(position, read_len, reverse)
        self.total[index] += 1
        self.mismatches[index] += 1

    def add_insertion(self, position, read_len, reverse):
        index = self.get_index(position, read_len, reverse)
        self.total[index] += 1
        self.insertions[index] += 1

    def add_deletion(self, position, read_len, reverse):
        index = self.get_index(position, read_len, reverse)
        self.total[index] += 1
        self.deletions[index] += 1


# [start, end]
def check_homopolymer(ref_seq, start, end):
    if start < 0 or end >= len(ref_seq):
        return False, None
    seq = ref_seq[start:end + 1]
    c = seq[0]
    if all(x == c for x in seq):
        return True, c
    return False, None


def check_homopolymer_mismatch(read_seq, read_pos, ref_seq, ref_pos):
    read_nucl = read_seq[read_pos]
    for i in range(3):
        is_homopolymer, c = check_homopolymer(ref_seq, ref_pos-2+i, ref_pos+i)
        if is_homopolymer and read_nucl != c:
            return True
    return False


def check_homopolymer_deletion(ref_seq, ref_pos):
    for i in range(3):
        is_homopolymer, _ = check_homopolymer(ref_seq, ref_pos-2+i, ref_pos+i)
        if is_homopolymer:
            return True
    return False


def check_homopolymer_insertion(read_seq, read_pos, ref_seq, ref_pos):
    read_nucl = read_seq[read_pos]
    for i in range(4):
        is_homopolymer, c = check_homopolymer(ref_seq, ref_pos-3+i, ref_pos-1+i)
        if is_homopolymer and read_nucl == c:
            return True
    return False

    
def check_polya(sequence):
    len_slice = 12
    n = len(sequence)
    if n < len_slice:
        return -1
    for i in range(n - len_slice):
        slice_ = sequence[i: i + len_slice]
        k = slice_.count('A')
        if (k >= 0.75 * len(slice_)):
            return i + slice_.find("AA")
    return -1


def get_sequence_to_check(alignment):
    cigar_tuples = alignment.cigartuples
    clipped_len = -1
    sequence_to_check = ''
    if len(cigar_tuples) > 1 and cigar_tuples[0][0] == 5 and cigar_tuples[1][0] == 4:
        # hard clipped
        clipped_len = cigar_tuples[1][1]
    elif cigar_tuples[0][0] == 4:
        # soft clipped
        clipped_len = cigar_tuples[0][1]
    if (clipped_len != -1):
        sequence_to_check = alignment.seq[:clipped_len]  # str(Seq.Seq(alignment.seq[:clipped_len]).reverse_complement()).upper()

    clipped_len = -1
    if len(cigar_tuples) > 1 and cigar_tuples[-1][0] == 5 and cigar_tuples[-2][0] == 4:
        # hard clipped
        clipped_len = cigar_tuples[-2][1]
    elif cigar_tuples[-1][0] == 4:
        # soft clipped
        clipped_len = cigar_tuples[-1][1]

    sequence_to_check_end = ''
    if (clipped_len != -1):
        sequence_to_check_end = str((alignment.seq[-clipped_len:]).upper())

    return sequence_to_check, sequence_to_check_end


def get_direction(alignment):
    sequence_to_check_start, sequence_to_check_end = get_sequence_to_check(alignment)
    sequence_to_check_start = str(Seq.Seq(sequence_to_check_start).reverse_complement()).upper()
    if check_polya(sequence_to_check_start) != -1:
        return 1
    elif check_polya(sequence_to_check_end) != -1:
        return 0
    return -1


def clip_seq(seq, cigar):
    start_clip = 0
    if cigar[0][0] == 4:
        start_clip = cigar[0][1]
    end_clip = 0
    if cigar[-1][0] == 4:
        end_clip = cigar[-1][1]
    return seq[start_clip:len(seq)-end_clip]


def increment_position(cigar_event, ref_index, read_index):
    if cigar_event == 0:
        ref_index += 1
        read_index += 1
    elif cigar_event == 1:
        read_index += 1
    elif cigar_event in [2,3]:
        ref_index += 1
    return ref_index, read_index


def process_alignment_pair(alignment_record1, alignment_record2, fasta_records, stats1, stats2,
                           homopolymer_stats1=None, homopolymer_stats2=None):
    if alignment_record1.is_secondary or alignment_record1.is_supplementary or \
            alignment_record2.is_secondary or alignment_record2.is_supplementary or \
            alignment_record1.seq is None or alignment_record2.seq is None:
        return
    if alignment_record1.reference_id == -1 or alignment_record2.reference_id == -1:
        stats1.add_unmapped_pair()
        stats2.add_unmapped_pair()
        return
    if alignment_record1.reference_id != alignment_record2.reference_id:
        stats1.add_mismapped_pair()
        stats2.add_mismapped_pair()
        return

    ref_start1 = alignment_record1.reference_start
    ref_length1 = alignment_record1.reference_length
    ref_end1 = ref_start1 + ref_length1
    cigar1 = alignment_record1.cigartuples
    ref_start2 = alignment_record2.reference_start
    ref_length2 = alignment_record2.reference_length
    ref_end2 = ref_start2 + ref_length2
    cigar2 = alignment_record2.cigartuples

    if ref_length1 is None or cigar1 is None or ref_length2 is None or cigar2 is None:
        stats1.add_unmapped_pair()
        stats2.add_unmapped_pair()
        return

    ref_region_start = min(ref_start1, ref_start2)
    ref_region_end = max(ref_end1, ref_end2)
    true_fragment = fasta_records[alignment_record1.reference_name][ref_region_start:ref_region_end]
    read_shift1 = ref_start1 - ref_region_start
    read_shift2 = ref_start2 - ref_region_start
    assert read_shift2 == 0 or read_shift1 == 0

    direction1 = get_direction(alignment_record1)
    if (direction1 == -1):
        stats1.add_no_polya()
        return
    elif (direction1 == 1):
        reverse_read1 = True
    else:
        reverse_read1 = False

    direction2 = get_direction(alignment_record2)
    if (direction2 == -1):
        stats2.add_no_polya()
        return
    elif (direction2 == 1):
        reverse_read2 = True
    else:
        reverse_read2 = False

    seq1 = clip_seq(alignment_record1.seq, cigar1)
    seq2 = clip_seq(alignment_record2.seq, cigar2)

    current_cigar_index1 = 0
    current_cigar_index2 = 0
    # skip clipped
    while cigar1[current_cigar_index1][0] in [4, 5]:
        current_cigar_index1 += 1
    while cigar2[current_cigar_index2][0] in [4, 5]:
        current_cigar_index2 += 1

    cigar1_consumed = 0
    cigar2_consumed = 0
    read_index1 = 0
    read_index2 = 0

    ref_index = 0
    while ref_index < len(true_fragment):
        if read_index1 >= len(seq1) or read_index2 >= len(seq2) or \
                current_cigar_index1 >= len(cigar1) or current_cigar_index2 >= len(cigar2):
            break

        event1 = cigar1[current_cigar_index1][0]
        event2 = cigar2[current_cigar_index2][0]

        if event1 in [4,5] or event2 in [4,5]:
            logger.warning("Reached right-side clipping event")
            logger.info("Read index 1: %d, RL: %d, cigar index: %d, cigar tuples %s" %
                        (read_index1, len(seq1), current_cigar_index1, str(cigar1)))
            logger.info("Read index 2: %d, RL: %d, cigar index: %d, cigar tuples %s" %
                        (read_index2, len(seq2), current_cigar_index2, str(cigar2)))
            break

        if ref_index < read_shift1:
            # read 1 is not started yet
            ref_index, read_index2 = increment_position(event2, ref_index, read_index2)
            cigar2_consumed += 1

        elif ref_index < read_shift2:
            # read 2 is not started yet
            ref_index, read_index1 = increment_position(event1, ref_index, read_index1)
            cigar1_consumed += 1

        elif event1 == event2:
            cigar1_consumed += 1
            cigar2_consumed += 1
            if event1 == 0:
                # M
                if seq1[read_index1] != seq2[read_index2]:
                    # mismatch, compare with the reference
                    if seq1[read_index1] != true_fragment[ref_index]:
                        stats1.add_mismatch(read_index1, len(seq1), reverse_read1)
                        if homopolymer_stats1 and check_homopolymer_mismatch(seq1, read_index1, true_fragment, ref_index):
                            homopolymer_stats1.add_mismatch(read_index1, len(seq1), reverse_read1)
                    else:
                        stats1.add_match(read_index1, len(seq1), reverse_read1)
                    if seq2[read_index2] != true_fragment[ref_index]:
                        stats2.add_mismatch(read_index2, len(seq2), reverse_read2)
                        if homopolymer_stats2 and check_homopolymer_mismatch(seq2, read_index2, true_fragment, ref_index):
                            homopolymer_stats2.add_mismatch(read_index2, len(seq2), reverse_read2)
                    else:
                        stats2.add_match(read_index2, len(seq2), reverse_read2)
                else:
                    # same letter
                    stats1.add_match(read_index1, len(seq1), reverse_read1)
                    stats2.add_match(read_index2, len(seq2), reverse_read2)
                ref_index += 1
                read_index1 += 1
                read_index2 += 1
            elif event1 == 1:
                # I
                if seq1[read_index1] == seq2[read_index2]:
                    # same
                    stats1.add_match(read_index1, len(seq1), reverse_read1)
                    stats2.add_match(read_index2, len(seq2), reverse_read2)
                else:
                    # different
                    stats1.add_insertion(read_index1, len(seq1), reverse_read1)
                    stats2.add_insertion(read_index2, len(seq2), reverse_read2)
                    if homopolymer_stats1 and check_homopolymer_insertion(seq1, read_index1, true_fragment, ref_index):
                        homopolymer_stats1.add_insertion(read_index1, len(seq1), reverse_read1)
                    if homopolymer_stats2 and check_homopolymer_insertion(seq2, read_index2, true_fragment, ref_index):
                        homopolymer_stats2.add_insertion(read_index2, len(seq2), reverse_read2)
                read_index1 += 1
                read_index2 += 1
            elif event1 == 2:
                # D
                stats1.add_match(read_index1, len(seq1), reverse_read1)
                stats2.add_match(read_index2, len(seq2), reverse_read2)
                ref_index += 1
            elif event1 == 3:
                # N
                ref_index += 1

        else:
            # non equal events
            if event1 == 1:
                # I in read1
                stats1.add_insertion(read_index1, len(seq1), reverse_read1)
                if homopolymer_stats1 and check_homopolymer_insertion(seq1, read_index1, true_fragment, ref_index):
                    homopolymer_stats1.add_insertion(read_index1, len(seq1), reverse_read1)
                cigar1_consumed += 1
                read_index1 += 1
            elif event2 == 1:
                # I in read2
                stats2.add_insertion(read_index2, len(seq2), reverse_read2)
                if homopolymer_stats2 and check_homopolymer_insertion(seq2, read_index2, true_fragment, ref_index):
                    homopolymer_stats2.add_insertion(read_index2, len(seq2), reverse_read2)
                cigar2_consumed += 1
                read_index2 += 1
            elif event1 == 3:
                # N in read 1
                assert event2 != 1 # insertions should be processed earlier
                ref_index, read_index2 = increment_position(event2, ref_index, read_index2)
                # consume intronic base only if reference is shifted
                cigar1_consumed += 1
                cigar2_consumed += 1
            elif event2 == 3:
                assert event1 != 1 # insertions should be processed earlier
                # N in read 2
                ref_index, read_index1 = increment_position(event1, ref_index, read_index1)
                cigar1_consumed += 1
                cigar2_consumed += 1
            elif event1 == 0:
                # read 1 - match, read 2 - deletion
                assert event2 == 2
                if seq1[read_index1] != true_fragment[ref_index]:
                    stats1.add_mismatch(read_index1, len(seq1), reverse_read1)
                else:
                    stats1.add_match(read_index1, len(seq1), reverse_read1)
                stats2.add_deletion(read_index2, len(seq2), reverse_read2)
                if homopolymer_stats2 and check_homopolymer_deletion(true_fragment, ref_index):
                    homopolymer_stats2.add_deletion(read_index2, len(seq2), reverse_read2)
                read_index1 += 1
                ref_index += 1
                cigar1_consumed += 1
                cigar2_consumed += 1
            elif event2 == 0:
                # read 1 - deletion, read 2 - match
                assert event1 == 2
                if seq2[read_index2] != true_fragment[ref_index]:
                    stats2.add_mismatch(read_index2, len(seq2), reverse_read2)
                else:
                    stats2.add_match(read_index2, len(seq2), reverse_read2)
                stats1.add_deletion(read_index1, len(seq1), reverse_read1)
                if homopolymer_stats1 and check_homopolymer_deletion(true_fragment, ref_index):
                    homopolymer_stats1.add_deletion(read_index1, len(seq1), reverse_read1)
                read_index2 += 1
                ref_index += 1
                cigar1_consumed += 1
                cigar2_consumed += 1
            else:
                logger.error("Unexpected pair of events: ")
                cigar1_consumed += 1
                cigar2_consumed += 1

        if cigar1_consumed >= cigar1[current_cigar_index1][1]:
            current_cigar_index1 += 1
            cigar1_consumed = 0
        if cigar2_consumed >= cigar2[current_cigar_index2][1]:
            current_cigar_index2 += 1
            cigar2_consumed = 0

    if (read_index1 < len(seq1) and read_index2 < len(seq2)) or \
            (current_cigar_index1 < len(cigar1) and cigar1[current_cigar_index1][0] < 4 and
             current_cigar_index2 < len(cigar2) and cigar2[current_cigar_index2][0] < 4):
        logger.warning("Both reads is processed incompletely")
        logger.info("Read: %d / %d, ref: %d / % d, cigar %d (%d) / %d" %
                    (read_index1, len(seq1), ref_index, len(true_fragment), current_cigar_index1,
                     cigar1[current_cigar_index1][0], len(cigar1)))
        logger.info("Read: %d / %d, ref: %d / % d, cigar %d (%d) / %d" %
                    (read_index2, len(seq2), ref_index, len(true_fragment), current_cigar_index2,
                     cigar1[current_cigar_index2][0], len(cigar2)))


def process_alignment(alignment_record1, fasta_records, stats1, homopolymer_stats1=None):
    if alignment_record1.is_secondary or alignment_record1.is_supplementary or alignment_record1.seq is None:
        return
    if alignment_record1.reference_id == -1:
        stats1.add_unmapped_pair()
        return

    ref_start1 = alignment_record1.reference_start
    ref_length1 = alignment_record1.reference_length
    ref_end1 = ref_start1 + ref_length1
    cigar1 = alignment_record1.cigartuples

    if ref_length1 is None or cigar1 is None:
        stats1.add_unmapped_pair()
        return

    ref_region_start = ref_start1
    ref_region_end = ref_end1
    true_fragment = fasta_records[alignment_record1.reference_name][ref_region_start:ref_region_end]

    direction1 = get_direction(alignment_record1)
    if (direction1 == -1):
        #stats1.add_no_polya()
        #return
        reverse_read1 = False
    elif (direction1 == 1):
        reverse_read1 = True
    else:
        reverse_read1 = False

    seq1 = clip_seq(alignment_record1.seq, cigar1)
    current_cigar_index1 = 0
    # skip clipped
    while cigar1[current_cigar_index1][0] in [4, 5]:
        current_cigar_index1 += 1

    cigar1_consumed = 0
    read_index1 = 0
    ref_index = 0
    while ref_index < len(true_fragment) and read_index1 < len(seq1) and current_cigar_index1 < len(cigar1):
        event1 = cigar1[current_cigar_index1][0]
        if event1 in [4,5]:
            logger.warning("Reached right-side clipping event")
            logger.info("Read index 1: %d, RL: %d, cigar index: %d, cigar tuples %s" %
                        (read_index1, len(seq1), current_cigar_index1, str(cigar1)))
            break

        cigar1_consumed += 1
        if event1 == 0:
            # M
            if seq1[read_index1] != true_fragment[ref_index]:
                stats1.add_mismatch(read_index1, len(seq1), reverse_read1)
                if homopolymer_stats1 and check_homopolymer_mismatch(seq1, read_index1, true_fragment, ref_index):
                    homopolymer_stats1.add_mismatch(read_index1, len(seq1), reverse_read1)
            else:
                stats1.add_match(read_index1, len(seq1), reverse_read1)

            ref_index += 1
            read_index1 += 1
        elif event1 == 1:
            # I
            stats1.add_insertion(read_index1, len(seq1), reverse_read1)
            if homopolymer_stats1 and check_homopolymer_insertion(seq1, read_index1, true_fragment, ref_index):
                homopolymer_stats1.add_insertion(read_index1, len(seq1), reverse_read1)
            read_index1 += 1
        elif event1 == 2:
            # D
            stats1.add_deletion(read_index1, len(seq1), reverse_read1)
            if homopolymer_stats1 and check_homopolymer_deletion(true_fragment, ref_index):
                homopolymer_stats1.add_deletion(read_index1, len(seq1), reverse_read1)
            ref_index += 1
        elif event1 == 3:
            # N
            ref_index += 1

        if cigar1_consumed >= cigar1[current_cigar_index1][1]:
            current_cigar_index1 += 1
            cigar1_consumed = 0

    if read_index1 < len(seq1) or ref_index < len(true_fragment) or \
            (current_cigar_index1 < len(cigar1) and cigar1[current_cigar_index1][0] < 4):
        logger.warning("Read is processed incompletely")
        logger.info("Read: %d / %d, ref: %d / % d, cigar %d (%d) / %d" %
                    (read_index1, len(seq1), ref_index, len(true_fragment), current_cigar_index1,
                     cigar1[current_cigar_index1][0], len(cigar1)))


def error_rate_stats(read_pairs, bam_records1, bam_records2, chr_records):
    counter = 0
    stats1 = ErrorRateStat("Reads 1 3-way comparison")
    stats2 = ErrorRateStat("Reads 2 3-way comparison")
    hstats1 = ErrorRateStat("Reads 1 3-way hompolymer comparison")
    hstats2 = ErrorRateStat("Reads 2 3-way hompolymer comparison")
    stats1_only = ErrorRateStat("Reads 1 simple comparison")
    hstats1_only = ErrorRateStat("Reads 1 simple hompolymer comparison")
    stats2_only = ErrorRateStat("Reads 2 simple comparison")
    hstats2_only = ErrorRateStat("Reads 2 simple hompolymer comparison")

    for read_pair in read_pairs:
        if read_pair[0] not in bam_records1 or read_pair[1] not in bam_records2:
            continue
        bam_record1 = bam_records1[read_pair[0]]
        bam_record2 = bam_records2[read_pair[1]]
        if bam_record1.reference_name != bam_record2.reference_name:
            continue

        process_alignment_pair(bam_record1, bam_record2, chr_records, stats1, stats2, hstats1, hstats2)
        process_alignment(bam_record1, chr_records, stats1_only, hstats1_only)
        process_alignment(bam_record2, chr_records, stats2_only, hstats2_only)

        counter += 1
        if counter % 1000 == 0:
            logger.info("Processed %d read pairs (%0.1f%%)" % (counter, 100 * counter / len(read_pairs)))

    return [("paired_stats", [stats1, stats2]), ("paired_homopolymer_stats", [hstats1, hstats2]),
            ("reads1_single_stats", [stats1_only, hstats1_only]), ("reads2_single_stats", [stats2_only, hstats2_only])]


def error_rate_single(bam_records1, chr_records):
    counter = 0
    stats1_only = ErrorRateStat("Reads 1 simple comparison")
    hstats1_only = ErrorRateStat("Reads 1 simple hompolymer comparison")

    for read_id, bam_record1 in bam_records1.items():
        process_alignment(bam_record1, chr_records, stats1_only)

        counter += 1
        if counter % 1000 == 0:
            logger.info("Processed %d reads (%0.1f%%)" % (counter, 100 * counter / len(bam_records1)))

    return [("reads1_single_stats", [stats1_only])]


def load_bam(read_set, bamfile):
    bam_records = {}
    for r in pysam.AlignmentFile(bamfile, "rb").fetch():
        if not read_set or r.query_name in read_set:
            bam_records[r.query_name] = r
    return bam_records


def load_tsv(tsv_file):
    read_pairs = []
    for l in open(tsv_file):
        if l.startswith("#"):
            continue
        t = l.strip().split()
        if len(t) < 4:
            continue
        read_pairs.append((t[2], t[3]))
    return read_pairs


def parse_args(args=None, namespace=None):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("--output", "-o", help="output prefix", type=str)
    # REFERENCE
    parser.add_argument("--reference", "-r", help="reference genome in FASTA format", type=str)

    parser.add_argument('--bam_pb', type=str, help='sorted and indexed BAM file for PacBio')
    parser.add_argument('--bam_ont', type=str, help='sorted and indexed BAM file for ONT')
    parser.add_argument('--tsv', type=str, help='TSV with barcode and read ids')
    args = parser.parse_args(args)
    return args


def set_logger(logger_instance):
    logger_instance.setLevel(logging.INFO)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger_instance.addHandler(ch)


def run_pipeline(args):
    logger.info("Loading genome from " + args.reference)
    chr_records = SeqIO.to_dict(SeqIO.parse(args.reference, "fasta"))

    if not args.tsv:
        logger.info("Loading alignments from " + args.bam_pb)
        bam_records1 = load_bam({}, args.bam_pb)
        logger.info("Counting error rates...")
        stat_list = error_rate_single(bam_records1, chr_records)
    else:
        logger.info("Loading read pairs from " + args.tsv)
        read_pairs = load_tsv(args.tsv)
        logger.info("Loading alignments from " + args.bam_pb)
        bam_records1 = load_bam(set(map(lambda x: x[0], read_pairs)), args.bam_pb)
        logger.info("Loading alignments from " + args.bam_ont)
        bam_records2 = load_bam(set(map(lambda x: x[1], read_pairs)), args.bam_ont)
        logger.info("Counting error rates...")
        stat_list = error_rate_stats(read_pairs, bam_records1, bam_records2, chr_records)

    logger.info("Saving stats to " + args.output)
    for s in stat_list:
        stat_fname = args.output + s[0] + ".tsv"
        outf = open(stat_fname, "w")
        for stat_obj in s[1]:
            stat_obj.dump(outf)
        outf.close()

    logger.info("Done")


def main(args):
    args = parse_args(args)
    set_logger(logger)
    run_pipeline(args)


if __name__ == "__main__":
    # stuff only to run when not called via 'import' here
    try:
        main(sys.argv[1:])
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
