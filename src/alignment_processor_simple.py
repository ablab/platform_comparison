############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
import pysam
from Bio import SeqIO

from src.long_read_assigner import *
from src.long_read_simple_assigner import *

from src.long_read_profiles import *
from src.polya_finder import *
from src.polya_verification import *

logger = logging.getLogger('CSA')


class LongReadSimpleAlignmentProcessor:
    """ class for aggregating all assignment information

    Parameters
    ----------
    gene_info
    bams
    params
    printer
    counter
    """

    def __init__(self, gene_info, bams, params, chr_record=None):
        self.gene_info = gene_info
        self.bams = bams
        self.params = params
        self.chr_record = chr_record

        self.assigner = LongReadSimpleAssigner(self.gene_info, self.params)
        self.profile_constructor = CombinedProfileConstructor(gene_info, params)
        self.polya_finder = PolyAFinder(self.params.polya_window, self.params.polya_fraction)
        self.polya_fixer = PolyAFixer(self.params)
        self.assignment_storage = []
        self.gene_region = (gene_info.start, gene_info.end)

    def process(self, intron_printer):
        self.assignment_storage = []
        self.gene_info.all_read_region_start = self.gene_info.start
        self.gene_info.all_read_region_end = self.gene_info.end

        for b in self.bams:
            self.process_single_file(b, intron_printer)

        if self.params.needs_reference:
            self.gene_info.all_read_region_start -= self.params.upstream_region_len
            self.gene_info.all_read_region_end += self.params.upstream_region_len
            self.gene_info.reference_region = \
                str(self.chr_record[self.gene_info.all_read_region_start - 1:self.gene_info.all_read_region_end + 1].seq)
            self.gene_info.canonical_sites = {}
        return self.assignment_storage

    def process_single_file(self, bamfile_in, intron_printer):
        processed_reads = set()
        for genic_region in self.gene_info.genic_regions:
            for alignment in bamfile_in.fetch(self.gene_info.chr_id, genic_region[0], genic_region[1]):
                read_id = alignment.query_name
                if alignment.reference_id == -1:
                    self.assignment_storage.append(ReadAssignment(read_id, None))
                    continue
                if alignment.is_supplementary:
                    continue
                if self.params.no_secondary and alignment.is_secondary:
                    continue

                logger.debug("=== Processing read " + read_id + " ===")

                # concat indels
                concat_blocks = concat_gapless_blocks(sorted(alignment.get_blocks()), alignment.cigartuples)
                if not concat_blocks:
                    logger.warning("Read %s has no aligned exons" % read_id)
                    continue
                # correct coordinates to GTF style (inclusive intervals)
                sorted_blocks = correct_bam_coords(concat_blocks)
                read_start = sorted_blocks[0][0]
                read_end = sorted_blocks[-1][1]
                read_tuple = (read_id, read_start, read_end)
                if read_tuple in processed_reads:
                    continue
                processed_reads.add(read_tuple)
                logger.debug("Read exons: " + str(sorted_blocks))
                if self.params.needs_reference:
                    if read_start < self.gene_info.all_read_region_start:
                        self.gene_info.all_read_region_start = read_start
                    if read_end > self.gene_info.all_read_region_end:
                        self.gene_info.all_read_region_end = read_end

                #polya_info = self.polya_finder.detect_polya(alignment)
                #sorted_blocks, polya_info, exon_changed = self.polya_fixer.correct_read_info(sorted_blocks, polya_info)
                polya_info = PolyAInfo(-1, -1, -1, -1)

                combined_profile = self.profile_constructor.construct_profiles(sorted_blocks, polya_info)
                read_assignment = self.assigner.assign_to_isoform(read_id, combined_profile)

                #if exon_changed:
                #    read_assignment.add_match_attribute(MatchEvent(MatchEventSubtype.aligned_polya_tail))
                read_assignment.polyA_found = (polya_info.external_polya_pos != -1 or
                                               polya_info.external_polyt_pos != -1 or
                                               polya_info.internal_polya_pos != -1 or
                                               polya_info.internal_polyt_pos != -1)
                read_assignment.polya_info = polya_info
                read_assignment.exons = sorted_blocks
                read_assignment.mapped_strand = "-" if alignment.is_reverse else "+"
                read_assignment.chr_id = self.gene_info.chr_id
                read_assignment.multimapper = alignment.is_secondary

                if intron_printer and self.chr_record and not alignment.is_secondary:
                    chr_id = self.gene_info.chr_id
                    read_start = sorted_blocks[0][0] - 10
                    read_end = sorted_blocks[-1][1] + 10
                    ref_region = str(self.chr_record[read_start - 1:read_end + 1].seq)
                    gene_profile_index = 0
                    gene_profile = combined_profile.read_intron_profile.gene_profile

                    for i, intron in enumerate(combined_profile.read_intron_profile.read_features):
                        if combined_profile.read_intron_profile.read_profile[i] != 1:
                            intron_printer.add_intron_info(read_id, chr_id, ".", intron, (0, 0), "novel", 0, 0, 0, 0, 0,
                                                           0)
                            continue

                        while gene_profile_index < len(gene_profile) and (
                                gene_profile[gene_profile_index] != 1 or not overlaps(intron,
                                                                                      self.gene_info.intron_profiles.features[
                                                                                          gene_profile_index])):
                            gene_profile_index += 1
                        if gene_profile_index >= len(gene_profile) or gene_profile[gene_profile_index] != 1:
                            logger.info("profile %s, index %d" % (str(gene_profile), gene_profile_index))
                            logger.info("read profile %s, index %d" % (
                            str(combined_profile.read_intron_profile.read_profile), i))

                        reference_intron = self.gene_info.intron_profiles.features[gene_profile_index]
                        gene_profile_index += 1
                        is_consistent = reference_intron == intron
                        left_diff = intron[0] - reference_intron[0]
                        right_diff = intron[1] - reference_intron[1]

                        ref_strand = self.check_canonical(reference_intron, ref_region, read_start)
                        read_strand = self.check_canonical(intron, ref_region, read_start)
                        if ref_strand is None:
                            # reference intron is non-canonical
                            if read_strand is None:
                                # read intron is non-canonical as well
                                intron_printer.add_intron_info(read_id, chr_id, ".", intron, reference_intron,
                                                               "both_noncanonical", 0, 0, 0, 0, left_diff, right_diff)
                            else:
                                # read is canonical
                                strand, donor_up, donor_down, acceptor_up, acceptor_down = \
                                    self.analyse_intron_sites(intron, ref_region, read_start, read_strand)
                                if read_strand == '+':
                                    donor_diff = left_diff
                                    acc_diff = right_diff
                                else:
                                    donor_diff = - right_diff
                                    acc_diff = - left_diff
                                intron_printer.add_intron_info(read_id, chr_id, read_strand, intron, reference_intron,
                                                               "reference_noncanonical",
                                                               donor_up, donor_down, acceptor_up,
                                                               acceptor_down, donor_diff, acc_diff)
                        else:
                            strand, donor_up, donor_down, acceptor_up, acceptor_down = \
                                self.analyse_intron_sites(intron, ref_region, read_start, ref_strand)
                            if ref_strand == '+':
                                donor_diff = left_diff
                                acc_diff = right_diff
                            else:
                                donor_diff = - right_diff
                                acc_diff = - left_diff
                            if ref_strand != read_strand:
                                intron_printer.add_intron_info(read_id, chr_id, ref_strand, intron, reference_intron,
                                                               "read_noncanonical",
                                                               donor_up, donor_down, acceptor_up,
                                                               acceptor_down, donor_diff, acc_diff)
                            else:
                                intron_type = "consistent" if is_consistent else "incosistent"
                                intron_printer.add_intron_info(read_id, chr_id, strand, intron, reference_intron,
                                                               intron_type,
                                                               donor_up, donor_down, acceptor_up, acceptor_down,
                                                               donor_diff, acc_diff)

                if self.params.sqanti_output:
                    indel_count, junctions_with_indels = self.count_indel_stats(alignment)
                    read_assignment.set_additional_info("indel_count", indel_count)
                    read_assignment.set_additional_info("junctions_with_indels", junctions_with_indels)
                    read_assignment.introns_match = \
                        all(e == 1 for e in combined_profile.read_intron_profile.read_profile)

                self.assignment_storage.append(read_assignment)
                logger.debug("=== Finished read " + read_id + " ===")

    def check_canonical(self, intron, ref_region, read_start, strand=None):
        intron_left_pos = intron[0] - read_start
        intron_right_pos = intron[1] - read_start
        left_site = ref_region[intron_left_pos:intron_left_pos + 2]
        right_site = ref_region[intron_right_pos - 1:intron_right_pos + 1]
        if (left_site == "GT" and right_site == "AG") and strand != '-':
            return '+'
        elif (left_site == "CT" and right_site == "AC") and strand != '+':
            return '-'
        else:
            return None

    def analyse_intron_sites(self, intron, ref_region, read_start, strand=None):
        seq_size = 10
        intron_left_pos = intron[0] - read_start
        intron_right_pos = intron[1] - read_start

        if strand is None:
            left_site = ref_region[intron_left_pos:intron_left_pos + 2]
            right_site = ref_region[intron_right_pos - 1:intron_right_pos + 1]
            if left_site == "GT" and right_site == "AG":
                strand = '+'
            elif left_site == "CT" and right_site == "AC":
                strand = '-'
            else:
                return None, None, None, None, None

        if strand not in ['+', '-']:
            return None, None, None, None, None

        left_upper = ref_region[intron_left_pos - seq_size:intron_left_pos]
        left_lower = ref_region[intron_left_pos + 2:intron_left_pos + seq_size + 2]
        right_upper = ref_region[intron_right_pos - seq_size - 1:intron_right_pos - 1]
        right_lower = ref_region[intron_right_pos + 1:intron_right_pos + seq_size + 1]

        # upstream and downstream here are relative to the genome
        if strand == "+":
            donor_upstream = left_upper.rfind("GT")
            donor_downstream = left_lower.find("GT")
            acc_upstream = right_upper.rfind("AG")
            acc_downstream = right_lower.find("AG")
        else:
            acc_upstream = left_upper.rfind("CT")
            acc_downstream = left_lower.find("CT")
            donor_upstream = right_upper.rfind("AC")
            donor_downstream = right_lower.find("AC")

        donor_upstream = seq_size - donor_upstream if donor_upstream != -1 else 0
        donor_downstream = 2 + donor_downstream if donor_downstream != -1 else 0
        acc_upstream = seq_size - acc_upstream if acc_upstream != -1 else 0
        acc_downstream = 2 + acc_downstream if acc_downstream != -1 else 0

        if strand == '+':
            return strand, donor_upstream, donor_downstream, acc_upstream, acc_downstream
        else:
            return strand, donor_downstream, donor_upstream, acc_downstream, acc_upstream

    def count_indel_stats(self, alignment):
        cigar_event_count = len(alignment.cigartuples)
        indel_events = [1, 2]
        indel_count = 0
        intron_cigar_positions = []
        for i in range(cigar_event_count):
            cigar = alignment.cigartuples[i]
            if cigar[0] in indel_events:
                indel_count += 1
            elif cigar[0] == 3:
                intron_cigar_positions.append(i)

        junctions_with_indels = 0
        for i in intron_cigar_positions:
            # indel right near intron
            if (i > 0 and alignment.cigartuples[i - 1][0] in indel_events) or \
                    (i < cigar_event_count - 1 and alignment.cigartuples[i + 1][0] in indel_events):
                junctions_with_indels += 1

            # indel separated by at most 'indel_near_splice_site_dist' matches from intron
            if (i > 1 and alignment.cigartuples[i - 2][0] in indel_events and
                alignment.cigartuples[i - 1][0] in [0, 7, 8] and
                alignment.cigartuples[i - 1][1] <= self.params.indel_near_splice_site_dist) or \
                    (i < cigar_event_count - 2 and alignment.cigartuples[i + 2][0] in indel_events and
                     alignment.cigartuples[i + 1][0]  in [0, 7, 8] and
                     alignment.cigartuples[i + 1][1] <= self.params.indel_near_splice_site_dist):
                junctions_with_indels += 1

        return indel_count, junctions_with_indels
