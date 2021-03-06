############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from functools import partial
from collections import defaultdict

from src.common import *

logger = logging.getLogger('CSA')


class MappedReadProfile:
    def __init__(self, gene_profile, read_profile, read_features):
        self.gene_profile = gene_profile
        self.read_profile = read_profile
        self.read_features = read_features


class CombinedReadProfiles:
    def __init__(self, read_intron_profile, read_exon_profile, read_split_exon_profile,
                 polya_info=None, cage_hits=-1, alignment=None):
        self.read_intron_profile = read_intron_profile
        self.read_exon_profile = read_exon_profile
        self.read_split_exon_profile = read_split_exon_profile
        self.alignment = alignment
        self.polya_info = polya_info
        self.cage_hits = cage_hits


# The following 2 classes are very similar, but lets keep them separately for now
# accepts sorted gapless alignment blocks
class OverlappingFeaturesProfileConstructor:
    # ignore_terminal -- bool flag, indicates whether to ignore leading and trailing -1s in the profile
    def __init__(self, known_features, gene_region,
                 comparator = partial(equal_ranges, delta=0),
                 absence_condition = contains,
                 delta=0):
        self.known_features = known_features
        self.gene_region = gene_region
        self.comparator = comparator
        self.absence_condition = absence_condition
        self.delta = delta

    def construct_intron_profile(self, sorted_blocks, polya_position=-1, polyt_position=-1):
        mapped_region = (sorted_blocks[0][0], sorted_blocks[-1][1])
        read_introns = junctions_from_blocks(sorted_blocks)
        return self.construct_profile_for_features(read_introns, mapped_region, polya_position, polyt_position)

    def construct_exon_profile(self, sorted_blocks, polya_position=-1, polyt_position=-1):
        mapped_region = (sorted_blocks[0][1] + self.delta, sorted_blocks[-1][0] - self.delta)
        return self.construct_profile_for_features(sorted_blocks, mapped_region, polya_position, polyt_position)

    def match_delta(self, feature1, feature2):
        return abs(feature1[0] - feature2[0]) + abs(feature1[1] - feature2[1])

    def match_genomic_features(self, read_features):
        matched_features = defaultdict(list)

        # TODO: starting value can be detected using binary search for long profiles
        gene_pos = 0
        read_pos = 0
        while gene_pos < len(self.known_features) and read_pos < len(read_features):
            if self.comparator(read_features[read_pos], self.known_features[gene_pos]):
                matched_features[read_pos].append(gene_pos)
                gene_pos += 1
            elif overlaps(read_features[read_pos], self.known_features[gene_pos]):
                gene_pos += 1
            elif left_of(read_features[read_pos], self.known_features[gene_pos]):
                read_pos += 1
            else:
                gene_pos += 1
        # eliminating non unique features
        for read_pos in matched_features.keys():
            if len(matched_features[read_pos]) > 1:
                deltas = [self.match_delta(read_features[read_pos], self.known_features[gene_pos])
                          for gene_pos in matched_features[read_pos]]
                best_match = min(deltas)
                filtered_features = []
                for i, d in enumerate(deltas):
                    if d == best_match:
                        filtered_features.append(matched_features[read_pos][i])
                matched_features[read_pos] = filtered_features

        corrected_features = []
        for i, intron in enumerate(read_features):
            if i in matched_features:
                # known feature
                # corresponding known feature, always take first for now
                corrected_features.append(self.known_features[matched_features[i][0]])
            else:
                corrected_features.append(intron)
        return corrected_features

    def construct_profile_for_features(self, read_features, mapped_region=(0, 0), polya_position=-1, polyt_position=-1):
        read_profile = [0] * (len(read_features))
        intron_profile = [0] * (len(self.known_features))
        matched_features = defaultdict(list)

        for i in range(len(intron_profile)):
            if self.absence_condition(mapped_region, self.known_features[i]):
                intron_profile[i] = -1
        for i in range(len(read_profile)):
            if self.absence_condition(self.gene_region, read_features[i]):
                read_profile[i] = -1

        # TODO: starting value can be detected using binary search for long profiles
        gene_pos = 0
        read_pos = 0
        while gene_pos < len(self.known_features) and read_pos < len(read_features):
            if self.comparator(read_features[read_pos], self.known_features[gene_pos]):
                intron_profile[gene_pos] = 1
                read_profile[read_pos] = 1
                matched_features[read_pos].append(gene_pos)
                gene_pos += 1
            elif overlaps(read_features[read_pos], self.known_features[gene_pos]):
                intron_profile[gene_pos] = -1
                gene_pos += 1
            elif left_of(read_features[read_pos], self.known_features[gene_pos]):
                if read_profile[read_pos] == 0 and gene_pos > 0:
                    read_profile[read_pos] = -1
                read_pos += 1
            else:
                if read_pos > 0:
                    intron_profile[gene_pos] = -1
                gene_pos += 1

        # eliminating non unique features
        for read_pos in matched_features.keys():
            if len(matched_features[read_pos]) > 1:
                deltas = [self.match_delta(read_features[read_pos], self.known_features[gene_pos])
                          for gene_pos in matched_features[read_pos]]
                best_match = min(deltas)
                for i in range(len(matched_features[read_pos])):
                    if deltas[i] > best_match:
                        intron_profile[matched_features[read_pos][i]] = -1

        corrected_features = []
        for i, v in enumerate(read_profile):
            if v == 1:
                # known feature
                # corresponding known feature, always take first for now
                corrected_features.append(self.known_features[matched_features[i][0]])
            else:
                corrected_features.append(read_features[i])

        # making everying beyond polyA tail as outside feature
        if polya_position != -1:
            for i in range(len(self.known_features)):
                # feature is surely beyond polyA tail
                if self.known_features[i][0] > polya_position + self.delta:
                    intron_profile[i] = -2
        if polyt_position != -1:
            for i in range(len(self.known_features)):
                # feature is surely before polyT tail
                if self.known_features[i][1] < polyt_position - self.delta:
                    intron_profile[i] = -2

        return MappedReadProfile(intron_profile, read_profile, read_features)


# accepts sorted gapless alignment blocks
class NonOverlappingFeaturesProfileConstructor:
    def __init__(self, known_exons, comparator=overlaps, delta=0):
        self.known_exons = known_exons
        self.comparator = comparator
        self.delta = delta

    def construct_profile(self, sorted_blocks, polya_position=-1, polyt_position=-1):
        exon_profile = [0] * (len(self.known_exons))
        read_profile = [0] * (len(sorted_blocks))
        read_exons = sorted_blocks
        gene_pos = 0
        read_pos = 0

        while gene_pos < len(self.known_exons) and read_pos < len(read_exons):
            if self.comparator(read_exons[read_pos], self.known_exons[gene_pos]):
                exon_profile[gene_pos] = 1
                read_profile[read_pos] = 1
                if read_exons[read_pos][1] < self.known_exons[gene_pos][1]:
                    read_pos += 1
                else:
                    gene_pos += 1
            elif overlaps(read_exons[read_pos], self.known_exons[gene_pos]):
                if read_exons[read_pos][1] < self.known_exons[gene_pos][1]:
                    read_pos += 1
                else:
                    gene_pos += 1
            elif left_of(read_exons[read_pos], self.known_exons[gene_pos]):
                if gene_pos > 0 and read_profile[read_pos] == 0:
                    read_profile[read_pos] = -1
                read_pos += 1
            else:
                if read_pos > 0 and exon_profile[gene_pos] == 0:
                    exon_profile[gene_pos] = -1
                gene_pos += 1

        # making everying beyond polyA tail as outside feature
        if polya_position != -1:
            for i in range(len(self.known_exons)):
                # feature is surely beyond polyA tail
                if self.known_exons[i][0] > polya_position + self.delta:
                    exon_profile[i] = -2
        if polyt_position != -1:
            for i in range(len(self.known_exons)):
                # feature is surely before polyT tail
                if self.known_exons[i][1] < polyt_position - self.delta:
                    exon_profile[i] = -2

        return MappedReadProfile(exon_profile, read_profile, read_exons)


class CombinedProfileConstructor:
    def __init__(self, gene_info, params):
        self.gene_info = gene_info
        self.params = params

        gene_region = (gene_info.start, gene_info.end)
        self.intron_profile_constructor = \
            OverlappingFeaturesProfileConstructor(self.gene_info.intron_profiles.features, gene_region,
                                                  comparator=partial(equal_ranges, delta=self.params.delta),
                                                  absence_condition=partial(overlaps_at_least, delta=self.params.minimal_intron_absence_overlap),
                                                  delta=self.params.delta)
        self.exon_profile_constructor = \
            OverlappingFeaturesProfileConstructor(self.gene_info.exon_profiles.features, gene_region,
                                                  comparator=partial(equal_ranges, delta=self.params.delta),
                                                  delta=self.params.delta)
        self.split_exon_profile_constructor = \
            NonOverlappingFeaturesProfileConstructor(self.gene_info.split_exon_profiles.features,
                                                     comparator=partial(overlaps_at_least,
                                                                        delta=self.params.minimal_exon_overlap),
                                                     delta=self.params.delta)

    def construct_profiles(self, sorted_blocks, polya_info):
        intron_profile = self.intron_profile_constructor.construct_intron_profile(sorted_blocks,
                                                                                  polya_info.external_polya_pos,
                                                                                  polya_info.external_polyt_pos)
        exon_profile = None
        split_exon_profile = self.split_exon_profile_constructor.construct_profile(sorted_blocks,
                                                                                   polya_info.external_polya_pos,
                                                                                   polya_info.external_polyt_pos)
        return CombinedReadProfiles(intron_profile, exon_profile, split_exon_profile,
                                    polya_info=polya_info)