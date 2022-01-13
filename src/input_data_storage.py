############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
import logging

logger = logging.getLogger('CSA')


class SampleData:
    def __init__(self, file_list, label, out_dir):
        # list of lists, since each sample may contain several libraries, and each library may contain 2 files (paired)
        self.file_list = file_list
        self.label = label
        self.out_dir = out_dir
        self.aux_dir = os.path.join(self.out_dir, "aux")
        self._init_paths()

    def _make_path(self, name):
        return os.path.join(self.out_dir, name)

    def _make_aux_path(self, name):
        return os.path.join(self.aux_dir, name)

    def _init_paths(self):
        self.out_assigned_tsv = self._make_path(self.label + ".read_assignments.tsv")
        self.out_raw_file = self._make_aux_path(self.label + ".save")
        # self.out_mapped_bed = self._make_path(self.label + ".mapped_reads.bed")
        self.out_corrected_bed = self._make_path(self.label + ".corrected_reads.bed")
        self.out_alt_tsv = self._make_path(self.label + ".SQANTI-like.tsv")
        self.intron_info_tsv = self._make_path(self.label + ".intron_info.")


class InputDataStorage:
    def __init__(self, args):
        # list of SampleData
        self.samples = []
        self.input_type = ""
        self.readable_names_dict = {}
        sample_files = []

        self.input_type = "bam"
        sample_files.append([[args.bam]])
        labels = self.get_labels(sample_files)

        for i in range(len(sample_files)):
            self.samples.append(SampleData(sample_files[i], labels[i], os.path.join(args.output, labels[i])))

    def get_labels(self, sample_files):
        labels = []
        for i in range(len(sample_files)):
            labels.append('{:02d}'.format(i) + "_" + os.path.splitext(os.path.basename(sample_files[i][0][0]))[0])
        return labels




