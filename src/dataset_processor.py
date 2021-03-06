############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import _pickle as pickle
import gc
import glob
import gzip
import os
from multiprocessing import Pool
from collections import namedtuple
from concurrent import futures

from src.input_data_storage import *
from src.alignment_processor import *
from src.assignment_io import *
from src.multimap_resolver import *
from src.stats import *
from src.alignment_processor_simple import *


logger = logging.getLogger('CSA')


class GeneClusterConstructor:
    MAX_GENE_CLUSTER = 50
    MAX_GENE_LEN = 100000

    def __init__(self, gene_db):
        self.gene_db = gene_db
        self.gene_sets = None

    def get_gene_sets(self):
        if self.gene_sets is None:
            self.gene_sets = self.fill_gene_sets()
        return self.gene_sets

    def fill_gene_sets(self):
        gene_sets = []
        current_gene_db_list = []
        for g in self.gene_db.features_of_type('gene', order_by=('seqid', 'start')):
            gene_name = g.id
            gene_db = self.gene_db[gene_name]

            if len(current_gene_db_list) > 0 and \
                    (all(not genes_overlap(cg, gene_db) for cg in current_gene_db_list) or
                     (len(current_gene_db_list) > self.MAX_GENE_CLUSTER and
                      all(not genes_contain(cg, gene_db) for cg in current_gene_db_list))):
                gene_sets.append(current_gene_db_list)
                current_gene_db_list = []

            if g.end - g.start > self.MAX_GENE_LEN:
                gene_sets.append([gene_db])
            else:
                current_gene_db_list.append(gene_db)

        if current_gene_db_list:
            gene_sets.append(current_gene_db_list)
        return gene_sets



def assign_reads_in_parallel(sample, chr_id, cluster, args, current_chr_record):
    tmp_printer = TmpFileAssignmentPrinter("{}_{}".format(sample.out_raw_file, chr_id), args)
    intron_info_printer = IntronInfoPrinter(sample.intron_info_tsv + chr_id + ".tsv", args) if args.do_assignment else None
    processed_reads = []
    bam_files = list(map(lambda x: x[0], sample.file_list))
    bam_file_pairs = [(pysam.AlignmentFile(bam, "rb"), bam) for bam in bam_files]
    logger.info("Processing chromosome " + chr_id)
    gffutils_db = gffutils.FeatureDB(args.genedb, keep_order=True)
    for g in cluster:
        if len(g) > 100:
            logger.debug("Potential slowdown in %s due to large gene cluster of size %d" % (chr_id, len(g)))
        gene_info = GeneInfo(g, gffutils_db, args.delta)
        if args.intron_stats:
            alignment_processor = LongReadSimpleAlignmentProcessor(gene_info, bam_files, args,
                                                                   current_chr_record)
            assignment_storage = alignment_processor.process(intron_info_printer)
        else:
            alignment_processor = LongReadAlignmentProcessor(gene_info, bam_file_pairs, args,
                                                             current_chr_record)
            assignment_storage = alignment_processor.process()
        gene_info.db = None
        tmp_printer.add_gene_info(gene_info)
        for read_assignment in assignment_storage:
            tmp_printer.add_read_info(read_assignment)
            processed_reads.append(BasicReadAssignment(read_assignment, gene_info))
    logger.info("Finished processing chromosome " + chr_id)
    return processed_reads


class ReadAssignmentLoader:
    def __init__(self, save_file_name, gffutils_db, multimappers_disct):
        logger.info("Loading read assignments from " + save_file_name)
        assert os.path.exists(save_file_name)
        self.save_file_name = save_file_name
        self.unpickler = pickle.Unpickler(open(save_file_name, "rb"), fix_imports=False)
        self.current_gene_info_obj = None
        self.is_updated = False
        self.gffutils_db = gffutils_db
        self.multimapped_reads = multimappers_disct


def load_assigned_reads(save_file_name, gffutils_db, multimapped_reads):
    gc.disable()
    logger.info("Loading read assignments from " + save_file_name)
    assert os.path.exists(save_file_name)
    save_file_name = save_file_name
    unpickler = pickle.Unpickler(open(save_file_name, "rb"), fix_imports=False)
    read_storage = []
    current_gene_info = None

    while True:
        try:
            obj = unpickler.load()
            if isinstance(obj, ReadAssignment):
                read_assignment = obj
                assert current_gene_info is not None
                read_assignment.gene_info = current_gene_info
                if read_assignment.read_id in multimapped_reads:
                    resolved_assignment = None
                    for a in multimapped_reads[read_assignment.read_id]:
                        if a.start == read_assignment.start() and a.end == read_assignment.end() and \
                                a.gene_id == current_gene_info.gene_db_list[0].id and \
                                a.chr_id == read_assignment.chr_id:
                            if resolved_assignment is not None:
                                logger.warning("Duplicate read: %s %s %s" % (read_assignment.read_id, a.gene_id, a.chr_id))
                            resolved_assignment = a

                    if not resolved_assignment:
                        logger.warning("Incomplete information on read %s" % read_assignment.read_id)
                        continue
                    elif resolved_assignment.assignment_type == ReadAssignmentType.noninformative:
                        continue
                    else:
                        read_assignment.assignment_type = resolved_assignment.assignment_type
                        read_assignment.multimapper = resolved_assignment.multimapper
                read_storage.append(read_assignment)
            elif isinstance(obj, GeneInfo):
                if current_gene_info and read_storage:
                    yield current_gene_info, read_storage
                read_storage = []
                current_gene_info = obj
                current_gene_info.db = gffutils_db
            else:
                raise ValueError("Read assignment file {} is corrupted!".format(save_file_name))
        except EOFError:
            break
    gc.enable()
    if current_gene_info and read_storage:
        yield current_gene_info, read_storage


# Class for processing all samples against gene database
class DatasetProcessor:
    def __init__(self, args):
        self.args = args
        logger.info("Loading gene database from " + self.args.genedb)
        self.gffutils_db = gffutils.FeatureDB(self.args.genedb, keep_order=True)
        if self.args.needs_reference:
            logger.info("Loading reference genome from " + self.args.reference)
            _, outer_ext = os.path.splitext(self.args.reference)
            if outer_ext.lower() in ['.gz', '.gzip']:
                with gzip.open(self.args.reference, "rt") as handle:
                    self.reference_record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
            else:
                self.reference_record_dict = SeqIO.to_dict(SeqIO.parse(self.args.reference, "fasta"))
        else:
            self.reference_record_dict = None
        self.gene_cluster_constructor = GeneClusterConstructor(self.gffutils_db)
        self.gene_clusters = self.gene_cluster_constructor.get_gene_sets()
        self.io_support = IOSupport(self.args)
        self.multimapped_reads = defaultdict(list)

    def process_all_samples(self, input_data):
        logger.info("Processing " + proper_plural_form("sample", len(input_data.samples)))
        for sample in input_data.samples:
            self.process_sample(sample)
        logger.info("Processed " + proper_plural_form("sample", len(input_data.samples)))

    # Run though all genes in db and count stats according to alignments given in bamfile_name
    def process_sample(self, sample):
        logger.info("Processing sample " + sample.label)
        logger.info("Sample has " + proper_plural_form("BAM file", len(sample.file_list)) + ": " + ", ".join(map(lambda x: x[0], sample.file_list)))
        self.multimapped_reads = defaultdict(list)

        self.assign_reads(sample)
        saves_file = sample.out_raw_file
        logger.info('Read assignments files saved to {}*.'.
                    format(sample.out_raw_file))
        total_alignments, polya_found = self.load_read_info(saves_file)

        polya_fraction = polya_found / total_alignments if total_alignments > 0 else 0.0
        logger.info("Total alignments processed: %d, polyA tail detected in %d (%.1f%%)" %
                    (total_alignments, polya_found, polya_fraction * 100.0))
        self.args.needs_polya_for_construction = polya_fraction >= 0.7

        self.process_assigned_reads(sample, saves_file)
        for f in glob.glob(saves_file + "_*"):
            os.remove(f)
        logger.info("Processed sample " + sample.label)

    def assign_reads(self, sample):
        logger.info('Assigning reads to isoforms')
        chrom_clusters = self.get_chromosome_gene_clusters()
        pool = Pool(self.args.threads)
        processed_reads = pool.starmap(assign_reads_in_parallel, [(sample, chr_id, c, self.args,
                                                                   (self.reference_record_dict[
                                                                        chr_id] if self.reference_record_dict else None))
                                                                  for (chr_id, c) in chrom_clusters],
                                       chunksize=1)
        pool.close()
        pool.join()

        logger.info("Resolving multimappers")
        self.multimapped_reads = defaultdict(list)
        for storage in processed_reads:
            for read_assignment in storage:
                self.multimapped_reads[read_assignment.read_id].append(read_assignment)

        multimap_resolver = MultimapResolver(self.args.multimap_strategy)
        multimap_pickler = pickle.Pickler(open(sample.out_raw_file + "_multimappers", "wb"),  -1)
        multimap_pickler.fast = True
        total_assignments = 0
        polya_assignments = 0
        for read_id in list(self.multimapped_reads.keys()):
            assignment_list = self.multimapped_reads[read_id]
            if len(assignment_list) == 1:
                total_assignments += 1
                polya_assignments += 1 if assignment_list[0].polyA_found else 0
                del self.multimapped_reads[read_id]
                continue
            multimap_resolver.resolve(assignment_list)
            multimap_pickler.dump(assignment_list)
            for a in assignment_list:
                if a.assignment_type != ReadAssignmentType.noninformative:
                    total_assignments += 1
                    if a.polyA_found:
                        polya_assignments += 1

        info_pickler = pickle.Pickler(open(sample.out_raw_file + "_info", "wb"),  -1)
        info_pickler.dump(total_assignments)
        info_pickler.dump(polya_assignments)
        if total_assignments == 0:
            logger.warning("No reads were assigned to isoforms, check your input files")
        else:
            logger.info('Finishing read assignment, total assignments %d, polyA percentage %.1f' %
                        (total_assignments, 100 * polya_assignments / total_assignments))

    def process_assigned_reads(self, sample, dump_filename):
        logger.info("Processing assigned reads " + sample.label)
        self.create_aggregators(sample)

        for chr_id in self.get_chromosome_list():
            chr_dump_file = dump_filename + "_" + chr_id
            chr_record = self.reference_record_dict[chr_id] if self.reference_record_dict else None
            # future_list = []

            for gene_info, assignment_storage in load_assigned_reads(chr_dump_file, self.gffutils_db, self.multimapped_reads):
                for read_assignment in assignment_storage:
                    self.pass_to_aggregators(read_assignment)

        self.finalize_aggregators(sample)

    def load_read_info(self, dump_filename):
        gc.disable()
        if not self.multimapped_reads:
            multimap_unpickler = pickle.Unpickler(open(dump_filename + "_multimappers", "rb"), fix_imports=False)
            while True:
                try:
                    obj = multimap_unpickler.load()
                    if isinstance(obj, list):
                        assignment_list = obj
                        read_id = assignment_list[0].read_id
                        self.multimapped_reads[read_id] = assignment_list
                    else:
                        raise ValueError("Multimap assignment file {} is corrupted!".format(dump_filename))
                except EOFError:
                    break

        info_unpickler = pickle.Unpickler(open(dump_filename + "_info", "rb"), fix_imports=False)
        total_assignments = info_unpickler.load()
        polya_assignments = info_unpickler.load()

        gc.enable()
        return total_assignments, polya_assignments

    def get_chromosome_list(self):
        chromosomes = []
        current_chromosome = ""
        for g in self.gene_clusters:
            chr_id = g[0].seqid
            if chr_id != current_chromosome:
                chromosomes.append(chr_id)
                current_chromosome = chr_id
        return chromosomes

    def get_chromosome_gene_clusters(self):
        chrom_clusters = []
        cur_cluster = []
        current_chromosome = ""
        for g in self.gene_clusters:
            chr_id = g[0].seqid
            if chr_id != current_chromosome:
                if cur_cluster:
                    chrom_clusters.append((current_chromosome, cur_cluster))
                    cur_cluster = []
                current_chromosome = chr_id
            cur_cluster.append(g)
        if cur_cluster:
            chrom_clusters.append((current_chromosome, cur_cluster))

        # chromosomes with large clusters take more time and should go first
        chrom_clusters = sorted(chrom_clusters, key=lambda x: sum(len(cluster) ** 2 for cluster in x[1]), reverse=True)
        return chrom_clusters

    def create_aggregators(self, sample):
        self.read_stat_counter = EnumStats()
        printer_list = []
        self.basic_printer = None
        if self.args.do_assignment:
            self.basic_printer = BasicTSVAssignmentPrinter(sample.out_assigned_tsv, self.args, self.io_support)
            printer_list.append(self.basic_printer)
        if self.args.do_correction:
            self.corrected_bed_printer = BEDPrinter(sample.out_corrected_bed, self.args, print_corrected=True)
            printer_list.append(self.corrected_bed_printer)
        if self.args.sqanti_output:
            self.sqanti_printer = SqantiTSVPrinter(sample.out_alt_tsv, self.args, self.io_support)
            printer_list.append(self.sqanti_printer)
        self.global_printer = ReadAssignmentCompositePrinter(printer_list)


    def pass_to_aggregators(self, read_assignment):
        if read_assignment is None:
            return
        self.read_stat_counter.add(read_assignment.assignment_type)
        self.global_printer.add_read_info(read_assignment)

    def finalize_aggregators(self, sample):
        if self.basic_printer:
            logger.info("Read assignments are stored in " + self.basic_printer.output_file_name)
        self.read_stat_counter.print_start("Read assignment statistics")

