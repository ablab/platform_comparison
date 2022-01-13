#!/usr/bin/env python3
#
# ############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

from gtf2db import *
from src.dataset_processor import *

logger = logging.getLogger('CSA')

PACBIO_CCS_DATA = 'pacbio_ccs'
NANOPORE_DATA = 'nanopore'
DATATYPE = {PACBIO_CCS_DATA, NANOPORE_DATA}


def parse_args(args=None, namespace=None):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("--output", "-o", help="output folder, will be created automatically [default=output]",
                        type=str, default="output")

    parser.add_argument("--genedb", "-g", help="gene database in gffutils DB format (use gft2db.py to convert)", type=str,
                        required=True)
    parser.add_argument("--reference", "-r", help="reference genome in FASTA format, "
                                                  "should be provided to compute some additional stats and "
                                                  "when reads in FASTA/FASTQ are used as an input", type=str,
                        required=True)

    parser.add_argument('--bam', type=str, help='sorted and indexed BAM file')
    parser.add_argument("--data_type", "-d", type=str, required=True, choices=DATATYPE,
                        help="type of data to process, supported types are: " + ", ".join(DATATYPE))

    parser.add_argument("--delta", type=int, default=None,
                        help="delta for inexact splice junction comparison, chosen automatically based on data type")
    parser.add_argument("--matching_strategy", choices=["exact", "precise", "default", "loose"],
                        help="read-to-isoform matching strategy from the most strict to least", type=str, default=None)

    parser.add_argument("--threads", "-t", help="number of threads to use", type=int, default="16")
    parser.add_argument('--check_canonical', action='store_true', default=False,
                        help="report whether splice junctions are canonical (requires reference genome)")
    parser.add_argument("--sqanti_output", help="produce SQANTI-like TSV output (requires more time)",
                        action='store_true', default=False)

    parser.add_argument('--simple_intron_comparison', action='store_true', default=False,
                        help="use simple intron chain comparison, report simple stats")
    parser.add_argument('--intron_stats', action='store_true', default=False,
                        help="count intron stats, will force --simple_intron_comparison and --check_canonical")

    args = parser.parse_args(args, namespace)

    if os.path.exists(args.output):
        # logger is not defined yet
        print("WARNING! Output folder already exists, some files may be overwritten")
    else:
        os.makedirs(args.output)

    if not check_params(args):
        parser.print_usage()
        exit(-1)
    return args


# Check user's params
def check_params(args):
    args.input_data = InputDataStorage(args)
    if args.input_data.input_type == "fastq" and args.reference is None and args.index is None:
        print("ERROR! Reference genome or index were not provided, reads cannot be processed")
        return False
    if args.data_type is None:
        print("ERROR! Data type is not provided, choose one of " + " ".join(DATATYPE))
        return False
    elif args.data_type not in DATATYPE:
        print("ERROR! Unsupported data type " + args.data_type + ", choose one of: " + " ".join(DATATYPE))
        return False

    check_input_files(args)
    return True


def check_input_files(args):
    for sample in args.input_data.samples:
        for lib in sample.file_list:
            for in_file in lib:
                if not os.path.isfile(in_file):
                    print("ERROR! Input file " + in_file + " does not exist")
                    exit(-1)
                if args.input_data.input_type == "bam":
                    # TODO: sort and index file if needed
                    bamfile_in = pysam.AlignmentFile(in_file, "rb")
                    if not bamfile_in.has_index():
                        print("ERROR! BAM file " + in_file + " is not indexed, run samtools sort and samtools index")
                        exit(-1)
                    bamfile_in.close()

    if args.genedb is not None:
        if not os.path.isfile(args.genedb):
            print("ERROR! Gene database " + args.genedb + " does not exist")
            exit(-1)


def create_output_dirs(args):
    for sample in args.input_data.samples:
        sample_dir = sample.out_dir
        if os.path.exists(sample_dir):
            logger.warning(sample_dir + " folder already exists, some files may be overwritten")
        else:
            os.makedirs(sample_dir)
        sample_aux_dir = sample.aux_dir
        if os.path.exists(sample_aux_dir):
            logger.warning(sample_aux_dir + " folder already exists, some files may be overwritten")
        else:
            os.makedirs(sample_aux_dir)


def set_logger(args, logger_instance):
    logger_instance.setLevel(logging.INFO)
    log_file = os.path.join(args.output, "csa.log")
    f = open(log_file, "w")
    f.write("CMD: " + ' '.join(sys.argv) + '\n')
    f.close()
    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.INFO)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    logger_instance.addHandler(fh)
    logger_instance.addHandler(ch)


def set_data_dependent_options(args):
    matching_strategies = {PACBIO_CCS_DATA: "precise", NANOPORE_DATA: "default"}
    if args.matching_strategy is None:
        args.matching_strategy = matching_strategies[args.data_type]

    args.resolve_ambiguous = 'default'
    args.needs_polya_for_construction = False
    args.no_secondary = False


def set_matching_options(args):
    MatchingStrategy = namedtuple('MatchingStrategy',
                                  ('delta', 'max_intron_shift', 'max_missed_exon_len', 'max_fake_terminal_exon_len',
                                   'max_suspicious_intron_abs_len', 'max_suspicious_intron_rel_len',
                                   'resolve_ambiguous', 'correct_minor_errors'))

    strategies = {
        'exact':   MatchingStrategy(0, 0, 0, 0, 0, 0.0, 'monoexon_only', False),
        'precise': MatchingStrategy(4, 30, 50, 20, 0, 0.0, 'monoexon_and_fsm', True),
        'default': MatchingStrategy(6, 60, 100, 40, 60, 1.0, 'monoexon_and_fsm', True),
        'loose':   MatchingStrategy(12, 60, 100, 40, 60, 1.0, 'all',  True),
    }

    strategy = strategies[args.matching_strategy]

    args.delta = args.delta if args.delta is not None else strategy.delta
    args.minor_exon_extension = 50
    args.major_exon_extension = 300
    args.max_intron_shift = strategy.max_intron_shift
    args.max_missed_exon_len = strategy.max_missed_exon_len
    args.max_fake_terminal_exon_len = strategy.max_fake_terminal_exon_len
    # short introns that are actually long deletions, fix minimaps logic
    args.max_suspicious_intron_abs_len = strategy.max_suspicious_intron_abs_len
    args.max_suspicious_intron_rel_len = strategy.max_suspicious_intron_rel_len
    args.min_abs_exon_overlap = 10
    args.min_rel_exon_overlap = 0.2
    args.micro_intron_length = 50
    args.max_intron_abs_diff = min(30, args.max_intron_shift)
    args.max_intron_rel_diff = 0.2
    args.apa_delta = args.minor_exon_extension
    args.minimal_exon_overlap = 5
    args.minimal_intron_absence_overlap = 20
    args.polya_window = 16
    args.polya_fraction = 0.75
    if args.resolve_ambiguous == 'default':
        args.resolve_ambiguous = strategy.resolve_ambiguous
    if args.resolve_ambiguous not in AmbiguityResolvingMethod.__dict__:
        logger.error("Incorrect resolving ambiguity method: " + args.resolve_ambiguous + ", default will be used")
        args.resolve_ambiguous = strategy.resolve_ambiguous
    args.resolve_ambiguous = AmbiguityResolvingMethod[args.resolve_ambiguous]
    args.correct_minor_errors = strategy.correct_minor_errors

    updated_strategy = MatchingStrategy(args.delta, args.max_intron_shift, args.max_missed_exon_len,
                                        args.max_fake_terminal_exon_len,
                                        args.max_suspicious_intron_abs_len, args.max_suspicious_intron_rel_len,
                                        args.resolve_ambiguous, args.correct_minor_errors)
    logger.debug('Using %s strategy. Updated strategy: %s.' % (args.matching_strategy, updated_strategy))


def set_additional_params(args):
    set_data_dependent_options(args)
    set_matching_options(args)

    args.print_additional_info = True
    args.no_polya = False

    args.indel_near_splice_site_dist = 10
    args.upstream_region_len = 20

    if args.intron_stats:
        args.check_canonical = True
        args.simple_intron_comparison = True

    args.multimap_strategy = "take_best"
    multimap_strategies = {}
    for e in MultimapResolvingStrategy:
        multimap_strategies[e.name] = e.value
    args.multimap_strategy = MultimapResolvingStrategy(multimap_strategies[args.multimap_strategy])

    args.needs_reference = True
    if args.needs_reference and not args.reference:
        logger.warning("Reference genome is not provided! This may affect quality of the results!")
        args.needs_reference = False

    args.do_correction = False
    args.do_assignment = True


def run_pipeline(args):
    logger.info(" === Read correction started === ")
    dataset_processor = DatasetProcessor(args)
    dataset_processor.process_all_samples(args.input_data)
    logger.info(" === Read correction finished === ")


def main(args):
    args = parse_args(args)
    set_logger(args, logger)
    create_output_dirs(args)
    set_additional_params(args)
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
