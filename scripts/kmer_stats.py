import argparse

import pysam
from itertools import groupby
import gffutils
import sys
from Bio import Seq, SeqIO

from src.common import correct_bam_coords, concat_gapless_blocks

is_hp = 0 ## use homopolymer compression
k = 14


def parse_args(args=None, namespace=None):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("--reference", "-r", help="reference genome in FASTA format", type=str)
    parser.add_argument('--bam_pb', type=str, help='sorted and indexed BAM file for PacBio')
    parser.add_argument('--bam_ont', type=str, help='sorted and indexed BAM file for ONT')
    parser.add_argument('--tsv', type=str, help='TSV with barcode and read ids')

    parser.add_argument("--genedb", "-g", help="gene database in gffutils DB format or GTF/GFF format", type=str)

    # parser.add_argument('--assign_pb', type=str, help='read assignments file for PacBio')
    # parser.add_argument('--assign_ont', type=str, help='read assignments file for ONT')

    # parser.add_argument('--fastq_pb', type=str, help='FASTQ file for PacBio')
    # parser.add_argument('--fastq_ont', type=str, help='FASTQ file for ONT')
    args = parser.parse_args(args)
    return args


def rev_comp(seq):
    c = dict(zip('ATCGNatcgn', 'TAGCNtagcn'))
    return ''.join(c.get(nucleotide, '') for nucleotide in reversed(seq))


def compress_hp(seq):
    if not is_hp: return seq
    return ''.join(x[0] for x in groupby(list(seq)))


def load_bam(read_set, bamfile):
    bam_records = {}
    for r in pysam.AlignmentFile(bamfile, "rb").fetch():
        if not r.is_supplementary and not r.is_secondary and (not read_set or r.query_name in read_set):
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
        sequence_to_check = alignment.seq[
                            :clipped_len]  # str(Seq.Seq(alignment.seq[:clipped_len]).reverse_complement()).upper()

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
    cigar_clip = cigar.copy()
    if cigar[0][0] == 4:
        start_clip = cigar[0][1]
        cigar_clip = cigar_clip[1:]
    end_clip = 0
    if cigar[-1][0] == 4:
        end_clip = cigar[-1][1]
        cigar_clip = cigar_clip[:-1]
    return seq[start_clip:len(seq)-end_clip], cigar_clip


def canon_kmer(kmer):
    kmer = ''.join(x[0] for x in groupby(list(kmer)))
    return min(kmer, rev_comp(kmer))

args = parse_args(sys.argv[1:])
#db = gffutils.FeatureDB(args.genedb, keep_order = True)
chr_records = SeqIO.to_dict(SeqIO.parse(args.reference, "fasta"))

#pb_seqs = SeqIO.parse(args.fastq_pb, format="fastq")
#ont_seqs = SeqIO.parse(args.fastq_ont, format="fastq")

#pb_seq_dict = {seq.id: str(seq.seq) for seq in pb_seqs}
#ont_seq_dict = {seq.id: str(seq.seq) for seq in ont_seqs}

exon_lens = [21, 50, 100, 200, 20000]
exon_titles = ['<' + str(exon_lens[0])] + ['%d-%d' % (exon_lens[i]+1, exon_lens[i+1]) for i in range(len(exon_lens)-2)] + ['>' + str(exon_lens[-2]+1)]

ont_filt_data = []
ont_data = []
pb_data = []

def extract_exons(bam_record):
    #ref_region_start = bam_record1.reference_start
    #ref_region_end = ref_region_start + bam_record1.reference_length
    # direction1 = get_direction(bam_record1)
    cigar1 = bam_record.cigartuples

    read_seq, cigar1 = clip_seq(bam_record.seq, cigar1)
    read_exons = []
    cigar_index = 0
    read_index1, read_index2 = 0, 0
    while cigar_index < len(cigar1):
        # init new block
        cigar_event, cigar_len = cigar1[cigar_index]
        if cigar_event == 0:
            read_index2 += cigar_len
            cigar_event += 1
        elif cigar_event == 1:
            read_index2 += cigar_len
            cigar_event += 1
        elif cigar_event == 2:
            cigar_event += 1
        elif cigar_event == 3:
            read_exons.append((read_index1, read_index2))
            read_index1 = read_index2
        cigar_index+=1
    if read_index1 != read_index2:
        read_exons.append((read_index1, read_index2))
    ref_exons = correct_bam_coords(concat_gapless_blocks(sorted(bam_record1.get_blocks()), bam_record1.cigartuples))
    exon_scores = []
    for read_exon, ref_exon in zip(read_exons[1:-1], ref_exons[1:-1]):
        read_exon_seq = compress_hp(read_seq[read_exon[0]:read_exon[1]])
        read_kmers = set([str(read_exon_seq[i:i+k]) for i in range(len(read_exon_seq)-k+1)])
        read_kmers_rev = set([rev_comp(read_exon_seq[i:i+k]) for i in range(len(read_exon_seq)-k+1)])
        exon_seq = str(chr_records[bam_record1.reference_name][ref_exon[0]:ref_exon[1]+1].seq)
        exon_len = len(exon_seq)
        exon_seq = compress_hp(exon_seq)
        exon_kmers = set([str(exon_seq[i:i+k]) for i in range(len(exon_seq)-k+1)])
        if not exon_kmers: continue
        intsect = len(exon_kmers.intersection(read_kmers))*1.0/len(exon_kmers)
        intsect_rev = len(exon_kmers.intersection(read_kmers_rev))*1.0/len(exon_kmers)
        exon_scores.append((exon_len,max(intsect, intsect_rev)))
    return exon_scores

read_pairs = load_tsv(args.tsv)
bam_records1 = load_bam(set(map(lambda x: x[0], read_pairs)), args.bam_pb)
bam_records2 = load_bam(set(map(lambda x: x[1], read_pairs)), args.bam_ont)
for read_pair in read_pairs:
    if read_pair[0] not in bam_records1 or read_pair[1] not in bam_records2:
        continue
    bam_record1 = bam_records1[read_pair[0]]
    bam_record2 = bam_records2[read_pair[1]]
    if bam_record1.reference_name != bam_record2.reference_name:
        continue

    pb_exon_scores = extract_exons(bam_record1)
    if pb_exon_scores: pb_data.extend(pb_exon_scores)
    ont_exon_scores = extract_exons(bam_record2)
    if ont_exon_scores: ont_data.extend(ont_exon_scores)

pb_filt_data = []
pb_exons = []
for exon_len, kmer_pct in pb_data:
    for i, l in enumerate(exon_lens):
        if exon_len <= l:
            pb_filt_data.append(kmer_pct*100)
            pb_exons.append(exon_titles[i])
            break

ont_exons = []
for exon_len, kmer_pct in ont_data:
    for i, l in enumerate(exon_lens):
        if exon_len <= l: 
            ont_filt_data.append(kmer_pct*100)
            ont_exons.append(exon_titles[i])
            break

with open("kmer_stats_compress%d_k%d.txt" % (is_hp, k), "w") as f:
    f.write(" ".join([str(s) for s in (pb_filt_data+ont_filt_data)]))
    f.write("\n")
    f.write(" ".join(['PacBio'] * len(pb_filt_data) + ['ONT'] * len(ont_filt_data)))
    f.write("\n")
    f.write(" ".join([str(s) for s in (pb_exons+ont_exons)]))

