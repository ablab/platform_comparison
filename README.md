# Scripts for PacBio vs ONT platform comparison paper

Scripts used for data analysis are in `scripts/` folder:
- Aligned lengths analysis: `aligned_lengths.py` 
- Aligned lengths analysis (Smith-Waterman alignment): `sw_alignment.py` 
- Errors analysis:  `compare_error_rate.py` 
- K-mer identity analysis: `kmer_stats.py` 
- TSS/Poly-A sites analysis: `find_tss_polya.py` 
- Intron chains comparison: `compare_introns.py` 
- Alternative splicing analysis: `alt_splicing.py` 
- Analysis of introns w/respect to read quality: `intron_chain_q.py` and `intron_splice_quality.py`

Scripts located in the root directory:
- Assign reads to isoforms: `assign_reads_to_isoforms.py`
- Correct spliced alignments: `correct_spliced_alignments.py`
- Convert GTF to gffutils database: `gtf2db.py`
