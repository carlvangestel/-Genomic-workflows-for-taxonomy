##TEST##



#8 Inspect MAPQ distribution before deciding threshold
#  view MAPQ stats: samtools view input.bam | awk '{print $5}' | sort -n | uniq -c
#  prints a histogram of MAPQ values to assist in picking a sensible cutoff (e.g., if 90% are ≥30, that’s safe, DISCARD reads MAPQ=0 before calculating percentage).
run_MAPQthreshold.sh

#8 Remove unmapped, secondary and supplementary alignments (-F 2308);  not properly paired reads -f 2; and retain only with a mapping quality (MAPQ) >30 (-q 30)
#  Sort and index bam file
run_filter_bam.sh

#9 count how many reads (percentage) retained after filtering
run_readsRetained.sh
