#!/usr/bin/env bash
# Miscellaneous commands used to analyze introns. Requires mawk 1.3.4 2015-05-03
cat all_SRA_introns.tsv | ./mawk-1.3.4-20150503/mawk -F ',' '{print $0 "\t" (NF - 1)/2 + 1}' | ./mawk-1.3.4-20150503/mawk '{print $1 "\t" $2 "\t" $3 "\t" $9}' >introns_and_sample_counts.tsv
cat all_SRA_introns.tsv | ./mawk-1.3.4-20150503/mawk '{print $1 "\t" $2 "\t" $3 "\t," $8}' | ./mawk-1.3.4-20150503/mawk -F ',' '{sum=0;for(i=2;i<=NF;i++){sum += $i} print $1 sum}' >introns_and_coverages.tsv
sort -k4,4nr -o introns_and_sample_counts.tsv introns_and_sample_counts.tsv
sort -k4,4nr -o introns_and_sample_counts.tsv introns_and_coverages.tsv