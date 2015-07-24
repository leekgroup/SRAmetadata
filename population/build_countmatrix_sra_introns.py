#!/usr/bin/python

# AUTOR: Jose Alquicira Hernandez
# PURPOSE: Get count matrix for introns
# DATE: Jan 10th, 2015

# PARAMETERS:
# 1st = tsv output file from Rail RNA

import sys

tsv = sys.argv[1] 
reads = []
n = 21504

file = open(tsv, 'r')
for l in file:
    l = l.strip('\n').split('\t')
    reads = map(int, l[7].split(','))
    if max(reads) >= 5 and float(sum(reads))/len(reads) > 0.95:
    	all_sample_reads = ['0'] * n
    	indexes = map(int, l[6].split(','))
    	j = 0
    	for i in indexes:
    		all_sample_reads[i] = reads[j]
    		j = j + 1
    	print l[0] + "_" + l[1] + "_" + l[2] + "_" + "("+ l[3]+ ")",
    	for r in all_sample_reads:
    		print "\t" + str(r),
    	print ""
file.close()
