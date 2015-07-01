"""
annotation_precision_and_recall.py

Computes the precision and recall of SRA samples' .

All data is read from stdin. Expects SRA intron lines to have greater than
5 columns and annotation introns to have at most 5 columns. The first three
columns of each input line must be the chromosome, start position, and end
position of a given intron, respectively. All columns should be tab-separated.

Tab-separated output fields:
1. Sample index
2. Precision of SRA sample's introns by annotation's introns
3. Recall of SRA sample's introns by annotation's introns
4. Recall of SRA sample's spliced reads by annotation's introns
"""

import sys
from collections import defaultdict

# Suppress output if number of introns in sample is < 1000
_SUPPRESS = 1000

annotation_introns = set()
overlaps, totals = defaultdict(int), defaultdict(int)
weighted_overlaps, weighted_totals = defaultdict(int), defaultdict(int)

for line in sys.stdin:
	tokens = line.strip().split('\t')
	chrom, start, end = tokens[0], int(tokens[1]), int(tokens[2])
	if len(tokens) <= 5:
		annotation_introns.add((chrom, start, end))
	else:
		samples = [int(sample) for sample in tokens[-2].split(',')]
		coverages = [int(coverage) for coverage in tokens[-1].split(',')]
		for i, sample in enumerate(samples):
			totals[sample] += 1
			weighted_totals[sample] += coverages[i]
		if (chrom, start, end) in annotation_introns:
			overlaps[sample] += 1
			weighted_overlaps[sample] += coverages[i]

annotation_intron_count = len(annotation_introns)

for sample in totals:
	if totals[sample] >= _SUPPRESS:
		print '\t'.join([
				str(sample), float(overlaps[sample]) / annotation_intron_count,
				float(overlaps[sample]) / totals[sample],
				float(weighted_overlaps[sample]) / weighted_totals[sample]
			])