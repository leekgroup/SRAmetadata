#!/usr/bin/env python
"""
number_of_samples.py

Counts number of samples in which junctions appear.
"""
import sys

for line in sys.stdin:
	tokens = line.strip().split('\t')
	print '\t'.join(tokens[:3] + [str(tokens[-2].count(',') - 1)])