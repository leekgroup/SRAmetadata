#!/usr/bin/env python

import sys
relevant_fields = ['cell type', 'tissue', 'cell line', 'barcode',
	'flowcell', 'lane', 'library id', 'platform', 'seqc sample', 'site',
	'treatment', 'patient id', 'gender', 'age', 'subtype', 'disease state']
all_fields = {}
for line in sys.stdin:
	if '!Sample_characteristics_ch1' in line:
		rest, _ content = line.partition(':')[1].strip()
		field = rest.rpartition('=')[1].strip()
		all_fields[field] = content

to_print = []
for field in relevant_fields:
	if field in all_fields:
		to_print.append(all_fields[field])
	else:
		to_print.append('NA')

print '\t'.join(to_print)