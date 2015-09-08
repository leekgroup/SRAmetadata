#!/usr/bin/env python
"""
associate.py

Associates output of tag.py with SRR and SRP accession numbers from output
of ann.py. All input is read from stdin. cat in this order:
index_to_SRA_accession.tsv
output of tag.py
output of ann.py
"""
import sys
acc_to_sample = {}
sample_to_line = {}
captured_header_line = False
printed_header_line = False
for line in sys.stdin:
    tokens = line.strip().split('\t')
    if len(tokens) == 5 and tokens[1][1:3] == 'RP' and tokens[2][1:3] == 'RS':
        acc_to_sample[tokens[1]] = tokens[2]
        acc_to_sample[tokens[4]] = tokens[2]
    elif 'cell line' in line and 'small rna' in line and 'submission update' \
        in line and 'attributes' in line and not captured_header_line:
        header_line = line.strip()
        captured_header_line = True
    elif captured_header_line and len(tokens) > 10:
        sample_to_line[tokens[9]] = line.strip()
    elif not printed_header_line:
        print '\t'.join(['sample or project', 'accession', 'junction count',
                            'annotated junction count', 'read count',
                            'reads overlapping annotated junctions',
                            'proportion of junctions that are annotated',
                    'proportion of reads overlapping annotated junctions']) \
                + header_line
        printed_header_line = True
        print line.strip() + '\t' + sample_to_line[acc_to_sample[tokens[1]]]a
    else:
        print line.strip() + '\t' + sample_to_line[acc_to_sample[tokens[1]]]