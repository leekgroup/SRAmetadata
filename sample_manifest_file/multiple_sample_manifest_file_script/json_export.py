#!/usr/bin/env python
"""
json_export.py

Writes extended JSON file for loading by mongodb, one line per document
(intron).
"""
import sys

for line in sys.stdin:
    (chrom, start, end, strand,
        start_motif, end_motif,
        sample_indexes, coverages) = line.strip().split('\t')
    sample_indexes = sample_indexes.split(',')
    coverages = coverages.split(',')
    if end_motif == 'AG':
        if start_motif == 'GT':
            motif_number = 0
        elif start_motif == 'GC':
            motif_number = 1
    else:
        assert end_motif == 'AC'
        assert start_motif == 'AT'
        motif_number = 2
    print ('{{ chrom: {chrom}, start: {start}, end: {end}, forward: {forward}, '
           'motif: {motif}, coverages : [ {coverages} ]}}').format(
                    chrom=chrom,
                    start=start,
                    end=end,
                    forward=('true' if strand == '+' else 'false'),
                    motif=motif_number,
                    coverages=', '.join(('%s : %s' % (
                                                sample_indexes[i], coverages[i]
                                            )) for i in xrange(
                                                        len(sample_indexes)
                                                    ))
                )