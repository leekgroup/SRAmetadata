#!/usr/bin/env python
"""
heatmap.py
Abhi Nellore
9/4/2015

Constructs matrix for heatmap with SRA study accession ID on one axis and
sample filter threshold on y axis.

Reads introns from stdin (all_SRA_introns.tsv.gz); writes to stdout.
all_SRA_introns.tsv.gz should have the following tab-separated fields on each
line:
1. chromosome
2. start position
3. end position
4. strand
5. start motif (e.g., GT)
6. end motif (e.g., AG)
7. comma-separated list of sample indexes in which junction was found
8. comma-separated list of numbers of reads in each corresponding sample from
    field 7

We executed:
gzip -cd all_SRA_introns.tsv.gz | pypy heatmap.py \
    -s index_to_SRA_accession.tsv >table_for_heatmap.tsv
"""
import sys
from collections import defaultdict

if __name__ == '__main__':
    import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    # Add command-line arguments
    parser.add_argument('-s', '--sra', type=str,
        default='index_to_SRA_accession.tsv',
        help=('TSV file whose tab-separated fields are sample index, '
              'project accession number, and other accession numbers, '
              'respectively')
        )
    parser.add_argument('--threshold-interval', type=float,
        default=.005,
        help='space between successive proportions of samples in which '
             'junction is required to appear')
    parser.add_argument('--min-threshold', type=float,
        default=0,
        help='minimum threshold from description of --threshold-interval')
    parser.add_argument('--max-threshold', type=float,
        default=.1,
        help='maximum threshold from description of --threshold-interval')
    parser.add_argument('--sample-count', type=int,
        default=21504,
        help='total number of SRA samples analyzed')
    args = parser.parse_args()

    index_to_project = {}
    with open(args.sra) as sra_stream:
        for line in sra_stream:
            tokens = line.strip().split('\t')
            index_to_project[tokens[0]] = tokens[1]

    projects = defaultdict(dict)
    for line in sys.stdin:
        tokens = line.strip().split('\t')
        samples = tokens[-2].strip().split(',')
        sample_count = len(samples)
        projects_to_increment = set([index_to_project[sample]
                                        for sample in samples])
        for project in projects_to_increment:
            try:
                projects[project][sample_count] += 1
            except KeyError:
                projects[project][sample_count] = 1

    # Dump results; first comes the header line
    thresholds = [args.threshold_interval * i + args.min_threshold
        for i in xrange(int((args.max_threshold - args.min_threshold)
                                / args.threshold_interval) + 1)]
    sample_thresholds = [round(float(args.sample_count) * threshold)
                            for threshold in thresholds]
    print ('\t'.join([''] + ['%.5f' % threshold for threshold in thresholds]))
    for project in projects:
        print '\t'.join([project]
                        + [str(sum(intron_count
                                for sample_count, intron_count
                                in projects[project].items()
                                if sample_count >= threshold))
                            for threshold in sample_thresholds])
