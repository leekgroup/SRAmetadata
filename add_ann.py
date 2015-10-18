#!/usr/bin/env python
"""
add_ann.py
Abhi Nellore
10/18/2015

Adds fields to all_SRA_introns.tsv.gz, which is read from stdin.
annotation files from command-line parameters and annotation from
GTF files specified as arguments of --annotations; writes to stdout.
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

Fields added:
x,y,z ,where each is 1 or 0. x indicates whether the 5' splice site is in
annotation, y indicates whether the 3' splice site is in annotation, and
z indicates whether the junction is in annotation.

We executed:
gzip -cd all_SRA_introns.tsv.gz | pypy add_ann.py >[output_file]
"""
import os
import sys
import subprocess

if __name__ == '__main__':
    import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    # Add command-line arguments
    parser.add_argument('--extract-splice-sites-path', type=str,
        default=os.path.join(os.path.dirname(__file__),
                                'extract_splice_sites.py'),
        help=('path to extract_splice_sites.py from HISAT v0.1.6-beta.'))
    parser.add_argument('--annotations', type=str, required=True, nargs='+',
        help='paths to GTF files encoding known junctions')
    args = parser.parse_args()

    annotated_junctions = set()
    annotated_5p = set()
    annotated_3p = set()
    refs = ['chr' + str(i) for i in xrange(1, 23)] + ['chrM', 'chrX', 'chrY']
    for annotation in args.annotations:
        extract_process = subprocess.Popen([sys.executable,
                                                args.extract_splice_sites_path,
                                                annotation],
                                                stdout=subprocess.PIPE)
        for line in extract_process.stdout:
            tokens = line.strip().split('\t')
            tokens[1] = int(tokens[1]) + 2
            tokens[2] = int(tokens[2])
            if not tokens[0].startswith('chr'):
                tokens[0] = 'chr' + tokens[0]
            if tokens[0] in refs:
                annotated_junctions.add(tuple(tokens[:-1]))
                if tokens[3] == '+':
                    annotated_5p.add(tuple(tokens[:-2]))
                    annotated_3p.add(tuple(tokens[0], tokens[2]))
                else:
                    assert tokens[3] == '-'
                    annotated_5p.add(tuple(tokens[0], tokens[2]))
                    annotated_3p.add(tuple(tokens[:-2]))
        extract_process.stdout.close()
        exit_code = extract_process.wait()
        if exit_code != 0:
            raise RuntimeError(
                'extract_splice_sites.py had nonzero exit code {}.'.format(
                                                                    exit_code
                                                                )
            )

    for line in sys.stdin:
        tokens = line.strip().split('\t')
        left = int(tokens[1])
        right = int(tokens[2])
        if tokens[3] == '+':
            threep = right
            fivep = left
        else:
            fivep = right
            threep = left
        if (tokens[0], fivep) in annotated_5p:
            x = '1'
        else:
            x = '0'
        if (tokens[0], threep) in annotated_3p:
            y = '1'
        else:
            y = '0'
        if (tokens[0], left, right) in annotated_junctions:
            z = '1'
        else:
            z = '0'
        print '\t'.join([line.strip(), x, y, z])