#!/usr/bin/env python
"""
ann.py
Abhi Nellore
9/8/2015

Writes stats on proportions of detected junctions found in annotation. Requires
extract_splice_sites.py from HISAT v0.1.6-beta.

Reads introns from stdin (all_SRA_introns.tsv.gz), 
annotation files from command-line parameters and annotation from
GTF file specified as argument of --annotation; writes to stdout.
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

Tab-separated output written to stdout:
1. "project" or "sample", depending on whether junctions under consideration
are across a project or from a given sample
2. SRP (project) or SRR (run-level, but called "sample" here) accession number
3. number of junctions found across project or sample
4. number of junctions from project or sample found in annotation
5. number of reads overlapping junctions in project or sample
6. number of reads from project or sample overlapping annotated junctions
7. proportion of junctions from project or sample overlapping annotated
    junctions
8. proportion of reads from project or sample overlapping annotated junctions
9- SHARQ metadata, if available for SRR number

We executed:
gzip -cd all_SRA_introns.tsv.gz | pypy rank.py \
    -s index_to_SRA_accession.tsv --annotation [GTF file]
    >[output file]
"""
import sys
from collections import defaultdict
import os
import subprocess

if __name__ == '__main__':
    import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    # Add command-line arguments
    parser.add_argument('-s', '--sra', type=str,
        default=os.path.join(os.path.dirname(__file__),
                                'index_to_SRA_accession.tsv'),
        help=('TSV file whose tab-separated fields are sample index, '
              'project accession number, and other accession numbers, '
              'respectively')
        )
    parser.add_argument('--extract-splice-sites-path', type=str,
        default=os.path.join(os.path.dirname(__file__),
                                'extract_splice_sites.py'),
        help=('path to extract_splice_sites.py from HISAT v0.1.6-beta.'))
    parser.add_argument('--annotations', type=str, required=True, nargs='+',
        help='paths to GTF files encoding known junctions')
    parser.add_argument('--sharq', type=str, required=False,
        default=None,
        help='path to SHARQ metadata if available; '
             'we used sra-all-fields-2015-9-17.txt')
    args = parser.parse_args()

    index_to_project, index_to_sample, sample_to_metadata = {}, {}, {}
    with open(args.sra) as sra_stream:
        for line in sra_stream:
            tokens = line.strip().split('\t')
            index_to_project[tokens[0]] = tokens[1]
            index_to_sample[tokens[0]] = tokens[4]

    if args.sharq is not None:
        with open(args.sharq) as sharq_stream:
            header = sharq_stream.readline().partition(',')[2].replace(
                        ',', '\t'
                    ).strip()
            for line in sharq_stream:
                partitioned = line.partition(',')
                sample_to_metadata[partitioned[0]] = partitioned[2].replace(
                        ',', '\t'
                    ).strip()

    annotated_junctions = set()
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
        extract_process.stdout.close()
        exit_code = extract_process.wait()
        if exit_code != 0:
            raise RuntimeError(
                'extract_splice_sites.py had nonzero exit code {}.'.format(
                                                                    exit_code
                                                                )
            )

    (project_junctions_ann, project_reads_ann, sample_junctions_ann,
        sample_reads_ann, project_junctions, project_reads,
        sample_junctions, sample_reads) = [defaultdict(int) for _ in xrange(8)]
    for line in sys.stdin:
        tokens = line.strip().split('\t')
        junction = tuple([tokens[0], int(tokens[1]), int(tokens[2])])
        if junction in annotated_junctions:
            annotated = True
        else:
            annotated = False
        samples = tokens[-2].strip().split(',')
        sample_increments = [int(increment) for increment
                                in tokens[-1].strip().split(',')]
        project_increments = defaultdict(int)
        for i, sample_index in enumerate(samples):
            sample = index_to_sample[sample_index]
            project_increments[index_to_project[sample_index]] \
                += sample_increments[i]
            sample_junctions[sample] += 1
            sample_reads[sample] += sample_increments[i]
            if annotated:
                sample_junctions_ann[sample] += 1
                sample_reads_ann[sample] += sample_increments[i]
        for project in project_increments:
            project_junctions[project] += 1
            project_reads[project] += project_increments[project]
            if annotated:
                project_junctions_ann[project] += 1
                project_reads_ann[project] += project_increments[project]

    # Dump results
    if args.sharq is None:
        print '\t'.join(['type', 'accession', 'junctions',
                            'annotated junctions', 'overlap instances',
                            'annotated overlap instances',
                            'proportion of junctions '
                            'that are annotated',
                            'proportion of overlap instances '
                            'that are annotated'])
        for sample in sample_junctions:
            print '\t'.join(['sample', sample, str(sample_junctions[sample]),
                                str(sample_junctions_ann[sample]),
                                str(sample_reads[sample]),
                                str(sample_reads_ann[sample]),
                                '%.10f' % (float(sample_junctions_ann[sample])
                                            / sample_junctions[sample]),
                                '%.10f' % (float(sample_reads_ann[sample])
                                            / sample_reads[sample])])
        for project in project_junctions:
            print '\t'.join(['project', project,
                            str(project_junctions[project]),
                            str(project_junctions_ann[project]),
                            str(project_reads[project]),
                            str(project_reads_ann[project]),
                            '%.10f' % (float(project_junctions_ann[project])
                                        / project_junctions[project]),
                            '%.10f' % (float(project_reads_ann[project])
                                        / project_reads[project])])
    else:
        print '\t'.join(['type', 'accession', 'junctions',
                            'annotated junctions', 'overlap instances',
                            'annotated overlap instances',
                            'proportion of junctions '
                            'that are annotated',
                            'proportion of overlap instances '
                            'that are annotated', header])
        na_line = '\t'.join(['NA']*(header.count('\t')+1))
        for sample in sample_junctions:
            print '\t'.join(['sample', sample, str(sample_junctions[sample]),
                                str(sample_junctions_ann[sample]),
                                str(sample_reads[sample]),
                                str(sample_reads_ann[sample]),
                                '%.10f' % (float(sample_junctions_ann[sample])
                                            / sample_junctions[sample]),
                                '%.10f' % (float(sample_reads_ann[sample])
                                            / sample_reads[sample]),
                                na_line if sample not in sample_to_metadata
                                else sample_to_metadata[sample]])
        for project in project_junctions:
            print '\t'.join(['project', project,
                            str(project_junctions[project]),
                            str(project_junctions_ann[project]),
                            str(project_reads[project]),
                            str(project_reads_ann[project]),
                            '%.10f' % (float(project_junctions_ann[project])
                                        / project_junctions[project]),
                            '%.10f' % (float(project_reads_ann[project])
                                        / project_reads[project]), na_line])
