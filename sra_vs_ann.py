#!/usr/bin/env python
"""
sra_vs_ann.py

Generates TSVs for plotting SRA vs annotation.

Reads introns from stdin (all_SRA_introns.tsv.gz), 
annotation files from command-line parameters and annotation from
GTF file(s) specified as argument(s) of --annotations; writes to stdout.
all_SRA_introns.tsv.gz should have AT LEAST the following tab-separated fields
on each line:
1. chromosome
2. start position
3. end position
...
next-to-last field: comma-separated list of sample indexes
last field: comma-separated list of read coverages corresponding to sample
    indexes

Output:
    [basename].coverage.tsv, with the following tab-separated fields:
    1. sample index
    2. number of annotated junction overlaps in sample
    3. number of unannotated junction overlaps in sample
    4. total number of junction overlaps (sum of fields 2 and 3)
    5. number of annotated junctions found in sample
    6. number of unannotated junctions found in sample
    7. total number of junctions (sum of fields 5 and 6)

    [basename].ann.tsv, with the following tab-separated fields:
    1. number of samples
    2. number of annotated junctions found in >= (field 1) samples
    3. number of unannotated junctions found in >= (field 1) samples
    4. total number of junctions found in >= (field 1) samples
    5. total number of junctions in annotation (same for all lines)

We used the annotations:
    gencode.v19.annotation.gtf.gz
        (from ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/)
    Homo_sapiens.GRCh37.75.gtf.gz
        (from ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/)
    refGene.gtf
        (downloaded according to instructions at 
         https://groups.google.com/a/soe.ucsc.edu/forum/#!msg/genome/
         bYEoa_hrSiI/cJ8WjnqXhlIJ ; uses command
         mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -N \ 
          -e "select * from refGene;" hg19 | cut -f2- | genePredToGtf
          file stdin refGene.gtf

    see also: http://genomewiki.ucsc.edu/index.php/Genes_in_gtf_or_gff_format)

We first unzipped all gtf.gz files above and then executed:
    1) gzip -cd all_SRA_introns.tsv.gz | pypy sra_vs_ann.py
        --annotation gencode.v19.annotation.gtf --basename gencode_v19
    2) gzip -cd all_SRA_introns.tsv.gz | pypy sra_vs_ann.py
        --annotation Homo_sapiens.GRCh37.75.gtf --basename ensembl_v75
    3) gzip -cd all_SRA_introns.tsv.gz | pypy sra_vs_ann.py
        --annotation refGene.gtf --basename refseq
    4) gzip -cd all_SRA_introns.tsv.gz | pypy sra_vs_ann.py
        --annotation refGene.gtf gencode.v19.annotation.gtf
        Homo_sapiens.GRCh37.75.gtf --basename union
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
    parser.add_argument('--extract-splice-sites-path', type=str,
        default=os.path.join(os.path.dirname(__file__),
                                'extract_splice_sites.py'),
        help=('path to extract_splice_sites.py from HISAT v0.1.6-beta.'))
    parser.add_argument('--annotations', type=str, required=True, nargs='+',
            help='space-separated paths to GTF files encoding known junctions'
        )
    parser.add_argument('--basename', type=str, required=True,
            help='basename for output files'
        )
    args = parser.parse_args()

    annotated_junctions = set()
    for annotation in args.annotations:
        extract_process = subprocess.Popen([sys.executable,
                                                args.extract_splice_sites_path,
                                                annotation],
                                                stdout=subprocess.PIPE)
        for line in extract_process.stdout:
            tokens = line.strip().split('\t')
            tokens[1] = str(int(tokens[1]) + 2)
            tokens[2] = str(int(tokens[2]))
            annotated_junctions.add(tuple(tokens[:-1]))
        extract_process.stdout.close()
        exit_code = extract_process.wait()
        if exit_code != 0:
            raise RuntimeError(
                'extract_splice_sites.py had nonzero exit code {}.'.format(
                                                                    exit_code
                                                                )
            )

    annotated_coverage = defaultdict(int)
    unannotated_coverage = defaultdict(int)
    sample_read_annotated = defaultdict(int)
    sample_junction_annotated = defaultdict(int)
    sample_read_unannotated = defaultdict(int)
    sample_junction_unannotated = defaultdict(int)
    for line in sys.stdin:
        tokens = line.strip().split('\t')
        coverages = [int(el) for el in tokens[-1].split(',')]
        sample_count = len(coverages)
        if tuple(tokens[:3]) in annotated_junctions:
            annotated_coverage[sample_count] += 1
            for index, sample_index in enumerate(
                    [int(el) for el in tokens[-2].split(',')]
                ):
                sample_read_annotated[sample_index] += coverages[index]
                sample_junction_annotated[sample_index] += 1
        else:
            unannotated_coverage[sample_count] += 1
            for index, sample_index in enumerate(
                    [int(el) for el in tokens[-2].split(',')]
                ):
                sample_read_unannotated[sample_index] += coverages[index]
                sample_junction_unannotated[sample_index] += 1
    sample_indexes = sorted(list(set(sample_read_unannotated.keys()
                                + sample_read_annotated.keys()
                                + sample_junction_annotated.keys()
                                + sample_junction_unannotated.keys())))
    with open(args.basename + '.coverage.tsv', 'w') as coverage_stream:
        for sample_index in sample_indexes:
            print >>coverage_stream, '%d\t%d\t%d\t%d\t%d\t%d\t%d' % (
                    sample_index, sample_read_annotated[sample_index],
                    sample_read_unannotated[sample_index],
                    sample_read_annotated[sample_index]
                    + sample_read_unannotated[sample_index],
                    sample_junction_annotated[sample_index],
                    sample_junction_unannotated[sample_index],
                    sample_junction_annotated[sample_index]
                    + sample_junction_unannotated[sample_index]
                )
    max_coverage = max(annotated_coverage.keys()
                        + unannotated_coverage.keys())
    annotated_junction_total, unannotated_junction_total = 0, 0
    total_annotated_junctions = len(annotated_junctions)
    with open(args.basename + '.ann.tsv', 'w') as ann_stream:
        for coverage in xrange(max_coverage, 0, -1):
            annotated_junction_total += annotated_coverage[coverage]
            unannotated_junction_total += unannotated_coverage[coverage]
            print >>ann_stream, '%d\t%d\t%d\t%d\t%d' % (
                    coverage, annotated_junction_total,
                    unannotated_junction_total,
                    annotated_junction_total + unannotated_junction_total,
                    total_annotated_junctions
                )
