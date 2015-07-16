#!/usr/bin/env bash
# Downloads all files with collected_introns results from all-of-SRA runs and merges them into one file
# This script requires the AWS CLI, PyPy
# Used PyPy v2.4.0
# $1: path to destination directory
# $2: path to Bowtie index for hg19
MANIFESTS=/scratch0/langmead-fs1/SRAmetadata/sample_manifest_file/multiple_sample_manifest_file_script # Path to batch manifest files
BUCKET=s3://rail-eu-west-1 # Where batch results ended up
OUTPUTDIR=$1
BOWTIEIDX=$2
mkdir -p ${OUTPUTDIR}
cd $1
# Download
for i in {0..41}; do aws s3 cp ${BUCKET}/sra_batch_${i}_sample_size_500_itn/collected_introns/collected_introns.tsv.gz ./batch_${i}.tsv.gz; done
# Exception: last batch has sample size 506
aws s3 cp ${BUCKET}/sra_batch_42_sample_size_506_itn/collected_introns/collected_introns.tsv.gz ./batch_42.tsv.gz
# Sort files, prepending batch number to each line
for i in {0..42}; do (gzip -cd batch_${i}.tsv.gz | sort -k1,1 -k2,2n -k3,3n | awk -v r=${i} '{print r "\t" $0}' | gzip >batch_${i}.sorted.tsv.gz) & done
# Merge sorted files, adjusting sample indexes so they're (batch number * 500 + original sample index).
# Note that two sample indexes will be missing because they were removed from batch manifest files. See NOTES for more information.
cmd="sort -m -k2,2 -k3,3n -k4,4n"; for i in {0..42}; do cmd="$cmd <(gzip -cd batch_${i}.sorted.tsv.gz) "; done
eval "$cmd" | gzip >unmerged_intron_lines.tsv.gz
# Note the --fix-batch-9 command-line parameter below! Batch 9 was preprocessed with the manifest file sra_batch_9_sample_size_500.txt
# but aligned with the manifest file sra_batch_9_sample_size_500_old.txt, which has an extra sample that wasn't found on the server.
# (See NOTES for which sample it was and how it was removed.) The --fix-batch-9 command-line parameter uses the old manifest file
# to extract the proper sample indexes.
gzip -cd unmerged_intron_lines.tsv.gz | pypy ${MANIFESTS}/combine.py --bowtie-idx ${BOWTIEIDX} --fix-batch-9 | sort -k1,1 -k2,2n -k3,3n | gzip >all_SRA_introns.tsv.gz