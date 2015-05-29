#!/usr/bin/env bash
### Collections introns across 3000 SRA samples
### Should be run on Rail-RNA v0.1.6a
### Uses manifest file from version of SRAmetadata repo tagged "intron_detection_experiment"
INDIR=s3://rail-eu-west-1/SRA3kprep
OUTDIR=s3://rail-eu-west-1/SRA3kintronexperiment1kcores
SRAMETADATAREPO=/Users/eterna/SRAmetadata # Location of SRA metadata repo. Its HEAD should be at "intron_detection_experiment".

rail-rna align elastic -m $SRAMETADATAREPO/sample_manifest_file/sample_size_3000.txt -o $OUTDIR -i $INDIR -a hg19 --region eu-west-1 -c 250 --core-instance-type c3.8xlarge --master-instance-type c3.8xlarge --core-instance-bid-price 0.52 --master-instance-bid-price 0.52 --no-consistent-view --deliverables itn,idx --max-task-attempts 6 --ec2-key-name raileuwest1