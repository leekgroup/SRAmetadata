#!/usr/bin/env bash
### Collections introns across 3000 SRA samples
### Should be run on Rail-RNA v0.1.6a
### Uses manifest file from version of SRAmetadata repo tagged "intron_detection_experiment"
INDIR=s3://rail-eu-west-1/SRA3kp
OUTDIR=s3://rail-eu-west-1/SRA3ka4
SRAMETADATAREPO=/Users/eterna/SRAmetadata # Location of SRA metadata repo. Its HEAD should be at "intron_detection_experiment".

python ~/rail/src align elastic -m $SRAMETADATAREPO/sample_manifest_file/sample_size_3000.txt -o $OUTDIR -i $INDIR -a hg19 --region eu-west-1 -c 1500 --core-instance-type c3.2xlarge --master-instance-type c3.2xlarge --core-instance-bid-price 0.13 --master-instance-bid-price 0.13 --no-consistent-view --deliverables itn,idx --max-task-attempts 6 --ec2-key-name raileuwest1 --json