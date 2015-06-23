#!/usr/bin/env bash
### Preprocesses 3000 SRA samples to obtain introns and their "initial" coverage across samples
### Should be run on Rail-RNA v0.1.7a
### Uses manifest file from version of SRAmetadata repo tagged "intron_detection_experiment"
OUTDIR=s3://rail-eu-west-1/SRA3kp
SRAMETADATAREPO=/Users/eterna/SRAmetadata # Location of SRA metadata repo. Its HEAD should be at "intron_detection_experiment".

rail-rna prep elastic -m $SRAMETADATAREPO/sample_manifest_file/sample_size_3000.txt -o $OUTDIR -c 20 --region eu-west-1 --master-instance-type c3.2xlarge --core-instance-type c3.2xlarge --no-consistent-view --master-instance-bid-price 0.11 --core-instance-bid-price 0.11 --ec2-key-name raileuwest1 --do-not-check-manifest --do-not-bin-quals --skip-bad-records --max-task-attempts 8