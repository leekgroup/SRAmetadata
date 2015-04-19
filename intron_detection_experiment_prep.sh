#!/usr/bin/env bash
### Preprocesses 3000 SRA samples to obtain introns and their "initial" coverage across samples
### Should be run on version of Rail-RNA repo tagged "intron_detection_experiment"
### Uses manifest file from version of SRAmetadata repo tagged "intron_detection_experiment"
RAIL_EXE=/Users/eterna/rail/src # The "executable" for Rail is its source directory
OUT=s3://rail-eu-west-1/prep3000
SRAMETADATAREPO=/Users/eterna/SRAmetadata # Location of SRA metadata repo. Its HEAD should be at "intron_detection_experiment".

python $RAIL_EXE prep elastic -m $SRAMETADATAREPO/sample_manifest_file/sample_size_3000.txt -o $OUT --region eu-west-1 -c 20 --core-instance-bid-price 0.11 --master-instance-bid-price 0.11 --no-consistent-view --do-not-check-manifest