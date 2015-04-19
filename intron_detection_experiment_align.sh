#!/usr/bin/env bash
### Collections introns across 3000 SRA samples
### Should be run on version of Rail-RNA repo tagged "intron_detection_experiment"
### Uses manifest file from version of SRAmetadata repo tagged "intron_detection_experiment"
RAIL_EXE=/Users/eterna/rail/src # The "executable" for Rail is its source directory
OUT=s3://rail-eu-west-1/collectintrons3000
SRAMETADATAREPO=/Users/eterna/SRAmetadata # Location of SRA metadata repo. Its HEAD should be at "intron_detection_experiment".

python $RAIL_EXE align elastic -m $SRAMETADATAREPO/sample_manifest_file/sample_size_3000.txt -o $OUT --region eu-west-1 -c 60 --core-instance-type c3.8xlarge --master-instance-type c3.8xlarge --core-instance-bid-price 0.35 --master-instance-bid-price 0.35 --no-consistent-view