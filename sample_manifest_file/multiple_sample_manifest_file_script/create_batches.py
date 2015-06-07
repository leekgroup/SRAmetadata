#!/usr/bin/env python
"""
create_batches.py

Divides all-of-Illumina human RNA-seq SRA manifest file into batches of 
approximately 500 reproducibly and creates scripts for running job flows on
EMR. Use this rather than the R script. Use Rail-RNA v0.1.7 to reproduce
results.
Usage: cat all_of_illumina_sra_manifest | python create_batches.py
"""
import random
import sys
import os
# This will need changing if the user seeks to reproduce results
_s3_bucket = 's3://rail-eu-west-1'
_region = 'eu-west-1'
_c3_2xlarge_bid_price = '0.13'
_c3_8xlarge_bid_price = '0.60'
_key_name = 'raileuw1'
_isofrag_idx_path = (
        's3://rail-eu-west-1/SRA3kat/cross_sample_results/isofrags.tar.gz'
    ) # from intron detection experiments

lines = sys.stdin.read().strip().split('\n')
random.seed(lines[0])
random.shuffle(lines)
manifests = [lines[i:i+500] for i in xrange(0, len(lines), 500)]
# Merge last two batches
manifests[-2] += manifests[-1]
manifests.pop()
for i, manifest in enumerate(manifests):
    manifest_filename = 'sra_batch_{}_sample_size_{}.txt'.format(
                                                            i, len(manifest)
                                                        )
    prep_output_dir = os.path.join(
                        _s3_bucket,
                        'sra_batch_{}_sample_size_{}_prep'.format(
                                                        i, len(manifest)
                                                    )
                    )
    align_output_dir = os.path.join(
                        _s3_bucket,
                        'sra_batch_{}_sample_size_{}_align'.format(
                                                        i, len(manifest)
                                                    )
                    )
    with open(manifest_filename, 'w') as manifest_stream:
        manifest_stream.write('\n'.join(manifest) + '\n')
    with open('sra_batch_{}_sample_size_{}_prep.sh'.format(
                                                        i, len(manifest)
                                                    ), 'w') as prep_stream:
        prep_stream.write(
"""#!/usr/bin/env bash
rail-rna prep elastic -m {manifest} -o {output_dir} -c 20 --region {region} --master-instance-type c3.2xlarge --core-instance-type c3.2xlarge --no-consistent-view --master-instance-bid-price {bid_price} --core-instance-bid-price {bid_price} --ec2-key-name {key_name} --do-not-check-manifest --do-not-bin-quals --skip-bad-records --max-task-attempts 8
""".format(manifest=manifest_filename, output_dir=prep_output_dir, region=_region, 
            bid_price=_c3_2xlarge_bid_price, key_name=_key_name)
)
    with open('sra_batch_{}_sample_size_{}_align.sh'.format(
                                                        i, len(manifest)
                                                    ), 'w') as align_stream:
        align_stream.write(
"""#!/usr/bin/env bash
rail-rna align elastic -m {manifest} -i {input_dir} -o {output_dir} -a hg19 --isofrag-idx {isofrag_idx} --region {region} -c 60 --core-instance-type c3.8xlarge --master-instance-type c3.8xlarge --core-instance-bid-price {bid_price} --master-instance-bid-price {bid_price} --no-consistent-view --deliverables bed,tsv,bw --max-task-attempts 6 --ec2-key-name {key_name}
""".format(manifest=manifest_filename, input_dir=prep_output_dir,
            output_dir=align_output_dir, region=_region, 
            bid_price=_c3_8xlarge_bid_price, key_name=_key_name,
            isofrag_idx=_isofrag_idx_path)
)

