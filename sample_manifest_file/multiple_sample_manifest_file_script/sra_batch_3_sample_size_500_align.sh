#!/usr/bin/env bash
rail-rna align elastic -m sra_batch_3_sample_size_500.txt -i s3://rail-eu-west-1/sra_batch_3_sample_size_500_prep -o s3://rail-eu-west-1/sra_batch_3_sample_size_500_align -a hg19 --isofrag-idx s3://rail-eu-west-1/SRA3kat/cross_sample_results/isofrags.tar.gz --region eu-west-1 -c 60 --core-instance-type c3.8xlarge --master-instance-type c3.8xlarge --core-instance-bid-price 0.60 --master-instance-bid-price 0.60 --no-consistent-view --deliverables bed,tsv,bw --max-task-attempts 6 --ec2-key-name raileuw1
