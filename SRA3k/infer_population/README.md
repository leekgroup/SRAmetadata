
# Get Files from 1000 genomes project

```
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/phase1_sites_missing_in_phase3/*

mkdir indels_data_reference

mv phase1_sites_missing_in_phase3.20141215.vcf* indels_data_reference

gzip -d phase1_sites_missing_in_phase3.20141215.vcf.gzip

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5a.20130502.sites.vcf.gz -o download_allwgs

cat <(grep "#" ALL.wgs.phase3_shapeit2_mvncall_integrated_v5a.20130502.sites.vcf) 
<(grep "INDEL" ALL.wgs.phase3_shapeit2_mvncall_integrated_v5a.20130502.sites.vcf) > ALL_wgs_phase3_indels.vcf
```

