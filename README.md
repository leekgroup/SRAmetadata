# SRAmetadata

This repository contains code, documentation and presentations on how to get metadata from SRA (Sequence Read Archive).

---

# Getting and cleaning data

The Sqlite database for SRA was downloaded from the website:
http://gbnci.abcc.ncifcrf.gov/backup/SRAmetadb.sqlite.gz

```
wget http://gbnci.abcc.ncifcrf.gov/backup/SRAmetadb.sqlite.gz
--2015-04-17 00:58:29--  http://gbnci.abcc.ncifcrf.gov/backup/SRAmetadb.sqlite.gz
Resolving gbnci.abcc.ncifcrf.gov... 129.43.40.100, 2607:f220:41d:4f4d::92
Connecting to gbnci.abcc.ncifcrf.gov|129.43.40.100|:80... connected.
HTTP request sent, awaiting response... 200 OK
Length: 836509538 (798M) [application/x-gzip]
Saving to: “SRAmetadb.sqlite.gz”

100%[======================================================================================================>] 836.509.538 20,3M/s   in 35s     

2015-04-17 00:59:17 (22,8 MB/s) - “SRAmetadb.sqlite.gz” saved [836509538/836509538]
```

```
gunzip SRAmetadb.sqlite.gz
```

---

## Create table for metadata and manifest file


The script __define_and_get_fields_SRA.R__ creates a table with relevant metadata from SRA and an associated manifest file. Both files can be 
linked, since the order of rows is the same.

```
Output files: 
- __all_illumina_sra_for_human.txt__
- __manifest_file_illumina_sra_human__


---

## Generate random sample

Then, a sample of 3000 runs without replacement was made with the script __sample_manifest_file.R__.
This script generates a file with the sample and a second file that maps the column number from the manifest file to the metadata fields 
in "sample_size_3000.txt". The value of the seed used to make the sample was 42.

Output files:
- __relationship_manifest_file-sample__
- __sample_size_3000.txt__

# To-do list:
Remove these fields from SRA metadata if present (they do not have registers):

SRR_bamFile

sequence_space

SRX_bamFile

SRX_fastqFTP

run_url_link

run_entrez_link

adapter_spec

instrument_name

sequence_space

base_caller

quality_scorer

number_of_levels

multiplier

qtype

experiment_entrez_link

anonymized_name

individual_name

sample_entrez_link

study_entrez_link



+ Have a controlled vocabulary for not reported data: "NAs"

* design_description:
 none provided,
 N/A

* library_selection:
 unespecified

* library_construction_protocol:
 none provided

# Notes
library_name has information about biological and technical replicates (manual curation??)
