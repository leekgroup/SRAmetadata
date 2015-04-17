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


The script __define_and_get_fields_SRA.R__ creates a table with relevant metadata from SRA and an associated manifest file. Both files can be 
linked, since the order of rows is the same.
 
Output file: 
- __all_illumina_sra_for_human.txt__
- __manifest_file_illumina_sra_human__

Then, a sample of 3000 runs without replacement was made with the script __sample_manifest_file.R__.
This script generates a file with the sample and a second file that maps the column number from the manifest file to the metadata fields 
in "sample_size_3000.txt"
