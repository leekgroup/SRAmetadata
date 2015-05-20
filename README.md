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

```R
# Load library
library('RSQLite')
library('magrittr')

# Define functions
"%p%"  <- function(x, y) paste0(x, y)

# Update sql if necessary:
sqlfile <- file.path('..', 'SRAmetadb.sqlite') # Using the file from another location to avoid having copies of this 2.4 Gb file

# Create connection
sra_con <- dbConnect(SQLite(),sqlfile)

# Query database
query <- "SELECT  
study_accession, sample_accession, experiment_accession, 
run_accession, submission_accession,
taxon_id,
common_name,
anonymized_name,
individual_name,
description,
sample_attribute,
study_title,
study_type,
study_abstract,
center_project_name,
study_description,
related_studies,
primary_study,
study_attribute,
sradb_updated,
SRR_bamFile,
instrument_name,
run_date,
run_center,
experiment_name,
run_attribute,
SRX_bamFile,
SRX_fastqFTP,
study_name,
design_description,
sample_name,
library_name,
library_strategy,
library_source,
library_selection,
library_layout,
library_construction_protocol,
spots,
adapter_spec,
read_spec,
platform,
instrument_model,
platform_parameters,
sequence_space,
base_caller,
quality_scorer,
number_of_levels,
multiplier,
qtype,
experiment_attribute
FROM sra
WHERE platform = 'ILLUMINA' AND
library_strategy = 'RNA-Seq' AND
taxon_id = 9606;"

selected <- dbGetQuery(sra_con, query)
query <- paste0("SELECT run_accession, file_name, md5, bytes, sradb_updated FROM fastq WHERE run_accession IN ('", 
                paste(selected$run_accession, collapse="', '"), "')")
fastq_namefiles <- dbGetQuery(sra_con, query)

# Disconnect from db
dbDisconnect(sra_con)

metadata <- merge(selected, fastq_namefiles, by = "run_accession")

metadata$date_download <- rep(Sys.Date(),nrow(metadata))

# Rename column name for file data
names(metadata)[names(metadata) == 'sradb_updated.y'] <- 'sradb_updated_file'


get.fastq.urls <- function(df) {
    url <- rep(NA, length(df$run_accession))
    for(i in 1:length(df$run_accession)) {
        run <- df$run_accession[i]
        filename <- df$file_name[i]
        if(nchar(run) < 10) {
            url[i] <- file.path('ftp://ftp.sra.ebi.ac.uk/vol1/fastq',
                                substring(run, 1, 6), run, filename)
        } else {
            dir2 <- paste( c(rep(x='0', 12-nchar(run)), substring(run, 10, 
                                                                  nchar(run))), collapse = '' )
            url[i] <- file.path('ftp://ftp.sra.ebi.ac.uk/vol1/fastq', 
                                substring(run, 1, 6), dir2, run, filename)
        }  
    }
    return(url)
}


metadata$URL <- get.fastq.urls(metadata)
metadata$layout <- unlist(lapply(strsplit(metadata$library_layout, " - "), `[`, 1))


paired_and_forward <- intersect(grep("_1", metadata$file_name), 
                                grep("PAIRED", metadata$layout))

paired_and_reverse <- intersect(grep("_2", metadata$file_name), 
                                grep("PAIRED", metadata$layout))

ids_filename <- c("run_accession", "study_accession", "sample_accession", 
                  "experiment_accession", "submission_accession", "file_name", "URL", "md5")

PF <- metadata[paired_and_forward, ids_filename]
PR <- metadata[paired_and_reverse, ids_filename]

PF <- cbind(PF, paired_and_forward)
PR <- cbind(PR, paired_and_reverse)
paired_experiments <- merge(PF, PR, by = head(ids_filename, -3))

colnames(paired_experiments) <- c("run_accession", "study_accession", "sample_accession", "experiment_accession",
                                  "submission_accession", "file_name.f", "URL.f", "md5.f","paired_and_forward", "file_name.r",
                                  "URL.r", "md5.r","paired_and_reverse")


reverse_url <- rep("NA", nrow(metadata))
reverse_url[paired_experiments$paired_and_forward] <- paired_experiments$URL.r
metadata$URL_R <- reverse_url

reverse_md5 <- rep("NA", nrow(metadata))
reverse_md5[paired_experiments$paired_and_forward] <- paired_experiments$md5.r
metadata$md5_R <- reverse_md5

# cbind(metadata$URL, metadata$URL_R)[paired_and_forward,]
# cbind(metadata$URL, metadata$URL_R)[paired_and_reverse,]

metadata <- metadata[-paired_experiments$paired_and_reverse, ]

search.field <- function(column, field) {
    temp <- rep('NA', length(column))
    pos <- grep(paste0('.*', field, ':.*'), column, ignore.case = FALSE)
    temp[pos]<- as.numeric(pos)
    for(x in temp){
        if(x != 'NA'){
            x <- as.numeric(x)
            f <- sub(".*" %p% field %p% ": ", "", grep('.*' %p% field %p% ':.*',
                                                       unlist(strsplit(column[x], "||", fixed = TRUE)), 
                                                       perl = TRUE, value = TRUE, ignore.case = TRUE)[1], 
                     perl = TRUE,  
                     ignore.case = TRUE)
            temp[x]<- sub(" $",'', f)        
        } 
    }
    return(unlist(temp))
}



# Get cell type
metadata$cell_type <- search.field(metadata$sample_attribute, "cell type")

# Get tissue
metadata$tissue <- search.field(metadata$sample_attribute, "tissue")

# Get cell line
metadata$cell_line <- search.field(metadata$sample_attribute, "cell line")

# Get strain
metadata$strain <- search.field(metadata$sample_attribute, "Strain")

# Get age
metadata$age <- search.field(metadata$sample_attribute, "age")

# Get disease
metadata$disease <- search.field(metadata$sample_attribute, "disease")



# Get population
metadata$population <- sub('hapmap ','', search.field(metadata$sample_attribute, "population"), 
                           ignore.case = TRUE)
# Get race
metadata$race <- search.field(metadata$sample_attribute, "race")


metadata$population <- ifelse((metadata$population == 'NA') & (metadata$race == 'NA'),'NA',
                              ifelse(metadata$population == 'NA', metadata$race, metadata$population))

metadata$race <- NULL

# Get sex of individuals
metadata$sex <- search.field(metadata$sample_attribute, "sex") %>%
    sub('^male$', 'M', ., ignore.case = TRUE) %>% 
    sub('^female$', 'F', . , ignore.case = TRUE) %>%
    sub('not_documented', 'NA', . , fixed = TRUE) %>%
    sub('not determined', 'NA', . , fixed = TRUE) %>%
    sub('not applicable', 'NA', . , fixed = TRUE) %>%
    sub('Unknown', 'NA', . , fixed = TRUE) %>%
    sub('not collected', 'NA', . , fixed = TRUE) %>%
    sub('U', 'NA', . , fixed = TRUE) %>%
    sub('missing','NA', . , fixed = TRUE) %>%
    sub('asexual', 'NA', . , fixed = TRUE) %>%
    sub('mixed sex', 'B', . , fixed = TRUE) %>%
    sub('1 Male, 2 Female', 'M, FF', . , fixed = TRUE) %>%
    sub('2 Male, 1 Female', 'MM, F', . , fixed = TRUE) %>%
    sub('3 Female', 'FFF', . , fixed = TRUE)

# Get gender
metadata$gender <- search.field(metadata$sample_attribute, "gender") %>%
    sub('^male$', 'M', ., ignore.case = TRUE) %>% 
    sub('^female$', 'F', . , ignore.case = TRUE) %>%
    sub('MAL', 'M', . , fixed = TRUE) %>%
    sub('FEM', 'F', . , fixed = TRUE) %>%
    sub('XY', 'M', . , fixed = TRUE) %>%
    sub('XX', 'F', . , fixed = TRUE) %>%
    sub('not_documented', 'NA', . , fixed = TRUE)  %>%
    sub('--', 'NA', . , fixed = TRUE) %>%
    sub('3 Female', 'FFF', . , fixed = TRUE) %>%
    sub('N/A', 'NA', . , fixed = TRUE)


metadata$sex <- ifelse((metadata$sex == 'NA') & (metadata$gender == 'NA'),'NA',
                       ifelse(metadata$sex == 'NA', metadata$gender, metadata$sex))

metadata$gender <- NULL

# Cheking merged columns
# table(metadata$gender)[c(1,3,4)] + table(metadata$sex)[c(2,4,7)]

# Get source type of samples
metadata$source_name <- search.field(metadata$sample_attribute, "source_name")

# Removing new line characters
metadata$read_spec <- gsub("\n", "", metadata$read_spec)


# Create manifest file
labels <- c("study_accession", "sample_accession", 
            "experiment_accession", "run_accession")

sample_labels <- as.vector(apply(metadata[,labels], 1 ,
                                 function(x){paste( x ,collapse = "_")}))


s_i <- grep("SINGLE", metadata$layout)
p_i <- grep("PAIRED", metadata$layout)

single <- cbind(metadata$URL[s_i], metadata$md5[s_i], sample_labels[s_i] %p% "-1-1")
rownames(single) <- s_i
paired <- cbind(metadata$URL[p_i], metadata$md5[p_i], metadata$URL_R[p_i], metadata$md5_R[p_i], 
                sample_labels[p_i] %p% "-1-1")
rownames(paired) <- p_i
order_list <- c(p_i,s_i)
rownames(metadata) <- 1:nrow(metadata)

paired_as_single <- subset(paired, paired[,3] == 'NA' | paired[,4] == 'NA')
paired <- subset(paired, paired[,3] != 'NA' | paired[,4] != 'NA')
single <- rbind(single, paired_as_single[,c(1,2,5)])
order_list <- c(as.numeric(rownames(paired)),as.numeric(rownames(single)))

# metadata[p_i,]
# metadata[s_i,]
# metadata[as.numeric(rownames(paired)),]
# metadata[as.numeric(rownames(single)),]
# rownames(paired)
# rownames(single)

# Write table with all Illimina data
write.table(metadata[match(order_list, rownames(metadata)),], "all_illumina_sra_for_human.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(paired, "manifest_file_illumina_sra_human", sep = "\t", quote = FALSE,
            col.names = FALSE, row.names = FALSE)
write.table(single, "manifest_file_illumina_sra_human", sep = "\t", quote = FALSE, 
            col.names = FALSE, row.names = FALSE, append = TRUE)


# MANIFEST FILE FORMAT
# <FASTQ URL>(tab)<optional MD5>(tab)<sample label>
# <FASTQ URL 1>(tab)<optional MD5 1>(tab)<FASTQ URL 2>(tab)<optional MD5 2>(tab)<sample label>


```
Output files: 
- __all_illumina_sra_for_human.txt__
- __manifest_file_illumina_sra_human__


---

## Generate random sample

Then, a sample of 3000 runs without replacement was made with the script __sample_manifest_file.R__.
This script generates a file with the sample and a second file that maps the column number from the manifest file to the metadata fields 
in "sample_size_3000.txt". The value of the seed used to make the sample was 42.

```R
file <- file.path('..', 'manifest_file_illumina_sra_human')

com <- paste0("wc -l ", file, " | awk '{ print $1 }'")
n <- system(com, intern = TRUE)
n_size = 3000

file.create("relationship_manifest_file-sample")

###

set.seed(42)
line <- sample(1:n, n_size)
x <- rep("NA", n_size)

for(k in 1:length(line)){
    x[k] <- system(paste0("sed -n -e'", line[k], "p' ", file),
                   intern = TRUE)
    write(line[k], "relationship_manifest_file-sample", append=TRUE)
}

write.table(x, file = paste0("sample_size_", n_size ,".txt"), quote = FALSE,
            sep = "\t", row.names = FALSE, col.names = FALSE)

```


Output files:
- __relationship_manifest_file-sample__
- __sample_size_3000.txt__

# Reproducibility

```R
> # Ensure reproducibility
> options(width = 120)
> devtools::session_info()
Session info -----------------------------------------------------------------------------------------------------------
 setting  value                                      
 version  R version 3.1.1 Patched (2014-10-16 r66782)
 system   x86_64, linux-gnu                          
 ui       X11                                        
 language (EN)                                       
 collate  en_US.UTF-8                                
 tz       <NA>                                       

Packages ---------------------------------------------------------------------------------------------------------------
 package    * version date       source                           
 DBI          0.3.1   2014-09-24 CRAN (R 3.1.1)                   
 devtools   * 1.7.0   2015-01-17 CRAN (R 3.1.1)                   
 magrittr     1.5     2015-04-14 Github (smbache/magrittr@89f143d)
 RSQLite      1.0.0   2014-10-25 CRAN (R 3.1.1)                   
 rstudioapi * 0.2     2014-12-31 CRAN (R 3.1.1)  

```
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
