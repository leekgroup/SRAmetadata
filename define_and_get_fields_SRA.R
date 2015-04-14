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
experiment_accession, study_accession, sample_accession, 
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
experiment_attribute,
submission_accession
FROM sra
WHERE platform = 'ILLUMINA' AND
study_type = 'Transcriptome Analysis' AND
taxon_id = 9606;"

selected <- dbGetQuery(sra_con, query)
query <- paste0("SELECT run_accession, file_name, md5, bytes,sradb_updated FROM fastq WHERE run_accession IN ('", 
                paste(selected$run_accession, collapse="', '"), "')")
fastq_namefiles <- dbGetQuery(sra_con, query)
metadata <- merge(selected, fastq_namefiles, by = "run_accession")

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

# Write table with all Illimina data
write.table(metadata, "all_illumina_sra_for_human.txt", sep = "\t", quote = FALSE,
            row.names = FALSE)

# Disconnect from db
dbDisconnect(sra_con)

