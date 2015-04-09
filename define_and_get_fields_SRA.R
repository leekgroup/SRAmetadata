# Load library
library('RSQLite')

# Define functions
"%p%"  <- function(x, y) paste0(x, y)

# Update sql if necessary:
sra.meta <- file.path('..', 'SRAmetadb.sqlite') # Using the file from another location to avoid having copies of this 2.4 Gb file
if(file.exists(sra.meta)) {
  sqlfile <- sra.meta
} else {
  sqlfile <- getSRAdbFile()
}

# Create connection
sra_con <- dbConnect(SQLite(),sqlfile)

# Query database
timeStart <- proc.time()

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
taxon_id = 9606"
selected <- dbGetQuery(sra_con, query)
print("Consumed time of query:")
proc.time() - timeStart

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

# Getting URL to fastq files
metadata$URL <- get.fastq.urls(metadata)

search.field <- function(column, field) {
  r <- paste0('.*\\|?\\|?', field,': (.*?) \\|?\\|.*|.*\\|?\\|?', field,
              ': (.*?)\\|?\\|?.*?')
  res <- sub(r, "\\1", column, perl = TRUE, ignore.case = TRUE)
  unlist(lapply(res, function(x) 
    if(grepl(".*\\|\\|.*", x) == TRUE) {
      return(NA)
    } else {
      return(x)
    })
  )
}

# Get population
metadata$population <- search.field(metadata$sample_attribute, "population")

# Get cell line / ID for individuals
metadata$cell_line <- search.field(metadata$sample_attribute, "cell line")

# Get sex of individuals
metadata$sex <- search.field(metadata$sample_attribute, "sex")

# Get source type of samples
metadata$strain <- search.field(metadata$sample_attribute, "Strain")
metadata$source_name <- search.field(metadata$sample_attribute, "source_name")

# Removing new line characters
metadata$read_spec <- gsub("\n", "", metadata$read_spec)

# Write table with all Illimina data
write.table(metadata, "all_illumina_sra_for_human.txt", sep = "\t", quote = FALSE,
            row.names = FALSE)

# # Create new database
# meta_con <- dbConnect(SQLite(), dbname = "test.sqlite")
# 
# # Define parameters to create table for metadata
# table_name = "metadata"
# primary_key_name = "meta_id"
# fields <- sub("\\.", "__", colnames(metadata))
# 
# ids <- paste(fields[1:4], " varchar(255) not NULL,", sep = "", collapse = " ")
# registers <- paste(fields[5:length(fields)], "longtext default 'NA'", 
#                    sep = " ", collapse = ", ")
# 
# query <- "CREATE TABLE " %p% table_name %p% "(" %p% 
#   primary_key_name %p% " int unsigned auto_increment not NULL," %p%
#   ids %p% registers %p% ");"
# 
# dbGetQuery(meta_con, query)
# 
# for(i in 1:nrow(metadata)) {
#   x <- unlist(gsub("'", "prime", metadata[i,], fixed = TRUE))
#   values <- paste("'" %p% x, "'", sep = "", collapse = ",")
#   query <- "INSERT INTO " %p%  table_name %p% " VALUES (" %p% '"' %p%
#     i %p% '"' %p% "," %p% values %p% ");"
#   dbGetQuery(meta_con, query)
# }

dbDisconnect(sra_con)
# dbDisconnect(meta_con)
