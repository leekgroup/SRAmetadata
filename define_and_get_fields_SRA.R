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
query <- "SELECT exp.experiment_accession, exp.study_accession, exp.sample_accession, ru.run_accession,
sam.broker_name 'sam.broker_nam',
sam.center_name 'sam.center_name',
sam.taxon_id,
sam.scientific_name,
sam.common_name,
sam.anonymized_name,
sam.individual_name,
sam.description 'sam.description',
sam.sample_attribute,
sam.sradb_updated 'sam.sradb_updated',
stu.study_title,
stu.study_type,
stu.study_abstract,
stu.broker_name 'stu.broker_name',
stu.center_name 'stu.center_name',
stu.center_project_name,
stu.study_description,
stu.related_studies,
stu.primary_study,
stu.study_attribute,
stu.sradb_updated 'stu.sradb_updated',
ru.bamFile 'ru.bamFile',
ru.broker_name 'ru.broker_name',
ru.instrument_name,
ru.run_date,
ru.run_file,
ru.run_center,
ru.total_data_blocks,
ru.experiment_name,
ru.run_attribute,
ru.sradb_updated 'ru.sradb_updated',
exp.bamFile 'exp.bamFile',
exp.fastqFTP,
exp.broker_name 'exp.broker_name',
exp.center_name 'exp.center_name',
exp.title,
exp.study_name,
exp.design_description,
exp.sample_name,
exp.sample_member,
exp.library_name,
exp.library_strategy,
exp.library_source,
exp.library_selection,
exp.library_layout,
exp.targeted_loci,
exp.library_construction_protocol,
exp.spot_length,
exp.adapter_spec,
exp.read_spec,
exp.platform,
exp.instrument_model,
exp.platform_parameters,
exp.sequence_space,
exp.base_caller,
exp.quality_scorer,
exp.number_of_levels,
exp.multiplier,
exp.qtype,
exp.experiment_attribute,
exp.submission_accession,
exp.sradb_updated 'exp.sradb_updated'
FROM sample sam, experiment exp, study stu, run ru
WHERE sam.sample_accession = exp.sample_accession AND 
exp.study_accession= stu.study_accession AND 
exp.experiment_accession = ru.experiment_accession AND
exp.instrument_model LIKE '%illumina%'"
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
metadata$individual <- search.field(metadata$sample_attribute, "individual")

# Get sex of individuals
metadata$sex <- search.field(metadata$sample_attribute, "sex")

# Get source type of samples
metadata$strain <- search.field(metadata$sample_attribute, "Strain")
metadata$source_name <- search.field(metadata$sample_attribute, "source_name")

# Write table with all Illimina data
write.table(metadata, "all_illumina_sra.txt", sep = "\t", quote = FALSE,
    row.names = FALSE)

# Create new database
meta_con <- dbConnect(SQLite(), dbname = "test.sqlite")

# Define parameters to create table for metadata
table_name = "metadata"
primary_key_name = "meta_id"
fields <- sub("\\.", "__", colnames(metadata))

ids <- paste(fields[1:4], " varchar(255) not NULL,", sep = "", collapse = " ")
registers <- paste(fields[5:length(fields)], "longtext default 'NA'", 
    sep = " ", collapse = ", ")

query <- "CREATE TABLE " %p% table_name %p% "(" %p% 
    primary_key_name %p% " int unsigned auto_increment not NULL," %p%
    ids %p% registers %p% ");"

dbGetQuery(meta_con, query)

for(i in 1:nrow(metadata)) {
    x <- unlist(gsub("'", "prime", metadata[i,], fixed = TRUE))
    values <- paste("'" %p% x, "'", sep = "", collapse = ",")
    query <- "INSERT INTO " %p%  table_name %p% " VALUES (" %p% '"' %p%
        i %p% '"' %p% "," %p% values %p% ");"
    dbGetQuery(meta_con, query)
}

dbDisconnect(sra_con)
dbDisconnect(meta_con)
