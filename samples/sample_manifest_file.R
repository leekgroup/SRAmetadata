file <- file.path('..', '/manifest_file_illumina_sra_human')

com <- paste0("wc -l ", file, " | awk '{ print $1 }'")
n <- system(com, intern = TRUE)
n_samples = 1
n_size = 3000

file.create("relationship_manifest_file-sample")

###

set.seed(42)
line <- sample(1:n, n_size, replace = TRUE)
x <- rep("NA", n_size)

for(k in 1:length(line)){
    x[k] <- system(paste0("sed -n -e'", line[k], "p' ", file),
                   intern = TRUE)
    write(line[k], "relationship_manifest_file-sample", append=TRUE)
}

write.table(x, file = paste0("sample_size_", n_size ,".txt"), quote = FALSE,
            sep = "\t", row.names = FALSE, col.names = FALSE)
