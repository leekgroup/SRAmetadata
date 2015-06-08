library('R.utils')

file <- file.path('..', '..', 'manifest_file_illumina_sra_human')
n <- countLines(file)[1]
n_size = 500


###

initial_sample <- 1:n

for(i in seq_len(floor(length(initial_sample)/n_size))){
    cat("sample: ",i, "\n")
    file.create(paste0("relationship_manifest_file-sample_", i))
    set.seed(27)
    if(i == floor(n/n_size)){
        #cat("i == floor(length(initial_sample)/n_size", "\t", i, "==", floor(length(initial_sample)/n_size))      
        n_size <- length(initial_sample)
        line <- sample(initial_sample, n_size)
        x <- rep("NA", n_size)
    }else{
         line <- sample(initial_sample, n_size)
         x <- rep("NA", n_size)
    }
    for(k in 1:length(line)){
        x[k] <- scan(file, nlines = 1, skip = line[k] - 1, what = character(), sep = '\n', quiet = TRUE)
        write(line[k], paste0("relationship_manifest_file-sample_", i), append=TRUE)
    }
    
    write.table(x, file = paste0("sample_size_", n_size ,"_", i,".txt"), quote = FALSE,
                sep = "\t", row.names = FALSE, col.names = FALSE)
    
    initial_sample <- setdiff(initial_sample, line)
}
# Temporary solution for warning issue in R.utils package.
# See https://github.com/HenrikBengtsson/R.utils/issues/19#event-284604477
quit(save="no")

# Ensure reproducibility
options(width = 120)
devtools::session_info()
