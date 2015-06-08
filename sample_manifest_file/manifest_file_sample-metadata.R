setwd("/home/joseah/Documents/jeff_leek_lab/github/SRAmetadata/sample_manifest_file")

library('magrittr')
library('stringr')
md <- read.table('../all_illumina_sra_for_human.txt', 
                 sep = "\t",
                 header = TRUE, 
                 stringsAsFactors = FALSE, 
                 quote = "",
                 comment.char = "")


mf <- read.table('sample_size_3000.txt', 
                 sep = "\t",
                 header = FALSE, 
                 stringsAsFactors = FALSE, 
                 quote = "",
                 comment.char = "",
                 colClasses = c("character", rep("NULL", 4)), fill = TRUE)

runacc <- as.vector(mf$V1) %>% str_split("\\.") %>% lapply(`[[`, 5) %>% 
    unlist() %>% str_split("/") %>% lapply(`[[`, 6) %>% unlist() %>% 
    str_replace("_1", "")

md_sample <- md[match(runacc, md$run_accession),]

write.table(md_sample, "sample_size_3000_metadata.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)