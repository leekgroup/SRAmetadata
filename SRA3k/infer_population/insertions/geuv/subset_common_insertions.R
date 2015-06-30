library('magrittr')
library('stringr')
library('GenomicRanges')

load("insertions_countmatrix_geuvadis.rda")

ref <- read.delim("../ALL_wgs_phase3_indels_insertions.vcf", 
           header = FALSE, 
           col.names = c("chr", "start", "id"))

print(head(count_m$countDF))


chr_start_counts <- rownames(count_m$countDF) %>% str_split("-") %>% lapply(`[`, 1) %>% 
    unlist() %>% str_split(":")

print(head(chr_start_counts))


row.names(count_m) %>% head() %>% print()
row.names(count_m) %>% str_split("-") %>% head() %>% print()
row.names(count_m) %>% str_split("-") %>% lapply(`[`, 1) %>% head() %>% print()
row.names(count_m) %>% str_split("-") %>% lapply(`[`, 1) %>%
    unlist() %>% head() %>% print()
row.names(count_m) %>% str_split("-") %>% lapply(`[`, 1) %>%
    unlist() %>% str_split(":") %>% head() %>% print()


chr_start_counts <- data.frame(chr = lapply(chr_start_counts, `[`, 1) %>% unlist(), 
                        start = lapply(chr_start_counts, `[`, 2) %>% unlist())

print(head(chr_start_counts))

ref$chr <- paste0('chr', ref$chr)

select_dels_ref <- (chr_start_counts$chr %in% ref$chr) & (chr_start_counts$start %in% ref$start)
print(head(select_dels_ref))

countsMf_ref <- count_m$countDF[select_dels_ref, ]

save(countsMf_ref, file = "insertions_countmatrix_geuvadis_ref.rda")

#count_m
