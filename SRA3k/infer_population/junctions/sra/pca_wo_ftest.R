library('GenomicRanges')
library('limma')
library('magrittr')
library('stringr')
library('Cairo')

"%p%"  <- function(x, y) paste0(x, y)

n <- 95 

load("junctions_countmatrix_sra_filter0." %p%  n %p% ".rda")

print(dim(countsMF))


# Normalize data
countsMf <- log2(countsMF+1)


# Get PCA -----------------------------------------------------------------

pca <- countsMf %>% t() %>%  prcomp( center = TRUE,
                            scale. = TRUE)

save(pca, file = "pca_junctions_geuv_filter0." %p% n %p% " _wo_ftest.rda")


