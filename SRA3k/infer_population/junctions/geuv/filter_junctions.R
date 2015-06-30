library('GenomicRanges')

load("junctions_countmatrix_geuvadis.rda")

countsMf <- as.data.frame(count_m$countDF)

means <- Reduce("+", countsMf) / ncol(countsMf)
countsMF <- countsMf[means > 0.95,] 

save(countsMF, file = "junctions_countmatrix_geuvadis_filter0.95.rda")
