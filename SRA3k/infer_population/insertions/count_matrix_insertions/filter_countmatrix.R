# Filter count matrix by 0.2 presence of insertions in all samples

load("insertions_countmatrix_geuvadis.rda")
print(nrow(count_m$countDF))

means <- Reduce("+", count_m$countDF) / ncol(count_m$countDF)
countsMf <- count_m$countDF[means > 0.1,]

print(nrow(countsMf))

save(countsMf, file = "insertions_countmatrix_geuvadis_filter_0.1.rda")
