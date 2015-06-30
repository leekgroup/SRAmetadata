library('GenomicRanges')
library('magrittr')
library('stringr')
library('Cairo')


#load('./insertions/insertions_countmatrix_geuvadis_ref.rda')
#countsMf_ref <- as.data.frame(countsMf_ref)
#means <- Reduce("+", countsMf_ref) / ncol(countsMf_ref)
#insertions <- countsMf_ref[means > 0.1,] 


#load('./deletions/deletions_countmatrix_geuvadis_ref.rda')
#countsMf_ref <- as.data.frame(countsMf_ref)
#means <- Reduce("+", countsMf_ref) / ncol(countsMf_ref)
#deletions <- countsMf_ref[means > 0.1,]

#indels <- rbind(insertions, deletions)
#save(indels, file = "indels_geuvadis_ref_filter_0.1.rda")

#load("indels_geuvadis_ref.rda")
#load("indels_geuvadis_ref_filter_0.1.rda")


#print(dim(indels))

#countsMf <- log2(indels+1)

#pca <- countsMf %>% t() %>%  prcomp( center = TRUE,
#                            scale. = TRUE)

#save(pca, file = "pca_indels_ref_normalized.rda")

load("pca_indels_ref_normalized.rda")

pops <- row.names(pca$x) %>% 
    str_split("_") %>% 
    lapply(`[`, 1) %>% 
    unlist() %>% factor()


library(ggbiplot)

CairoPNG(file = "pca_indels_ref_12.png", 800, 800)
g <- ggbiplot(pca, choices = c(1,2),  obs.scale = 1, var.scale = 1, 
              groups = pops, ellipse = TRUE, 
              circle = FALSE, var.axes = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)
dev.off()

CairoPNG(file = "pca_indels_ref_13.png", 800, 800)
g <- ggbiplot(pca, choices = c(1,3),  obs.scale = 1, var.scale = 1, 
              groups = pops, ellipse = TRUE, 
              circle = TRUE, var.axes = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)
dev.off()

CairoPNG(file = "pca_indels_ref_23.png", 800, 800)
g <- ggbiplot(pca, choices = c(2,3),  obs.scale = 1, var.scale = 1, 
              groups = pops, ellipse = TRUE, 
              circle = TRUE, var.axes = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)
dev.off()


exp_var <- round(summary(pca)$importance[2,] * 100, 1)
paste0(names(exp_var[1]), " (",exp_var[1], "% explained var.)")

palette(c("firebrick1", "seagreen3", "yellowgreen", "turquoise3", "violetred"))

CairoPNG(file = "pca_indels_ref_12_base.png", 800, 800)
plot(pca$x[,1], pca$x[,2], 
     col = pops, 
     pch = 19,
     xlab = paste0(names(exp_var[1]), " (",exp_var[1], "% explained var.)"),
     ylab = paste0(names(exp_var[2]), " (",exp_var[2], "% explained var.)")
     )
legend("topleft", levels(pops), col = palette()[1:5], pch = 19, pt.cex = 1,
       title = "populations")
dev.off()

CairoPNG(file = "pca_indels_ref_13_base.png", 800, 800)
plot(pca$x[,1], pca$x[,3], 
     col = pops,
     xlab = paste0(names(exp_var[1]), " (",exp_var[1], "% explained var.)"),
     ylab = paste0(names(exp_var[3]), " (",exp_var[3], "% explained var.)"), 
     pch = 19)
legend("bottomleft", levels(pops), col = palette()[1:5], pch = 19, pt.cex = 1,
       title = "populations")
dev.off()

CairoPNG(file = "pca_indels_ref_23_base.png", 800, 800)
plot(pca$x[,2], pca$x[,3], 
     col = pops,
     xlab = paste0(names(exp_var[2]), " (",exp_var[2], "% explained var.)"),
     ylab = paste0(names(exp_var[3]), " (",exp_var[3], "% explained var.)"),
     pch = 19)
legend("bottomright", levels(pops), col = palette()[1:5], pch = 19, pt.cex = 1,
       title = "populations")
dev.off()
