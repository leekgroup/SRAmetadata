library('GenomicRanges')
library('limma')
library('magrittr')
library('stringr')
library('Cairo')

"%p%"  <- function(x, y) paste0(x, y)

n <- 95 

load("junctions_countmatrix_geuvadis_filter0." %p%  n %p% ".rda")

print(dim(countsMF))

pops <- countsMF %>% colnames() %>% str_split("_") %>% lapply(`[`, 1) %>% 
		unlist()
center <- countsMF %>% colnames() %>% str_split("_") %>% lapply(`[`, 3) %>%
                unlist()

# Normalize data
countsMf <- log2(countsMF+1)



# Get PCA -----------------------------------------------------------------

pca <- countsMf %>% t() %>%  prcomp( center = TRUE,
                            scale. = TRUE)

save(pca, file = "pca_junctions_geuv_filter0." %p% n %p% " _wo_ftest.rda")

# load("pca_junctions_geuv_filter0." %p% n %p% " _wo_ftest.rda")

pops <- row.names(pca$x) %>% 
    str_split("_") %>% 
    lapply(`[`, 1) %>% 
    unlist() %>% factor()

center <- row.names(pca$x) %>% 
    str_split("_") %>% 
    lapply(`[`, 3) %>% 
    unlist() %>% factor()


# Plot PCA ----------------------------------------------------------------

library(ggbiplot)

CairoPNG(file = "pca_junctions_geuv_filter0." %p% n %p% "_woftest_12.png", 800, 800)
g <- ggbiplot(pca, choices = c(1,2),  obs.scale = 1, var.scale = 1, 
              groups = pops, ellipse = TRUE, 
              circle = FALSE, var.axes = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)
dev.off()

CairoPNG(file = "pca_junctions_geuv_filter0." %p% n %p% "_woftest_13.png", 800, 800)
g <- ggbiplot(pca, choices = c(1,3),  obs.scale = 1, var.scale = 1, 
              groups = pops, ellipse = TRUE, 
              circle = TRUE, var.axes = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)
dev.off()

CairoPNG(file = "pca_junctions_geuv_filter0." %p% n %p% "_woftest_23.png", 800, 800)
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

CairoPNG(file = "pca_junctions_geuv_filter0." %p% n %p% "_woftest_12_base.png", 800, 800)
plot(pca$x[,1], pca$x[,2], 
     col = pops, 
     pch = 19,
     xlab = paste0(names(exp_var[1]), " (",exp_var[1], "% explained var.)"),
     ylab = paste0(names(exp_var[2]), " (",exp_var[2], "% explained var.)")
     )
legend("topleft", levels(pops), col = palette()[1:5], pch = 19, pt.cex = 1,
       title = "populations")
dev.off()

CairoPNG(file = "pca_junctions_geuv_filter0." %p% n %p% "_woftest_13_base.png", 800, 800)
plot(pca$x[,1], pca$x[,3], 
     col = pops,
     xlab = paste0(names(exp_var[1]), " (",exp_var[1], "% explained var.)"),
     ylab = paste0(names(exp_var[3]), " (",exp_var[3], "% explained var.)"), 
     pch = 19)
legend("bottomleft", levels(pops), col = palette()[1:5], pch = 19, pt.cex = 1,
       title = "populations")
dev.off()

CairoPNG(file = "pca_junctions_geuv_filter0." %p% n %p% "_woftest_23_base.png", 800, 800)
plot(pca$x[,2], pca$x[,3], 
     col = pops,
     xlab = paste0(names(exp_var[2]), " (",exp_var[2], "% explained var.)"),
     ylab = paste0(names(exp_var[3]), " (",exp_var[3], "% explained var.)"),
     pch = 19)
legend("bottomright", levels(pops), col = palette()[1:5], pch = 19, pt.cex = 1,
       title = "populations")
dev.off()


# center


CairoPNG(file = "pca_junctions_geuv_filter0." %p% n %p% "_woftest_12_center.png", 800, 800)
g <- ggbiplot(pca, choices = c(1,2),  obs.scale = 1, var.scale = 1,
              groups = center, ellipse = TRUE,
              circle = FALSE, var.axes = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal',
               legend.position = 'top')
print(g)
dev.off()

CairoPNG(file = "pca_junctions_geuv_filter0." %p% n %p%"_woftest_13_center.png", 800, 800)
g <- ggbiplot(pca, choices = c(1,3),  obs.scale = 1, var.scale = 1,
              groups = center, ellipse = TRUE,
              circle = TRUE, var.axes = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal',
               legend.position = 'top')
print(g)
dev.off()

CairoPNG(file = "pca_junctions_geuv_filter0." %p% n %p% "_woftest_23_center.png", 800, 800)
g <- ggbiplot(pca, choices = c(2,3),  obs.scale = 1, var.scale = 1,
              groups = center, ellipse = TRUE,
              circle = TRUE, var.axes = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal',
               legend.position = 'top')
print(g)
dev.off()

