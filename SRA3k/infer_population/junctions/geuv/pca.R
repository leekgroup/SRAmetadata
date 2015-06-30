library('GenomicRanges')
library('limma')
library('magrittr')
library('stringr')
library('Cairo')


load("junctions_countmatrix_geuvadis_filter0.6.rda")

print(dim(countsMF))

pops <- countsMF %>% colnames() %>% str_split("_") %>% lapply(`[`, 1) %>% 
		unlist()


# Create model ------------------------------------------------------------

m <- model.matrix(~pops)
fit <- lmFit(countsMF,m)

m0 <- model.matrix(~1, data = as.data.frame(t(countsMF)))
fit0 <- lmFit(countsMF, m0) 


# Get the f statistic from 2 lmFit objects --------------------------------

getF <- function(fit, fit0, theData) {
    
    rss1 <- rowSums((fitted(fit)-theData)^2)
    df1 <- ncol(fit$coef)
    rss0 <- rowSums((fitted(fit0)-theData)^2)
    df0 <- ncol(fit0$coef)
    
    fstat <- ((rss0-rss1)/(df1-df0))/(rss1/(ncol(theData)-df1))
    f_pval <- pf(fstat, df1-1, ncol(theData)-df1,lower.tail=FALSE)
    fout <- cbind(fstat,df1-1,ncol(theData)-df1,f_pval)
    colnames(fout)[2:3] = c("df1","df0")
    fout <- data.frame(fout)
    return(fout)
}

# Get F-statistic 
x <- getF(fit, fit0, countsMF)

# Filter junctions by p-value
countsMF_filter <- x[x$f_pval < 1e-5, ]
countsMF_pass <- countsMF_filter %>% rownames()
countsMF_pass <- countsMF[rownames(countsMF) %in% countsMF_pass,]

head(countsMF_pass)

# Normalize data
countsMf <- log2(countsMF_pass+1)

print(dim(countsMF))


# Get PCA -----------------------------------------------------------------

pca <- countsMf %>% t() %>%  prcomp( center = TRUE,
                            scale. = TRUE)

save(pca, file = "pca_junctions_geuv_filter0.6_ref_normalized.rda")

# load("pca_junctions_geuv_filter0.6_ref_normalized.rda")

pops <- row.names(pca$x) %>% 
    str_split("_") %>% 
    lapply(`[`, 1) %>% 
    unlist() %>% factor()


# Plot PCA ----------------------------------------------------------------

library(ggbiplot)

CairoPNG(file = "pca_junctions_geuv_filter0.6_ref_12.png", 800, 800)
g <- ggbiplot(pca, choices = c(1,2),  obs.scale = 1, var.scale = 1, 
              groups = pops, ellipse = TRUE, 
              circle = FALSE, var.axes = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)
dev.off()

CairoPNG(file = "pca_junctions_geuv_filter0.6_ref_13.png", 800, 800)
g <- ggbiplot(pca, choices = c(1,3),  obs.scale = 1, var.scale = 1, 
              groups = pops, ellipse = TRUE, 
              circle = TRUE, var.axes = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)
dev.off()

CairoPNG(file = "pca_junctions_geuv_filter0.6_ref_23.png", 800, 800)
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

CairoPNG(file = "pca_junctions_geuv_filter0.6_ref_12_base.png", 800, 800)
plot(pca$x[,1], pca$x[,2], 
     col = pops, 
     pch = 19,
     xlab = paste0(names(exp_var[1]), " (",exp_var[1], "% explained var.)"),
     ylab = paste0(names(exp_var[2]), " (",exp_var[2], "% explained var.)")
     )
legend("topleft", levels(pops), col = palette()[1:5], pch = 19, pt.cex = 1,
       title = "populations")
dev.off()

CairoPNG(file = "pca_junctions_geuv_filter0.6_ref_13_base.png", 800, 800)
plot(pca$x[,1], pca$x[,3], 
     col = pops,
     xlab = paste0(names(exp_var[1]), " (",exp_var[1], "% explained var.)"),
     ylab = paste0(names(exp_var[3]), " (",exp_var[3], "% explained var.)"), 
     pch = 19)
legend("bottomleft", levels(pops), col = palette()[1:5], pch = 19, pt.cex = 1,
       title = "populations")
dev.off()

CairoPNG(file = "pca_junctions_geuv_filter0.6_ref_23_base.png", 800, 800)
plot(pca$x[,2], pca$x[,3], 
     col = pops,
     xlab = paste0(names(exp_var[2]), " (",exp_var[2], "% explained var.)"),
     ylab = paste0(names(exp_var[3]), " (",exp_var[3], "% explained var.)"),
     pch = 19)
legend("bottomright", levels(pops), col = palette()[1:5], pch = 19, pt.cex = 1,
       title = "populations")
dev.off()
