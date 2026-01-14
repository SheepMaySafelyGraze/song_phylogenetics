# multivariate threshold model

library(ape)
library(phytools)
library(tidyverse)
library(taxize)
library(dplyr)
library(reshape2)
library(umap)
library(stringr)
library(corrplot)


trait_data <- read.csv('data/namefix_trait_data_100_1912.csv')
tree <- read.tree(file='data/namefix_tree_100_1912.tre')

# names have been matched
geiger::name.check(tree, data.names=trait_data$species)



# read results from threshML run
threshml_covmatrix <- read.csv('threshml_fit_2_covariance_matrix.csv', header = FALSE,
                               sep= ',')
threshml_covmatrix <- matrix(apply(threshml_covmatrix,1, as.numeric),byrow=TRUE,
                             8,8)
rownames(threshml_covmatrix) <- 0:7
colnames(threshml_covmatrix) <- 0:7

# colour palette
pal <- colorRampPalette(c("darkblue", "white", "darkred"))(100)


### comparing with acoustic distances between motif centroids in PCA space
# (we expect positive association, i.e. high distance low correlation)

PCA_means <- trait_data[,c(17:(17+37))] %>%
  group_by(gmm_cluster) %>%
  summarise_all(mean)

dist_PCA <- as.matrix(dist(PCA_means[,-1]))

mantel.test(dist_PCA, threshml_covmatrix, nperm=250000, graph=TRUE,
            alternative='two.sided')


# plotting, reordered by hclust

corrplot(threshml_covmatrix, method='color', title='Correlations between liabilities',
         mar=c(4,4,4,4), tl.col='black', tl.srt=45,
         col=colorRampPalette(c("darkblue", "white", "darkred"))(100),
         order='hclust',
         hclust.method='complete',
         cex.main=2.25,
         tl.cex=1.5,
         cl.cex=1.5)

pal <- colorRampPalette(c("darkblue", "white", "darkred"))(100)


d <- as.dist(1 - threshml_covmatrix)
hc_corr <- hclust(d, method = "average")
plot(hc_corr)

d_acoustic <- as.dist(dist_PCA)
hc_acoustic <- hclust(d_acoustic, method='average')
plot(hc_acoustic)

# ordering from hclust
perm = c(2,3,6,0,1,5,4,7)
dist_PCA_perm <- dist_PCA[perm+1, perm+1]


par(mar = c(6,6,6,6))
image(t(dist_PCA_perm[nrow(dist_PCA_perm):1, ]), col = colorRampPalette(c("white", "darkred"))(100),
      axes=FALSE, main='PCA centroid distances',
      cex.main=2)
axis(1, at = (0:7)/7, labels = perm, las = 2)
axis(2, at = (0:7)/7, labels = rev(perm))
box()

# TODO: add continuous colour bar legend
legend(grconvertX(0.5, "device"), grconvertY(1, "device"), 
       c(min(dist_PCA), 'max'), fill = colorRampPalette(c("white", "darkred"))(2), xpd = NA,
       cex=1.7)


