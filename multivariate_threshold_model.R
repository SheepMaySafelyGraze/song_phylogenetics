# multivariate threshold model

library(ape)
library(phytools)
library(tidyverse)
library(taxize)
library(dplyr)
library(reshape2)
library(mvMORPH)
library(umap)
library(stringr)
library(corrplot)

trait_data <- read.csv('phylogenetic_analysis/namefix_trait_data_100_1912.csv')
tree <- read.tree(file='phylogenetic_analysis/namefix_tree_100_1912.tre')

# names have been matched
geiger::name.check(tree, data.names=trait_data$species)

# adding categorical variables: motif i present in each species

for (i in 0:7){
  motif_i_present = (trait_data$gmm_cluster == i)
  trait_data[paste('motif_', i,'_present', sep='')] = motif_i_present
}

motif_df <- trait_data %>%
  group_by(species) %>%
  summarise(motif_0 = any(motif_0_present),
            motif_1 = any(motif_1_present),
            motif_2 = any(motif_2_present),
            motif_3 = any(motif_3_present),
            motif_4 = any(motif_4_present),
            motif_5 = any(motif_5_present),
            motif_6 = any(motif_6_present),
            motif_7 = any(motif_7_present))



# sort motif_df so that entries are in the same order as tree labels

motif_df <- data.frame(motif_df)
species <- tree$tip.label

rownames(motif_df) <- species

perm <- match(tree$tip.label, motif_df$species)

motif_df <- motif_df[perm,]

# convert to 1/0
motif_df <- lapply(motif_df, as.integer)
motif_df <- data.frame(motif_df)

### function to apply multivariate threshold model via MCMC


library(Rphylip)
# path to threshml exe file
path = r"{C:\Users\lfd20\Documents\Code\Proj1-Code\MotifAnalysis\phylip-3.698\exe}"

# casting motif_df to named array
motif_arr <- as.matrix(motif_df[,2:9])
rownames(motif_arr) <- species

# WARNING: very long runtime
fit_8dim <- Rthreshml(tree, motif_arr, types=rep('discrete',8), path=path, cleanup=FALSE)

print(fit_8dim$Covariance_matrix)
# the covariance matrix seems to be saved incorrectly, though outfile correct

# manually reading in
threshml_covmatrix <- read.csv('threshml_fit_1_covariance_matrix.csv', header = FALSE,
                               sep= ',')
threshml_covmatrix <- matrix(apply(threshml_covmatrix,1, as.numeric),byrow=TRUE,
                             8,8)
rownames(threshml_covmatrix) <- 0:7
colnames(threshml_covmatrix) <- 0:7

corrplot(threshml_covmatrix, method='color', title='Correlations between liabilities',
mar=c(2,2,2,2), tl.col='black', tl.srt=45,
col=colorRampPalette(c("darkblue", "white", "darkred"))(100))

pal <- colorRampPalette(c("darkblue", "white", "darkred"))(100)

# comparing with acoustic distances between PCA comps (we expect positive association, i.e. high distance low correlation)

PCA_means <- trait_data[,c(17:(17+37))] %>%
  group_by(gmm_cluster) %>%
  summarise_all(mean)

dist_PCA <- as.matrix(dist(PCA_means[,-1]))

mantel.test(dist_PCA, threshml_covmatrix, nperm=250000, graph=TRUE,
            alternative='two.sided')

par(mar = c(6,6,6,6))
image(t(dist_PCA[nrow(dist_PCA):1, ]), col = colorRampPalette(c("darkblue", "white", "darkred"))(100),
      axes=FALSE, main='PCA centroid distances')
axis(1, at = (0:7)/7, labels = 0:7, las = 2)   # column labels
axis(2, at = (0:7)/7, labels = rev(0:7))       # row labels
box()


legend(grconvertX(0.5, "device"), grconvertY(1, "device"), 
       c(min(dist_PCA), round(max(dist_PCA),1)), fill = pal[c(1,100)], xpd = NA)
