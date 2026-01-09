# PCA test for phylogenetic signal in binary motif data

library(ape)
library(phytools)
library(tidyverse)
library(taxize)
library(dplyr)
library(reshape2)
library(mvMORPH)
library(umap)
library(stringr)

# load data

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


# performing PCA on binary data

motif_pca <- princomp(data.frame(motif_df)[,2:9])
summary(motif_pca)

# first component explain 26% of variance in data

comp_1_data <- motif_pca$scores[,1]
comp_1_data <- setNames(comp_1_data, motif_df$species)

# testing for signal with pagel's lambda


sig_res <- phylosig(tree, comp_1_data, method='K', test=TRUE, nsim=2500)

class(tree)








