# scratch file testing issues with ancestral nodes under threshold model

library(ape)
library(phytools)
library(tidyverse)
library(taxize)
library(dplyr)
library(reshape2)
library(stringr)
library(adephylo)
library(mvtnorm)

trait_data <- read.csv('phylogenetic_analysis/namefix_trait_data_100_1912.csv')
tree <- read.tree(file='phylogenetic_analysis/namefix_tree_100_1912.tre')

# names have been matched
geiger::name.check(tree, data.names=trait_data$species)

# writing out to pass to birdtree.org
write.csv(data.frame(as.character(lapply(tree$tip.label, \(x) str_replace(x, '_', ' ')))),
          row.names=FALSE, file='birds.csv')



for (i in 0:7){
  motif_i_present = (trait_data$gmm_cluster == i)
  trait_data[paste('motif_', i,'_present', sep='')] = motif_i_present
}

motif_data <- trait_data %>%
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

motif_data <- data.frame(motif_data)
species <- tree$tip.label

perm <- match(species, motif_data$species)

motif_data <- motif_data[perm,]
motif_data <- motif_data[,2:9, drop=TRUE]

# convert to 1/0
motif_data <- data.frame(lapply(motif_data, as.integer))
rownames(motif_data) <- species


# read birdtree results (sampled trees from Ericson et al.) and select tree
birdtree_trees <- read.nexus(file='birdtree_trees.nex')

ind = 37

bt_tree <- birdtree_trees[[ind]]

bt_thresh_fit_near <- ancThresh(tree, setNames(motif_data[,8], species), control=list(sample=10000, burnin=5e5, propliab=1.5),
                           ngen=5e6, quiet=TRUE)




# fitting and appending over 8 chains, 1.5 million steps each
ngen =3e6
n_chains = 16


for (i in 0:7){
  print(i)
  
  fit <- ancThresh(tree, setNames(motif_data[,i+1], species),
                   ngen=ngen, control=list(print=FALSE, sample=12000))
  par <- fit$par
  liab <- fit$liab
  ace <- fit$ace
  for (j in 2:n_chains){
    fit <- ancThresh(tree, setNames(motif_data[,i+1], species),
                     ngen=ngen, control=list(print=FALSE, sample=12000))
    ace <- ace + fit$ace
    par <- rbind(par, fit$par)
    liab <- rbind(liab, fit$liab)
  }
  print(ace/n_chains)
  write.csv(ace/n_chains, file=paste('results/motif_', i,'_ace.cex', sep=''))
  write.csv(par, file=paste('results/motif_', i,'_par.cex', sep=''))
  write.csv(liab, file=paste('results/motif_', i,'_liab.cex', sep=''))
}




