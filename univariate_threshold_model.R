# threshold model for binary motif data

library(ape)
library(phytools)
library(tidyverse)
library(taxize)
library(dplyr)
library(reshape2)
library(mvMORPH)
library(umap)
library(stringr)

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

# for each motif, plot in which species this motif is observed
par(mar=margin(1.2,1.2,1.2,1.2))
for (i in 2:9){
  plot(tree, show.tip.label = F, main=paste('motif', i-2))
  
  motif_presence = motif_df[,i]
  tiplabels(bg = ifelse(motif_presence==1, 'red', 'blue'), pch=21, cex=2)
  legend('bottomleft', legend=c('present', 'absent'), col=c('red', 'blue'), cex=2, pt.cex=3, pch=20)
}

# function to fit threshBayes for a motif and return fitted object

one_hot_binary <- function(x) {
  # this function is needed since ancThresh doesn't seem to work witth named vectors
  if (is.null(names(x))) stop("x must be a named vector")
  if (!all(x %in% c(0, 1))) stop("x must contain only 0/1 values")
  
  m <- cbind(
    `0` = as.integer(x == 0),
    `1` = as.integer(x == 1)
  )
  
  rownames(m) <- names(x)
  m
}

fit_threshBayes <- function(i, n_samples=500, n_sims=20000, burnin=NULL, plot=FALSE){
  if(is.null(burnin)){
    burnin = 0.2 * n_sims
  }
  
  X <- as.numeric(motif_df[,i+1+1])
  names(X) <- tree$tip.label
  print('assuming order is same in data and tip labels')
  
  X <- one_hot_binary(X)
  print(dim(X))
  fit <- ancThresh(tree, X, types = 'disc', ngen = n_sims, sequence = c('0','1'),
                   control = list(sample = n_samples, plot=plot))
  return(fit)
}


for (i in 0:7){
fit <- fit_threshBayes(i, n_samples=500, n_sims=1000000, plot=FALSE)

par(mar=margin(1.2,1.2,1.2,1.2))
plot(
  tree,
  show.tip.label = FALSE,
  main=paste('Motif', i),
  cex.main=2
)

tip.cols <- ifelse(motif_df[,2+i] == 1, "firebrick", "steelblue")
tiplabels(
  pch = 21,
  bg  = tip.cols,
  cex = 2
)

node.ids <- (Ntip(tree) + 1):(Ntip(tree) + Nnode(tree))
th_ancr <- fit$ace
cex.scale <- 1 + 1 * th_ancr[,2]

light_red = alpha("firebrick", alpha=0.8)
light_blue = alpha('steelblue', alpha=0.8)

nodelabels(pie=fit$ace, piecol=c(light_blue, light_red), cex=0.5)
nodelabels(node=node.ids,pch=21, bg = ifelse(fit$ace[,2]>=0.5, light_red, light_blue), cex=0.8)

legend('bottomleft', legend=c('present', 'absent'), col=c('firebrick', 'steelblue'), cex=2, pt.cex=3, pch=20)
}
