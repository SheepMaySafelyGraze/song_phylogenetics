# PCA reconstruction in clades predicted to have a single motif 'birth'

library(ape)
library(phytools)
library(tidyverse)
library(dplyr)
library(stringr)
library(mvMORPH)
library(umap)

trait_data <- read.csv('phylogenetic_analysis/namefix_trait_data_100_1912.csv')
tree <- read.tree(file='phylogenetic_analysis/namefix_tree_100_1912.tre')

species <- tree$tip.label

PCA_data <- trait_data[,c(4,54,17:(17+36))]

# aggregating PCA values
agg_PCA <- PCA_data %>%
  group_by(species, gmm_cluster) %>%
  summarise(across(everything(), mean), .groups = 'drop')


####### functions to fit BM model and reconstruct ancestral values
reconstruct_PCA <- function(data, tree, method='PL-LOOCV', model='BM'){
  
  dat_list <- list(Y=as.matrix(data))
  
  fit <- mvgls(Y~1, data = dat_list, tree = tree, method=method, model=model)
  
  ancr <- ancestral(fit)
  
  return(list(ancr, fit))
}

get_ancestral_reconstruction <- function(motif_ind, data, tree, n_comp=2,
                                         method='PL-LOOCV', model='BM'){
  # expects aggregated PCA data with gmm_cluster and species columns
  
  # subset to motif
  agg_PCA_i <- data[data$gmm_cluster==motif_ind,]
  rownames(agg_PCA_i) <- agg_PCA_i$species
  
  # drop taxa for which motif does not occur
  drop_taxa <- setdiff(tree$tip.label, agg_PCA_i$species)
  tree_i <- drop.tip(tree, drop_taxa)
  
  # ensure data ordered like tree tips
  perm <- match(tree_i$tip.label, agg_PCA_i$species)
  agg_PCA_i <- agg_PCA_i[perm,]
  
  # drop to number of components, drop species and cluster
  agg_PCA_i <- agg_PCA_i[,3:(3 + n_comp - 1)]
  
  res <- reconstruct_PCA(agg_PCA_i, tree_i, method=method, model='BM')
  
  return(list(agg_PCA_i ,res[[1]], res[[2]], tree_i))
}


####### helper functions for plotting 
root_distances <- function(tree) {
  # Number of tips and nodes
  n_tips <- length(tree$tip.label)
  n_nodes <- tree$Nnode
  total_nodes <- n_tips + n_nodes
  
  # Initialize distance vector
  dist <- numeric(total_nodes)
  
  # Postorder traversal ensures parents are processed before children
  edge <- tree$edge
  edge.length <- tree$edge.length
  
  for (i in seq_len(nrow(edge))) {
    parent <- edge[i, 1]
    child  <- edge[i, 2]
    dist[child] <- dist[parent] + edge.length[i]
  }
  
  # Name the vector
  names(dist) <- c(tree$tip.label,
                   paste0("Node", seq_len(n_nodes)))
  
  return(dist)
}


plot_PCA <- function(trait_PC, anc_PC, pruned_tree, umap=FALSE, scale=FALSE){
  n = length(pruned_tree$tip.label) # number of tips
  M = length(pruned_tree$node.label) + n # total nodes
  
  dists_to_root <- root_distances(pruned_tree)
  dists_to_root <- dists_to_root/max(dists_to_root)
  anc_dists = dists_to_root[(n+1):M]
  
  dat <- rbind(trait_PC, anc_PC)
  if (scale){
    dat
  }
  
  if (umap | dim(anc_PC)[2] != 2) {
    print('using umap')
    trait.umap <- umap(rbind(trait_PC, anc_PC))
    coords = trait.umap$layout
  } else {
  coords = rbind(trait_PC, anc_PC)
  }
  
  
  plot(coords[1:n,],
       pch = 21,
       bg = rgb(red=0, green=0, blue=1, alpha=0.5),
       asp = 1,
       xlab='',
       ylab='',
       cex=3.5,
       xlim=c(min(coords[,1])-0.3, max(coords[,1])+0.3),
       ylim=c(min(coords[,2])-0.3, max(coords[,2])+0.3),
       cex.axis=1.5)
  
  points(coords[(n+1):M,],
         pch = 21,
         bg = rgb(red=1, green=0, blue=0, alpha=0.5),
         asp = 1,
         xlab='',
         ylab='',
         cex=2.5*(anc_dists) + 1)
  
  
  for(i in 1:nrow(pruned_tree$edge)) {
    parent <- pruned_tree$edge[i, 1]
    child  <- pruned_tree$edge[i, 2]
    
    segments(coords[parent, 1], coords[parent, 2],
             coords[child, 1],  coords[child, 2],
             col = rgb(red=0, green=0, blue=0, alpha=0.1))
  }
  
  if (umap) {
    return(trait.umap)
  }
}


get_living_fossil <- function(taxa_df, root_val){
  distances <- apply(taxa_df, 1, function(row) {
    sqrt(sum((row - root_val)^2))
  })
  
  return(which.min(distances))
}





# per motif presence reconstruction figures
### for motif 7, clade 110 seems to have only a single birth of the motif

clade_110 <- extract.clade(tree, 110)
species_110 <- clade_110$tip.label
agg_PCA_clade_110 <- data.frame(agg_PCA[agg_PCA$species %in% species_110,])

### for motif 2, clade 186 seems to have only a single birth of the motif

clade_186 <- extract.clade(tree, 186)
species_186 <- clade_186$tip.label
agg_PCA_clade_186 <- data.frame(agg_PCA[agg_PCA$species %in% species_186,])



# motif 7 within clade 110
motif_ind <- 7
max_PC <- 2
res <- get_ancestral_reconstruction(motif_ind, agg_PCA_clade_110, clade_110, n_comp=max_PC)
plot_map <- plot_PCA(data.frame(res[[1]]), res[[2]], res[[4]], umap=FALSE)

n <- nrow(res[[1]])  # number of tips
# identify and highlight current taxa most similar to root taxa
if (is.null(plot_map)){
  points(res[[2]][1,1], res[[2]][1,2],pch=3, col='maroon', cex=3)
  
  
  living_fossil_ind <- get_living_fossil(res[[1]], res[[2]][1,])
  living_fossil <- str_replace(rownames(res[[1]])[living_fossil_ind], '_', ' ')
  
  points(res[[1]][living_fossil_ind,1],res[[1]][living_fossil_ind,2],
         pch=16, col='aquamarine', cex=2)
} else {
  points(plot_map$layout[n+1,1], plot_map$layout[n+1,2],pch=3, col='maroon', cex=3)
  
  living_fossil_ind <- get_living_fossil(res[[1]], res[[2]][1,])
  living_fossil <- str_replace(rownames(res[[1]])[living_fossil_ind], '_', ' ')
  
  points(plot_map$layout[living_fossil_ind,1],plot_map$layout[living_fossil_ind,2],
         pch=16, col='aquamarine', cex=2)
}

legend('topright',pch=c(16,16,3,16), col=c(rgb(red=0, green=0, blue=1, alpha=0.5),rgb(red=1, green=0, blue=0, alpha=0.5),'maroon','aquamarine'),
       legend=c('extant','ancestral','root', living_fossil), cex=1.3)
title(paste('Reconstruction of first ', paste(max_PC, 'PCs'), '; motif ', motif_ind, sep=''), cex.main=2)


# motif 6 also has a single birth before clade 110
motif_ind <- 6
max_PC <- 2
res <- get_ancestral_reconstruction(6, agg_PCA_clade_110, clade_110, n_comp=max_PC)
plot_map <- plot_PCA(data.frame(res[[1]]), res[[2]], res[[4]], umap=FALSE)

n <- nrow(res[[1]])  # number of tips
# identify and highlight current taxa most similar to root taxa
if (is.null(plot_map)){
  points(res[[2]][1,1], res[[2]][1,2],pch=3, col='maroon', cex=3)
  
  
  living_fossil_ind <- get_living_fossil(res[[1]], res[[2]][1,])
  living_fossil <- str_replace(rownames(res[[1]])[living_fossil_ind], '_', ' ')
  
  points(res[[1]][living_fossil_ind,1],res[[1]][living_fossil_ind,2],
         pch=16, col='aquamarine', cex=2)
} else {
  points(plot_map$layout[n+1,1], plot_map$layout[n+1,2],pch=3, col='maroon', cex=3)
  
  living_fossil_ind <- get_living_fossil(res[[1]], res[[2]][1,])
  living_fossil <- str_replace(rownames(res[[1]])[living_fossil_ind], '_', ' ')
  
  points(plot_map$layout[living_fossil_ind,1],plot_map$layout[living_fossil_ind,2],
         pch=16, col='aquamarine', cex=2)
}

legend('topright',pch=c(16,16,3,16), col=c(rgb(red=0, green=0, blue=1, alpha=0.5),rgb(red=1, green=0, blue=0, alpha=0.5),'maroon','aquamarine'),
       legend=c('extant','ancestral','root', living_fossil), cex=1.3)
title(paste('Reconstruction of first ', paste(max_PC, 'PCs'), '; motif ', motif_ind, sep=''), cex.main=2)





# motif 2 within clade 186
max_PC <- 2
motif_ind <- 2
res <- get_ancestral_reconstruction(motif_ind, agg_PCA_clade_186, clade_186, n_comp=max_PC)
plot_map <- plot_PCA(data.frame(res[[1]]), res[[2]], res[[4]], umap=FALSE)

n <- nrow(res[[1]])  # number of tips
# identify and highlight current taxa most similar to root taxa
if (is.null(plot_map)){
  points(res[[2]][1,1], res[[2]][1,2],pch=3, col='maroon', cex=3)
  
  
  living_fossil_ind <- get_living_fossil(res[[1]], res[[2]][1,])
  living_fossil <- str_replace(rownames(res[[1]])[living_fossil_ind], '_', ' ')
  
  points(res[[1]][living_fossil_ind,1],res[[1]][living_fossil_ind,2],
         pch=16, col='aquamarine', cex=2.5)
} else {
  points(plot_map$layout[n+1,1], plot_map$layout[n+1,2],pch=3, col='maroon', cex=3)
  
  living_fossil_ind <- get_living_fossil(res[[1]], res[[2]][1,])
  living_fossil <- str_replace(rownames(res[[1]])[living_fossil_ind], '_', ' ')
  
  points(plot_map$layout[living_fossil_ind,1],plot_map$layout[living_fossil_ind,2],
         pch=16, col='aquamarine', cex=2)
}

legend('topright',pch=c(16,16,3,16), col=c(rgb(red=0, green=0, blue=1, alpha=0.5),rgb(red=1, green=0, blue=0, alpha=0.5),'maroon','aquamarine'),
       legend=c('extant','ancestral','root', living_fossil), cex=1.75)
title(paste('Reconstruction of first ', paste(max_PC, 'PCs'), '; motif ', motif_ind, sep=''), cex.main=2.25)










