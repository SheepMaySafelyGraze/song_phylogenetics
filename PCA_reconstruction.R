### reconstruction of PCA components

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
PCA_data <- trait_data[,c(4,54,17:(17+36))]  # PCA components and gmm_cluster (motif)

# helper function

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

# summary statitics

raw_counts <- trait_data %>%
  count(species, gmm_cluster) %>%
  complete(species, gmm_cluster, fill = list(n = 0))

ggplot(raw_counts, aes(x=gmm_cluster, group=gmm_cluster, y=n)) + geom_boxplot() +
  labs(title='Number of recordings per species', x='Motif') +
  theme_minimal() + scale_x_continuous(breaks = unique(trait_data$gmm_cluster))
  

motif_counts <- trait_data %>%
  count(gmm_cluster) %>%
  group_by(gmm_cluster)

barplot(height = motif_counts$n, names.arg=motif_counts$gmm_cluster,
        xlab='Motif', ylab='n', main='Total number of recordings')

# what do the correlations look like between PCA coefficients?
cormat <- cor(PCA_data[, 3:39])
cormat[lower.tri(cormat)] <- NA

melt_cormat <- melt(cormat, na.rm=TRUE)

ggplot(data = melt_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Correlation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed() + labs(title='Correlation between PC coefficients')


########
#
# plot of tree, colouring points which have a particular motif
#
#
########

for(motif_ind in 0:7){
  motif_dat <- data.frame(trait_data %>%
    group_by(species) %>%
    summarise(motif_present = any(gmm_cluster==motif_ind)))
  print(motif_dat)
  rownames(motif_dat) <- motif_dat$species
  
  # reorder
  motif_dat <- motif_dat[match(tree$tip.label,motif_dat$species),]
  print(identical(motif_dat$species, tree$tip.label))
  
  tip_col <- ifelse(motif_dat$motif_present==TRUE, 'red','black')
  plot(tree, show.tip.label = FALSE)
  tiplabels(
    pch = 16,
    col = tip_col,
    cex = 1
  )

  title(paste('Species having motif', motif_ind))
}

# hard to tell whether there is phylogenetic signal just by looking!


# using mvMORPH for phylogenetic reconstruction under a BM model

# prepare data (aggregate)
agg_PCA <- PCA_data %>%
  group_by(species, gmm_cluster) %>%
  summarise(across(everything(), mean), .groups = 'drop')

# for each motif, modelling evolution and reconstructing
analyse_motif <- function(motif_ind, max_PC=5, method='PL-LOOCV'){
  i_agg_PCA <- data.frame(agg_PCA[agg_PCA$gmm_cluster==motif_ind,])
  rownames(i_agg_PCA) <- i_agg_PCA$species

  i_agg_PCA <- i_agg_PCA[,3:(3+max_PC-1)]

  # pruning tree to those species with an instance of motif
  i_tree <- drop.tip(tree, setdiff(tree$tip.label,rownames(i_agg_PCA)))

  Y <- as.matrix(i_agg_PCA)
  dat = list(Y=Y)

  fit <- mvgls(Y ~ 1, data = i_agg_PCA, tree=i_tree, model='lambda', method=method)

  reconst_i <- ancestral(fit)
  return(list(i_agg_PCA, reconst_i, i_tree, fit))
}

########
#
# TODO: understand differences between U-map and t-sne
#        make t-sne/UMAP plot, showing current species and ancestral species, as in
#        phylomorphoplot (maybe see source code for how best to do)
#
#         make sure these models are legit, testing different fitting methods etc.
#         + do significance tests for phylogenetic signal (lambda) in each
########

# function to make t-SNE plot of reconstructed values
plot_tsne <- function(reconst, motif_ind, use_tsne=TRUE){

  pruned_tree <- reconst[[3]]
  
  n <- dim(reconst[[1]])[1]
  N <- n + dim(reconst[[2]])[1]
  
  dists_to_root <- root_distances(pruned_tree)
  dists_to_root <- dists_to_root/max(dists_to_root)
  anc_dists = dists_to_root[(n+1):N]
  
  if (use_tsne==TRUE){
    reconst_tsne <- tsne(rbind(reconst[[1]], reconst[[2]]))
  } else {
      reconst_tsne <- rbind(reconst[[1]][,1:2], reconst[[2]][,1:2])
      }

  tip_tsne <- reconst_tsne[1:n, ]
  anc_tsne <- reconst_tsne[(n+1):N,]
  tip_tsne <- tip_tsne[match(pruned_tree$tip.label, rownames(reconst[[1]])),]
  plot(tip_tsne,
       pch = 21,
       bg = rgb(red=0, green=0, blue=1, alpha=0.5),
       asp = 1,
       xlab='',
       ylab='',
       cex=1.5)
  
  points(anc_tsne,
       pch = 21,
       bg = rgb(red=1, green=0, blue=0, alpha=0.5),
       asp = 1,
       xlab='',
       ylab='',
       cex=(anc_dists) + 0.5)
  
  
  for(i in 1:nrow(pruned_tree$edge)) {
    parent <- pruned_tree$edge[i, 1]
    child  <- pruned_tree$edge[i, 2]
    
    segments(reconst_tsne[parent, 1], reconst_tsne[parent, 2],
             reconst_tsne[child, 1],  reconst_tsne[child, 2],
             col = rgb(red=0, green=0, blue=0, alpha=0.1))
  }
  }

# function to make UMAP plot
plot_umap <- function(trait_PC, anc_PC, pruned_tree){
  n = length(pruned_tree$tip.label) # number of tips
  M = length(pruned_tree$node.label) + n # total nodes
  
  dists_to_root <- root_distances(pruned_tree)
  dists_to_root <- dists_to_root/max(dists_to_root)
  anc_dists = dists_to_root[(n+1):M]
  
  trait.umap <- umap(rbind(trait_PC, anc_PC))
  coords = trait.umap$layout
  
  plot(coords[1:n,],
       pch = 21,
       bg = rgb(red=0, green=0, blue=1, alpha=0.5),
       asp = 1,
       xlab='',
       ylab='',
       cex=1.5,
       xlim=c(min(coords[,1])-0.3, max(coords[,1])+0.3),
       ylim=c(min(coords[,2])-0.3, max(coords[,2])+0.3))
  
  points(coords[(n+1):M,],
         pch = 21,
         bg = rgb(red=1, green=0, blue=0, alpha=0.5),
         asp = 1,
         xlab='',
         ylab='',
         cex=1*(anc_dists) + 0.5)
  
  
  for(i in 1:nrow(pruned_tree$edge)) {
    parent <- pruned_tree$edge[i, 1]
    child  <- pruned_tree$edge[i, 2]
    
    segments(coords[parent, 1], coords[parent, 2],
             coords[child, 1],  coords[child, 2],
             col = rgb(red=0, green=0, blue=0, alpha=0.1))
  }
  
  return(trait.umap)
}

for (motif_ind in 7:7){
max_PC=5
reconst <- analyse_motif(motif_ind, max_PC=max_PC, method='LL')
n = nrow(reconst[[1]])
lambda <- reconst[[4]]$param

# t-sne
# use_tsne=TRUE
# plot_tsne(reconst, motif_ind, use_tsne=use_tsne)
# title(paste('Plot of', ifelse(use_tsne, paste('tsne for', max_PC, 'PCs'), 'first 2 PCs'), '; motif', motif_ind))

# UMAP
set.seed(37)
plot_map <- plot_umap(reconst[[1]], reconst[[2]], reconst[[3]])
par(cex.main=1.3, cex.axis=1.3)
title(paste('UMAP of first', paste(max_PC, 'PCs'), '; motif', motif_ind, " \n Pagel's λ:", round(lambda,2)))

# identify and highlight current taxa most similar to root taxa
points(plot_map$layout[n+1,1], plot_map$layout[n+1,2],pch=3, col='maroon', cex=3)

knn_matrix <- plot_map$knn$indexes[n+1,-1]

root_neighbours <- knn_matrix[knn_matrix %in% 1:n]  # extant taxa which are neighbours to root
living_fossil_ind <- root_neighbours[1]
living_fossil <- str_replace(rownames(reconst[[1]])[living_fossil_ind], '_', ' ')

points(plot_map$layout[living_fossil_ind,1],plot_map$layout[living_fossil_ind,2],
       pch=16, col='aquamarine', cex=2)
print(plot_map)
legend('topright',pch=c(16,16,3,16), col=c(rgb(red=0, green=0, blue=1, alpha=0.5),rgb(red=1, green=0, blue=0, alpha=0.5),'maroon','aquamarine'),
       legend=c('extant','ancestral','root', living_fossil))
}



