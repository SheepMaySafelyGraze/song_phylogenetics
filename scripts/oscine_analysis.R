##########################
#
#
# script for producing oscine results and figures
#
#
##########################

library(pheatmap)
library(grid)
library(gridExtra)


# load utility functions
utils <- list.files('utils')
for (file in utils) {
  if (!grepl('.Rev', file)){
    source(paste('utils/', file, sep=''))
  }
}

# load data and tree
data_path <- 'data/final/max_genera_100_osc_dat.csv'
tree_path <- 'data/final/max_genera_100_tree.tree'

data <- read.csv(data_path)
tree <- read.tree(tree_path)

# convert gmm clusters binary motif presence, and align rows with tree tip labels
motif_data <- motif_data_to_binary(tree, data)[[2]]


### 1 - identify significant motifs via a permutation test with Moran's I
set.seed(2026)

p_vals <- numeric(8)

for (i in 0:7) {
  moran_res_i <- moran_I_permutation_test(tree, motif_data[, i+1], n_samps=10000)
  print(paste('Motif ', i, '  I: ', round(moran_res_i[[1]],3),
              '   p: ', round(moran_res_i[[2]],3), sep=''))
  p_vals[i+1] <- moran_res_i[[2]]
}


significant_motifs <- which(p_vals<0.05) - 1


### plotting recordings of motifs in PCA space

pca_inds <- 17:(17+37)  # columns of PCA data and gmm cluster
plot_PCA_centroids(data[,pca_inds], significant_motifs, main='Oscine Motifs',
                   cols=c(rainbow(4)[c(1,2)], 'grey50', rainbow(4)[3]))


### plotting correlation estimates and PCA distance matrices

# reading log of runs
runs <- read.csv('output/runs.csv', colClasses = c(motif="character"))

cov_mat <- read_cor_matrix(list(model_type='MV_threshold', motif=paste(significant_motifs, sep='', collapse=''), data=data_path), runs, join=TRUE)
pca_mat <- get_PCA_distance_matrix(data[,pca_inds], significant_motifs)

rownames(cov_mat) <- significant_motifs
colnames(cov_mat) <- significant_motifs

rownames(pca_mat) <- significant_motifs
colnames(pca_mat) <- significant_motifs



# plot cov and pca matrices via corrplot

corrplot(cov_mat, method='circle', type='upper', diag=FALSE, addCoef.col = 'black',
         number.cex=1.75, cl.cex=1.2, tl.cex=1.75, tl.col='black', tl.srt=30,
         title='Correlation estimates', mar=c(1,1,3.5,1), cex.main=2, col=COL2("RdBu", 200))

corrplot(pca_mat, is.corr=FALSE, method='circle', type='upper', diag=FALSE, addCoef.col = 'black',
         number.cex=1.75, cl.cex=1.2, tl.cex=1.75, tl.col='black', tl.srt=30,
         col=colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu"))[c(2,3,5,9)])(200), col.lim=c(0,120),
         title='PCA distances', mar=c(1,1,3.5,1), cex.main=2)





pca_mat[lower.tri(pca_mat, diag=TRUE)] <- NA
cov_mat[lower.tri(pca_mat, diag=TRUE)] <- NA

# define breaks to align colours 
breaks_cov <- seq(min(cov_mat[upper.tri(cov_mat)]), max(cov_mat[upper.tri(cov_mat)]), length.out = 100)
breaks_pca <- seq(min(pca_mat[upper.tri(pca_mat)]), max(pca_mat[upper.tri(pca_mat)]), length.out = 100)


# define heatmap objects
cov_mat <- pheatmap(cov_mat,
                    color = colorRampPalette(c("lightblue", "red"))(99),
                    breaks = breaks_cov,
                    display_numbers = TRUE, cluster_rows=FALSE, cluster_cols=FALSE,
                    fontsize=16, fontsize_number = 18, silent=TRUE,
                    number_color = 'white', na_col='white', main='Correlation',
                    main.cex=5,
                    number_format='%.2f')

pca_mat <- pheatmap(pca_mat,
                    color = colorRampPalette(c("red", "lightblue"))(99),
                    breaks = breaks_pca,
                    display_numbers = TRUE, cluster_rows=FALSE, cluster_cols=FALSE,
                    fontsize=16, fontsize_number = 18, silent=TRUE,
                    number_color='white', na_col='white', main='PCA distance',
                    main.cex=5,
                    number_format='%.1f')

# plot in grid
grid.arrange(cov_mat$gtable, pca_mat$gtable, ncol = 2, top=textGrob('Oscine motifs', gp=gpar(fontsize=24)))


# confirm this with correlation from Mk model
