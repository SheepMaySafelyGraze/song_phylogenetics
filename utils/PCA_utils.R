##########################
#
#
# PCA reconstruction utilities
#
#
##########################


get_PCA_distance_matrix <- function(PCA_df, motif_inds){
  # assumed PCA_df has final column for gmm cluster, and remaining are PCA components
  
  PCA_means <- PCA_df[,c(1:(ncol(PCA_df)))] %>%
    group_by(gmm_cluster) %>%
    summarise_all(mean)
  
  PCA_means <- PCA_means[PCA_means$gmm_cluster %in% motif_inds,2:ncol(PCA_means)]

  dist_PCA <- as.matrix(dist(PCA_means))
  
  print(PCA_means)
  
  return(dist_PCA)
}


plot_PCA_centroids <- function(PCA_df, motif_inds,
                               main='Distribution of Motifs', cols=NULL,
                               alpha=0.05){
  
  kdes = list()
  
  if (is.null(cols)){
    cols <- rainbow(length(motif_inds))
  }
  
  PCA_means <- PCA_df[,c(1:(ncol(PCA_df)))] %>%
    group_by(gmm_cluster) %>%
    summarise_all(mean)
  print(PCA_means)
  
  for (i in 1:length(motif_inds)){
    ind = motif_inds[i]
    
    PCA_ind <- PCA_df[PCA_df$gmm_cluster == ind, ]
    
    PCA_ind <- PCA_ind[, 1:2]

    kdes[[i]] <- MASS::kde2d(PCA_ind$PC1, PCA_ind$PC2)
  }
  
  xlim <- range(unlist(lapply(kdes, `[[`, "x")))
  ylim <- range(unlist(lapply(kdes, `[[`, "y")))
  
  contour(kdes[[1]]$x, kdes[[1]]$y, kdes[[1]]$z, col=cols[1],
          xlim = xlim, ylim = ylim, cex.axis=1.3, cex.lab = 1.3,
          xlab='PC1', ylab='PC2')
  PCA_ind <- PCA_df[PCA_df$gmm_cluster == motif_inds[1], ]
  points(PCA_ind$PC1, PCA_ind$PC2,bg=alpha(cols[[1]], alpha), col=alpha(cols[[1]], alpha), pch=21)
  
  for (i in 2:length(kdes)){
    
    contour(kdes[[i]]$x, kdes[[i]]$y, kdes[[i]]$z, col=cols[i],
            xlim = xlim, ylim = ylim, add=TRUE)
    
    PCA_ind <- PCA_df[PCA_df$gmm_cluster == motif_inds[i], ]
    points(PCA_ind$PC1, PCA_ind$PC2,bg=alpha(cols[[i]], alpha), pch=21, col=alpha(cols[[i]], alpha))
    
  }
  
  
  # add centroids
  for (i in 1:length(kdes)) {
    points(PCA_means[PCA_means$gmm_cluster == motif_inds[i],]$PC1,
           PCA_means[PCA_means$gmm_cluster == motif_inds[i],]$PC2, pch=23,
           bg=cols[[i]], col='black', cex=3)
  }
  
  legend("topright", legend = motif_inds, col = cols, lty = 1, title='Motif',
         cex=1.5)
  title(main, cex.main=1.85)
}


