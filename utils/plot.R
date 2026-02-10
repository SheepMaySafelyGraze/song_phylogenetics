##########################
#
#
# plotting functions for ancestral states
#
#
##########################

library(scales)
source('utils/data_handler.R')

thresh_ace_to_binary <- function(liab_df, thresh_df){
  # converts liabilities to binary traits
  
  out_df <- liab_df
  
  if (nrow(out_df) != nrow(thresh_df)){
    warning('differing number of samples of liabilities and thresholds')
    return(NULL)
  }
  
  for (i in 1:nrow(out_df)){
    th_i <- (thresh_df$threshold)[i]
    out_df[i, 4:ncol(out_df)] <- ifelse(out_df[i, 4:ncol(out_df)] > th_i, 1, 0)
  }

  # drop metadata columns
  out_df <- out_df[, c(5:ncol(out_df)), drop=TRUE]
  
  return(out_df)
}

plot_thresh_Mk <- function(df, motif_dat, tree){
  # assumes interior nodes are ordered in default RevBayes order
  # assumes df contains only sampled node states
  
  if ('Iteration' %in% colnames(df)){
    df <- df[,2:ncol(df)]
  }
  
  # reordering to match ape ordering
  perm <- match_ape_rb_nodes(tree)
  df <- df[,perm]

  ace <- colMeans(df)
  
  # colours
  present='red'
  absent='blue'
  
  plot(tree, show.tip.label=FALSE)
  
  tiplabels(pch=21, bg=ifelse(motif_dat == 1, present, absent), cex=2)
  
  node.ids <- (Ntip(tree) + 1):(Ntip(tree) + Nnode(tree))
  
  light_pres = alpha(present, alpha=0.8)
  light_abs = alpha(absent, alpha=0.8)
  
  nodelabels(pie=ace, piecol=c(light_pres, light_abs), cex=0.55)
  nodelabels(node=node.ids,pch=21, bg = ifelse(ace>=0.5, light_pres, light_abs), cex=1.2)
  
  return(ace)
}






















