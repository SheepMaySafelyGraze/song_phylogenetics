##########################
#
#
# plotting functions
#
#
##########################

library(scales)
source('utils/data_handler.R')

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


plot_motif_at_tips <- function(tree, dat, motif_ind){
  if (!('motif_0' %in% colnames(dat))){
    dat <- motif_data_to_binary(tree, dat)[[2]]
  }
  
  plot(tree, show.tip.label=FALSE)
  
  colours <- ifelse(dat[,motif_ind+1]==1, 'red', 'blue')
  
  tiplabels(pch=16, col=colours, cex=1.25)
  legend('bottomleft', legend=c('present', 'absent'), pch=16, col=c('red','blue'),
         cex=1.5)
}



















