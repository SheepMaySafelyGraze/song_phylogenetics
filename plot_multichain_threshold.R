# making 8-dimensional motif-presence plot

# TODO:
# - run reconstruction for each motif
# - make tree ultrametric for ease of plotting
# 
# - combine results for each node
# - for each node, figure out how to plot
# 
#
#

library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)
library(ggtree)

motif_0_prob_anc <- numeric(99)
motif_1_prob_anc <- numeric(99)
motif_2_prob_anc <- numeric(99)
motif_3_prob_anc <- numeric(99)
motif_4_prob_anc <- numeric(99)
motif_5_prob_anc <- numeric(99)
motif_6_prob_anc <- numeric(99)
motif_7_prob_anc <- numeric(99)

motif_0_prob_anc <- read.csv(paste('results/motif_', 0,'_ace.cex', sep=''))[,3]
motif_1_prob_anc <- read.csv(paste('results/motif_', 1,'_ace.cex', sep=''))[,3]
motif_2_prob_anc <- read.csv(paste('results/motif_', 2,'_ace.cex', sep=''))[,3]
motif_3_prob_anc <- read.csv(paste('results/motif_', 3,'_ace.cex', sep=''))[,3]
motif_4_prob_anc <- read.csv(paste('results/motif_', 4,'_ace.cex', sep=''))[,3]
motif_5_prob_anc <- read.csv(paste('results/motif_', 5,'_ace.cex', sep=''))[,3]
motif_6_prob_anc <- read.csv(paste('results/motif_', 6,'_ace.cex', sep=''))[,3]
motif_7_prob_anc <- read.csv(paste('results/motif_', 7,'_ace.cex', sep=''))[,3]


motif_probs_thresh <- data.frame(motif_0 = motif_0_prob_anc,
                                 motif_1 = motif_1_prob_anc,
                                 motif_2 = motif_2_prob_anc,
                                 motif_3 = motif_3_prob_anc,
                                 motif_4 = motif_4_prob_anc,
                                 motif_5 = motif_5_prob_anc,
                                 motif_6 = motif_6_prob_anc,
                                 motif_7 = motif_7_prob_anc)



motif_dat_tips <- data.frame(lapply(motif_data, as.numeric))

motif_probs_thresh <- rbind(motif_dat_tips, motif_probs_thresh)

# as a list
list_motif_probs_thresh <- list(199)
# tip values
for (i in 1:100){
  list_motif_probs_thresh[[i]] = matrix(motif_dat_tips[i,],
                                      2, 4, byrow=TRUE)
}
# ancestral nodes
for(i in 101:199){
  list_motif_probs_thresh[[i]] = matrix(motif_probs_thresh[i-100,],
                                            2,4, byrow=TRUE)
}

# as an expanded grid

motif_probs_expanded <- expand.grid(
  node = 1:(Ntip(tree) + tree$Nnode),
  motif = 0:7
)
motif_probs_expanded$prob <- numeric(8*199)
motif_probs_expanded$row <- numeric(8*199)
motif_probs_expanded$col <- numeric(8*199)

# populating with probabilities
cols <- rep(1:4, 2)
rows <- c(rep(1,4), rep(2,4))
for (i in 0:7){
  
  for (j in 1:199){
    motif_probs_expanded[j+i*199, 'row'] <- rows[i+1]
    motif_probs_expanded[j+i*199, 'col'] <- cols[i+1]
    motif_probs_expanded[j+i*199,'prob'] <- motif_probs_thresh[j,i+1]
  }
}


p <- ggtree(tree)

tree_df <- p$data

motif_prob_plot_data <- left_join(motif_probs_expanded,
                                               tree_df[, c("node", "x", "y")], by = "node")
# as a length-8 row
p +
  geom_tile(
    data = motif_prob_plot_data,
    aes(
      x = x + as.numeric(factor(motif)) * 0.3,
      y = y,
      fill = prob
    ),
    height = 1,
    width  = 0.3
  ) +
  scale_fill_viridis_c() +
  theme_tree2() +
  scale_y_reverse() +
  ggtitle('Motif presence reconstruction') +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 21))




### making individual motif plots


# motif 2 (fast trills)

i=2
par(mar=margin(1.2,1.2,1.2,1.2))
plot(
  tree,
  show.tip.label = FALSE,
  main=paste('Motif', i),
  cex.main=2
)

tip.cols <- ifelse(motif_dat_tips[,1+i] == 1, "firebrick", "steelblue")
tiplabels(
  pch = 21,
  bg  = tip.cols,
  cex = 2
)

node.ids <- (Ntip(tree) + 1):(Ntip(tree) + Nnode(tree))
th_ancr <- motif_probs_thresh[101:199,]
cex.scale <- 1 + 1 * th_ancr[,i+1]

light_red = alpha("firebrick", alpha=0.8)
light_blue = alpha('steelblue', alpha=0.8)

nodelabels(pie=th_ancr[,i+1], piecol=c(light_red, light_blue), cex=0.5)
nodelabels(node=node.ids,pch=21, bg = ifelse(th_ancr[,i+1]>=0.5, light_red, light_blue), cex=0.8)

legend('bottomleft', legend=c('present', 'absent'), col=c('firebrick', 'steelblue'), cex=2, pt.cex=3, pch=20)


# motif 7

i=7
par(mar=margin(1.2,1.2,1.2,1.2))
plot(
  tree,
  show.tip.label = FALSE,
  main=paste('Motif', i),
  cex.main=2,
  direction='rightwards'
)
axisPhylo()

tip.cols <- ifelse(motif_dat_tips[,1+i] == 1, "firebrick", "steelblue")
tiplabels(
  pch = 21,
  bg  = tip.cols,
  cex = 2
)

node.ids <- (Ntip(tree) + 1):(Ntip(tree) + Nnode(tree))
th_ancr <- motif_probs_thresh[101:199,]
cex.scale <- 1 + 1 * th_ancr[,i+1]

light_red = alpha("firebrick", alpha=0.8)
light_blue = alpha('steelblue', alpha=0.8)

nodelabels(pie=th_ancr[,i+1], piecol=c(light_red, light_blue), cex=0.5)
nodelabels(node=node.ids,pch=21, bg = ifelse(th_ancr[,i+1]>=0.5, light_red, light_blue), cex=0.8)

legend('bottomleft', legend=c('present', 'absent'), col=c('firebrick', 'steelblue'), cex=2, pt.cex=3, pch=20)







