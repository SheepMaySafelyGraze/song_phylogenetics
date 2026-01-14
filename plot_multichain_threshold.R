# making 8-dimensional motif-presence plot

library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)
library(ggtree)
library(ggnewscale)
library(viridis)


trait_data <- read.csv('phylogenetic_analysis/namefix_trait_data_100_1912.csv')
tree <- read.tree(file='phylogenetic_analysis/namefix_tree_100_1912.tre')

# names have been matched
geiger::name.check(tree, data.names=trait_data$species)

# create binary motif data
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


# sort motif_data so that entries are in the same order as tree labels
motif_data <- data.frame(motif_data)
species <- tree$tip.label

perm <- match(species, motif_data$species)

motif_data <- motif_data[perm,]
motif_data <- motif_data[,2:9, drop=TRUE]

# convert to 1/0
motif_data <- data.frame(lapply(motif_data, as.integer))
rownames(motif_data) <- species

# reading ancestral reconstruction results
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
colnames(motif_prob_plot_data)[3] <- 'Probability'

# numerically problematic nodes
problem_nodes <- c(146, 147, 148, 157, 158, 159, 167, 168)
grey_nodes <- data.frame(node=problem_nodes)
grey_plot_data <- left_join(grey_nodes, tree_df[tree_df$node %in% problem_nodes,
                                                c("node", "x", "y")], by="node")
grey_plot_data['fill'] = 'grey'



# plotting 8 together as a length-8 row
p +
  geom_tile(
    data = motif_prob_plot_data,
    aes(
      x = x + as.numeric(factor(motif)) * 0.3,
      y = y,
      fill = Probability
    ),
    height = 1,
    width  = 0.3
  ) +
  scale_fill_viridis_c() +
  theme_tree2() +
  scale_y_reverse() +
  ggtitle('Motif presence reconstruction') +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 21)) +
  new_scale_fill() +
  geom_tile(data=grey_plot_data,
            aes(x = x+1.35,
                y = y),
            fill='azure4',
            height = 1,
            width = 2.4,
            inherit.aes = FALSE)


### making individual motif plots

names = c('flat whistles',
          'slow trills',
          'fast trills',
          'chaotic songs',
          'ultrafast trills',
          'slow modulated whistles',
          'fast modulated whistles',
          'harmonic stacks')

for (i in 0:7) {
  par(mar=margin(1.2,1.2,1.2,1.2))
  plot(
    tree,
    show.tip.label = FALSE,
    main=paste('Motif ', i, ' (', names[i+1], ')', sep=''),
    cex.main=2,
    direction='rightwards'
  )
  
  present = 'firebrick'
  absent = 'steelblue'
  
  tip.cols <- ifelse(motif_dat_tips[,1+i] == 1, present, absent)
  tiplabels(
    pch = 21,
    bg  = tip.cols,
    cex = 2
  )
  
  node.ids <- (Ntip(tree) + 1):(Ntip(tree) + Nnode(tree))
  th_ancr <- motif_probs_thresh[101:199,]
  cex.scale <- 1 + 1 * th_ancr[,i+1]
  
  light_red = alpha(present, alpha=0.8)
  light_blue = alpha(absent, alpha=0.8)
  
  nodelabels(pie=th_ancr[,i+1], piecol=c(light_red, light_blue), cex=0.45)
  nodelabels(node=node.ids,pch=21, bg = ifelse(th_ancr[,i+1]>=0.5, light_red, light_blue), cex=1.8)
  
  nodelabels(
    pch = 21,
    bg  = "grey70",
    col = NA,
    cex = 4,
    node = problem_nodes
  )
  
  legend('bottomleft', legend=c('present', 'absent'), col=c(present, absent), cex=2, pt.cex=3, pch=20)
  
}













