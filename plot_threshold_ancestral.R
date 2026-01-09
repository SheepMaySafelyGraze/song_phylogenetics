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

motif_0_prob_anc <- numeric(99)
motif_1_prob_anc <- numeric(99)
motif_2_prob_anc <- numeric(99)
motif_3_prob_anc <- numeric(99)
motif_4_prob_anc <- numeric(99)
motif_5_prob_anc <- numeric(99)
motif_6_prob_anc <- numeric(99)
motif_7_prob_anc <- numeric(99)

motif_0_prob_anc <- ancThresh(tree, setNames(motif_data[,1], species),
          model='BM', control=list(print=FALSE))$ace[,2]
motif_1_prob_anc <- ancThresh(tree, setNames(motif_data[,2], species),
                              model='BM', control=list(print=FALSE))$ace[,2]
motif_2_prob_anc <- ancThresh(tree, setNames(motif_data[,3], species),
                              model='BM', control=list(print=FALSE))$ace[,2]
motif_3_prob_anc <- ancThresh(tree, setNames(motif_data[,4], species),
                              model='BM', control=list(print=FALSE))$ace[,2]
motif_4_prob_anc <- ancThresh(tree, setNames(motif_data[,5], species),
                              model='BM', control=list(print=FALSE))$ace[,2]
motif_5_prob_anc <- ancThresh(tree, setNames(motif_data[,6], species),
                              model='BM', control=list(print=FALSE))$ace[,2]
motif_6_prob_anc <- ancThresh(tree, setNames(motif_data[,7], species),
                              model='BM', control=list(print=FALSE))$ace[,2]
motif_7_prob_anc <- ancThresh(tree, setNames(motif_data[,8], species),
                              model='BM', control=list(print=FALSE))$ace[,2]


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
      x = x + as.numeric(factor(motif)) * 0.1,
      y = y,
      fill = prob
    ),
    height = 1,
    width  = 0.8
  ) +
  scale_fill_viridis_c() +
  theme_tree2()




tree_ult <- force.ultrametric(tree, method='extend')
p_ult <- ggtree(tree_ult)

tree_ult_df <- p_ult$data

motif_ult_prob_plot_data <- left_join(motif_probs_expanded,
                                  tree_ult_df[, c("node", "x", "y")], by = "node")

motif_ult_prob_plot_data <- motif_prob_plot_data %>%
                          mutate(cell = factor(interaction(row, col)))

# as a 2x4 matrix

cell_cols <- palette()[1:8]

p_ult +
  geom_tile(
    data = motif_ult_prob_plot_data,
    aes(
      x = x + (col-2.5)*0.4,
      y = y + (row-1.5)*0.5,
      fill = cell,
      alpha = prob
    ),
    height = 0.5,
    width  = 0.4
  ) +
  scale_fill_manual(values = cell_cols) +
  scale_alpha(range = c(0.2, 1)) +
  guides(fill = guide_legend(title = "Motif"),
         alpha = guide_legend(title = "Prob"))


# as a 2x4 matrix  -  ultrametric

cell_cols <- palette()[1:8]

p_ult +
  geom_tile(
    data = motif_ult_prob_plot_data,
    aes(
      x = x + (col-2.5)*0.4,
      y = y + (row-1.5)*0.5,
      fill = prob
    ),
    height = 0.5,
    width  = 0.4
  ) +
  scale_fill_viridis_c()




# as a 2x4 matrix  -  non-ultrametric

cell_cols <- palette()[1:8]

p +
  geom_tile(
    data = motif_prob_plot_data,
    aes(
      x = x + (col-2.5)*0.4,
      y = y + (row-1.5)*0.5,
      fill = prob
    ),
    height = 0.5,
    width  = 0.4
  ) +
  scale_fill_viridis_c()

