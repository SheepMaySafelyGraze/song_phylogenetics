# D statistic test for phylogenetic signal
# does it vary significantly across clades?

library(ape)
library(phytools)
library(tidyverse)
library(taxize)
library(dplyr)
library(reshape2)
library(stringr)
library(adephylo)
library(caper)


# computing phylo d statistic for each motif and plotting

results <- list()

# phylo d computation
for (i in 0:7){
  results[[i+1]] <- phylo.d(data.frame(dat=motif_data[,i+1], species=species),
                            tree_node_removed, names.col=species, binvar=dat, permut=2e5)
}


# plotting in grid
plot.new()
par(mfrow = c(2, 4),
    oma = c(0, 0, 3, 0),  # outer margins for super title
    mar = c(3, 3, 2, 1))  # inner margins for per-plot titles

titles <- list()
for (i in 0:7){
  titles[[i+1]] <- paste('Motif', i)
}


for (i in 1:8) {
  plot(results[[i]])
  title(main = titles[[i]])
}

mtext("Distribution of D Statistics", outer = TRUE, cex = 1.5, font = 2)

par(xpd = NA)

legend(
  x = grconvertX(0.125, "npc", "user"),
  y = grconvertY(1, "npc", "user") + 2.7,
  legend = c('Random', 'Threshold'),
  col = c("red", "blue"),
  lty = 1,
  horiz = TRUE
)
par(mfrow=c(1,1), xpd=FALSE)



# comparison of D statistics for clades

# threshold model reconstruction results seem to show
# greater uncertainty in clade beginning 143, and less in clade beginning 109
# we can test this

clade_109 <- extract.clade(tree, 109)
clade_142 <- extract.clade(tree, 142)

data_109 <- motif_data[clade_109$tip.label, ]
data_142 <- motif_data[clade_142$tip.label, ]

clade_109$node.label <- NULL
clade_142$node.label <- NULL

phy_d_109 <- phylo.d(data.frame(dat=data_109[,8], species=clade_109$tip.label),
                     clade_109, binvar=dat, names.col = species, permut=2e5)
phy_d_142 <- phylo.d(data.frame(dat=data_142[,8], species=clade_142$tip.label),
                     clade_142, binvar=dat, names.col = species, permut=2e5)

