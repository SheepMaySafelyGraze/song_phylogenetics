### Motif analysis
library(ape)
library(dplyr)
library(taxize)
library(phytools)
library(phangorn)
### load data and tree
trait_data <- read.csv('phylogenetic_analysis/namefix_trait_data_91.csv')
tree <- read.tree('phylogenetic_analysis/namefix_tree_91.tre')


### visualising data

n_species <- length(unique(trait_data$species))  # 100
n_families <- length(unique(trait_data$family))  # 80
n_genus <- length(unique(trait_data$gen))        # 95

# todo: resolve final 9 conflicts
print(length(intersect(tree$tip.label, trait_data$species)))

# for now trimming to common species
n_species = 91
tree <- drop.tip(tree, setdiff(tree$tip.label, trait_data$species))

drop_flag <- trait_data$species %in% tree$tip.label
trait_data <- trait_data[drop_flag, , drop=TRUE]

### plotting PCA components

# average trait data
trait_means <- trait_data %>%
  group_by(species) %>%
  summarise(PC1_mean=mean(PC1), PC2_mean=mean(PC2))

perm <- match(tree$tip.label, trait_means$species)

trait_means <- data.frame(trait_means[perm, ])
rownames(trait_means) <- trait_means$species

# plot for 2 selected clades

# root node = 92
# 'upper' clade node = 168
# chosen 'lower' clade node = 101

clade1 <- extract.clade(tree, 168)
clade2 <- extract.clade(tree, 101)
subsets <- c(clade1, clade2)
cols=c('red','blue')

par(mfrow = c(1, 2), oma = c(0, 0, 2, 0))
for (i in 1:2){
species_subset <- subsets[[i]]$tip.label

subset_tree <- keep.tip(tree, species_subset)
data_subset <- trait_means[species_subset, , drop=FALSE]
phylomorphospace(subset_tree,
                as.matrix(data_subset[,c('PC1_mean', 'PC2_mean')]),
                bty="l",ftype="off",node.by.map=FALSE,
                node.size=c(0,1.2),xlab="PC1",
                ylab="PC2", xlim=c(170, 350), ylim=c(-60, 50))
points(data_subset[,c('PC1_mean', 'PC2_mean')], col=cols[[i]], pch=20, cex=2)
  }
mtext("Phylomorphospace plots of 2 clades", outer = TRUE, cex = 1.5)
# TODO: plot subsets of data with phylomorphospace, to see if
# apparent PCA clusters correspond to species
## -> they do! with an interesting divergence of many species down to the 'bottom right'

### compute pagel's lambda and estimated covariance matrix

library(mvMORPH)

gmm_prob_dbl <- trait_data %>%
  group_by(species) %>%
  summarise(gmm_prob_0=mean(gmm_prob_0),
            gmm_prob_1=mean(gmm_prob_1),
            gmm_prob_2=mean(gmm_prob_2),
            gmm_prob_3=mean(gmm_prob_3),
            gmm_prob_4=mean(gmm_prob_4),
            gmm_prob_5=mean(gmm_prob_5),
            gmm_prob_6=mean(gmm_prob_6),
            gmm_prob_7=mean(gmm_prob_7),
            )

gmm_prob_dbl <- data.frame(gmm_prob_dbl)
rownames(gmm_prob_dbl) <- gmm_prob_dbl$species
gmm_prob_dbl <- gmm_prob_dbl[,2:8]
Y <- as.matrix(gmm_prob_dbl)
dat = list(Y=Y)

fit <- mvgls(Y~1, data=dat, tree, model='lambda', method='LL')









