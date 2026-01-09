### reconstruction of motif presence

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

# names have been matched
geiger::name.check(tree, data.names=trait_data$species)

# adding categorical variables: motif i present in each species

for (i in 0:7){
  motif_i_present = (trait_data$gmm_cluster == i)
  trait_data[paste('motif_', i,'_present', sep='')] = motif_i_present
}

motif_df <- trait_data %>%
  group_by(species) %>%
  summarise(motif_0 = any(motif_0_present),
            motif_1 = any(motif_1_present),
            motif_2 = any(motif_2_present),
            motif_3 = any(motif_3_present),
            motif_4 = any(motif_4_present),
            motif_5 = any(motif_5_present),
            motif_6 = any(motif_6_present),
            motif_7 = any(motif_7_present))

# sort motif_df so that entries are in the same order as tree labels

motif_df <- data.frame(motif_df)
species <- tree$tip.label

rownames(motif_df) <- species

perm <- match(tree$tip.label, motif_df$species)

motif_df <- motif_df[perm,]

# for each motif, plot in which species this motif is observed
for (i in 2:9){
  plot(tree, show.tip.label = F)
  
  motif_presence = motif_df[,i]
  tiplabels(bg = ifelse(motif_presence==1, 'red', 'blue'), pch=21, col='yellow')
}



motif_df <- motif_df[, 2:9 , drop=FALSE]
motif_df <- lapply(motif_df, as.integer)
motif_df <- data.frame(motif_df)



# fitting Mk model and inferring ancestral distribution

ordered_trans_mat <- matrix(c(
  0,1,
  2,0), 2,2,
  byrow=TRUE,
  dimnames = list(0:1,0:1))

motif_7_dat <- setNames(motif_df$motif_7,
                        species)
fit <- fitMk(tree, motif_7_dat,
      model = ordered_trans_mat,
      pi='estimated',
      rand.start=TRUE)

ancr_motif_7 <- ancr(fit)


cols<- c('red','blue')
plot(ancr_motif_7,args.nodelabels=list(piecol=cols)
     ,direction="upwards",
     args.plotTree=list(type="arc",arc_height=0.5))


# trying with a hidden rates model
fit_hrm_equal <- fitHRM(tree,motif_7_dat,ncat=2,
              pi="equal",parallel=TRUE,niter=30)

motif_7_hrm_ancr <- ancr(fit_hrm)$ace
motif_7_hrm_ancr_grouped <- data.frame('0' = motif_7_hrm_ancr[,1]
                                       + motif_7_hrm_ancr[,2],
                                       '1' = motif_7_hrm_ancr[,3]
                                       + motif_7_hrm_ancr[,4])

cols<-c(hcl.colors(6)[c(1,5)],
        c("0","1"))

ancr <- ancr(fit)
ancr$ace <- motif_7_hrm_ancr_grouped

plot(ancr,
     args.plotTree=list(type="arc", arc_height=0.5),
     args.nodelabels=list(piecol=cols, cex=0.25),
     args.tiplabels=list(piecol=cols))

# plot, with size showing probability
par(mar=margin(1.2,1.2,1.2,1.2))
plot(
  tree,
  show.tip.label = FALSE,
  main='Motif 7',
  cex.main=2
)

tip.cols <- ifelse(motif_7_dat == 1, "firebrick", "steelblue")
tiplabels(
  pch = 21,
  bg  = tip.cols,
  cex = 2
)

node.ids <- (Ntip(tree) + 1):(Ntip(tree) + Nnode(tree))
cex.scale <- 1 + 1 * motif_7_hrm_ancr_grouped[,2]

light_red = alpha("firebrick", alpha=0.5)
light_blue = alpha('steelblue', alpha=0.5)

nodelabels(node = node.ids,
           pch = 21,
           bg  = ifelse(motif_7_hrm_ancr_grouped[,2] >= 0.5, light_red, light_blue),
           cex = cex.scale)

legend('bottomleft', legend=c('present', 'absent'), col=c('firebrick', 'steelblue'), cex=2, pt.cex=3, pch=20)


# from this 2-state hidden rate model

# which parts are in which regime?
# -> plotting again and showing this
node.regimes <- data.frame('-' = motif_7_hrm_ancr[,1]
                                  + motif_7_hrm_ancr[,3],
                          '*' = motif_7_hrm_ancr[,2]
                                  + motif_7_hrm_ancr[,4])

plot(tree, show.tip.label = F)
nodelabels(node = node.ids, pch=21,
           bg=ifelse(node.regimes[,2]>=0.5, 'red', 'blue'))

# does the regime seem sensible in terms of two modes of evolution of a call?
# maybe ask the biologists on this

# how to assess the fit of this model? randomisation tests?


plot(tree, show.tip.label = F)
nodelabels(node = node.ids, pch=21,
           bg=ifelse(motif_7_hrm_ancr[,4]>=0.5, 'red', 'blue'))

tiplabels(bg = ifelse(motif_7_dat==1, 'red', 'blue'), pch=21, col='yellow')






# repeating analysis with probabilities (encorporating uncertainty)

# for each species, the probability of it having the motif is the probability that
# at leas one recording is in that motif, i.e. 1 less the product of 1-p_i for
# probability that record i is in the particular motif

print(colnames(trait_data)[55])
prob_data <- trait_data[,55:(55+7)]

combine_probabilities <- function(prob_vec){
  out <- 1
  for (p_i in prob_vec){
    out <- out*(1-p_i)
  }
  return(1 - out)
}

motif_probabilities <- trait_data %>%
  group_by(species) %>%
  summarise(motif_0_prob = combine_probabilities(gmm_prob_0),
            motif_1_prob = combine_probabilities(gmm_prob_1),
            motif_2_prob = combine_probabilities(gmm_prob_2),
            motif_3_prob = combine_probabilities(gmm_prob_3),
            motif_4_prob = combine_probabilities(gmm_prob_4),
            motif_5_prob = combine_probabilities(gmm_prob_5),
            motif_6_prob = combine_probabilities(gmm_prob_6),
            motif_7_prob = combine_probabilities(gmm_prob_7))

max_motif_probabilities <- trait_data %>%
  group_by(species) %>%
  summarise(motif_0_prob = max(gmm_prob_0),
            motif_1_prob = max(gmm_prob_1),
            motif_2_prob = max(gmm_prob_2),
            motif_3_prob = max(gmm_prob_3),
            motif_4_prob = max(gmm_prob_4),
            motif_5_prob = max(gmm_prob_5),
            motif_6_prob = max(gmm_prob_6),
            motif_7_prob = max(gmm_prob_7))

motif_probabilities <- data.frame(motif_probabilities)
max_motif_probabilities <- data.frame(max_motif_probabilities)

hist(motif_probabilities[,1])
