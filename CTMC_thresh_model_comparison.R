# model comparison CTMC models vs Threshold model

library(ape)
library(phytools)
library(tidyverse)
library(taxize)
library(dplyr)
library(reshape2)
library(stringr)
library(adephylo)
library(mvtnorm)

trait_data <- read.csv('phylogenetic_analysis/namefix_trait_data_100_1912.csv')
tree <- read.tree(file='phylogenetic_analysis/namefix_tree_100_1912.tre')

# names have been matched
geiger::name.check(tree, data.names=trait_data$species)


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


# sort motif_df so that entries are in the same order as tree labels

motif_data <- data.frame(motif_data)
species <- tree$tip.label

perm <- match(species, motif_data$species)

motif_data <- motif_data[perm,]
rownames(motif_data) <- species
motif_data <- motif_data[,2:9, drop=TRUE]

# convert to 1/0
motif_data <- data.frame(lapply(motif_data, as.integer))
motif_data <- data.frame(lapply(motif_data, as.character))



######### 
# FOR EACH MOTIF
# FITTING CTMC AND THRESHOLD MODELS, AND COMPARING AIC
# MODELS ARE: MK, HRM, AND THRESHOLD  (fit via threshBayes)
#
#
#
#
#
#
#
#
#########


# there is an issue: ancThresh may be providing logLik for reconstructed internal
# nodes aswell, so recomuting for just tip values to check

get_tip_logLik <- function(tree, tip_liab, thresh_fit){
  dists <- distRoot(tree)
  n <- length(tree$tip.label)

  start <- floor(thresh_fit$burnin / thresh_fit$sample) + 1
  n_samples <- (thresh_fit$ngen - thresh_fit$burnin) / thresh_fit$sample
  
  l <- numeric(n)
  for (tip in 1:n) {
    # averaging over liabilities after burnin
    if ("ancThresh" %in% class(thresh_fit)){
      liabs <- tip_liab[start:(start+n_samples), tip]
      mu_liab <- mean(liabs)
    } else {
      mu_liab <- tip_liab[tip]
    }
    
    l[tip] <- mu_liab
  }
  sigma <- vcv.phylo(tree, model='Brownian')
  detSigma <- determinant(sigma, logarithm = TRUE)$modulus[1]
  invSigma <- solve(sigma)
  
  logLik_out <- -l%*%invSigma%*%l/2 -nrow(sigma)*log(2*pi) - detSigma/2
  
  return(logLik_out)
}




motif_inds = 0:7
thresh_aic <- numeric(8)
sym_aic <- numeric(8)
ard_aic <- numeric(8)
HRM_aic <- numeric(8)


for (i in motif_inds){
  print(paste('motif', i))
  aics_i = numeric(4)
  
  motif_i_data <- motif_data[, i+1]
  names(motif_i_data) <- species
  
  # CTMC models
  Mk_sym_fit_i <- fitMk(tree, motif_i_data, model='SYM', pi='fitzjohn')
  Mk_ard_fit_i <- fitMk(tree, motif_i_data, model='ARD', pi='fitzjohn')

  HRM_ard_fit_i <- fitHRM(tree, motif_i_data, ncat=2, model='ARD', pi='fitzjohn')
  
  # thresh model
  X <- one_hot_binary(motif_i_data)

  # thresh_fit_i <- ancThresh(tree, motif_i_data, ngens=1e+05)
  thresh_fit_alt <- fitThresh(tree, motif_i_data)
  
  sym_aic[i+1] <- AIC(Mk_sym_fit_i)
  ard_aic[i+1] <- AIC(Mk_ard_fit_i)
  HRM_aic[i+1] <- AIC(HRM_ard_fit_i)
  thresh_aic[i+1] <- AIC(thresh_fit_alt)
}


######### 
#
#
#
#
#comparing distributions of AICs under shuffles
#
# i.e. does either model benefit from the phylogenetic signal
# or are they both just going of the proportions in the population
#
#
#
######### 

get_null_AIC_dist <- function(i, tree, n_samps){
  
  ARD_samps <- numeric(n_samps)
  thresh_samples <- numeric(n_samps)
  
  AIC_out = data.frame(ARD=NA, thresh=NA)
  
  for (j in 1:n_samps){
    data_i <- sample(motif_data[,i+1], size=100, replace=TRUE)
    names(data_i) <- species
    
    fit_ARD <- fitMk(tree, data_i, model='ARD', pi='fitzjohn')
    fit_thresh <- fitThresh(tree, data_i)
    
    AIC_out[j,] <- c(AIC(fit_ARD), AIC(fit_thresh))
  }
  
  return(AIC_out)
}

for (i in 1:7){
  print(i)
  AIC_dist_i <- get_null_AIC_dist(i, tree, n_samps=100)
  
  data_i <- motif_data[,i+1]
  names(data_i) <- species
  
  Mk_ard_fit_i <- fitMk(tree, data_i, model='ARD', pi='fitzjohn')
  thresh_fit_alt <- fitThresh(tree, data_i)
  
  print(paste('motif',i))
  print(paste('Threshold model:', round(mean(AIC_dist_i$thresh < AIC(thresh_fit_alt)),4)))
  print(paste('ARD with estimated root:', round(mean(AIC_dist_i$ARD < AIC(Mk_ard_fit_i)),4)))
}


Mk_ard_fit_i <- fitMk(tree, motif_0_data, model='ARD')
thresh_fit_alt <- fitThresh(tree, motif_0_data)

mean(AIC_dist_0$thresh < AIC(thresh_fit_alt))
mean(AIC_dist_0$ARD < AIC(Mk_ard_fit_i))

motif_7_data <- motif_data[,8]
names(motif_7_data) <- species
AIC_dist_0 <- get_null_AIC_dist(7, tree, n_samps=75)

Mk_ard_fit_i <- fitMk(tree, motif_7_data, model='ARD')
thresh_fit_alt <- fitThresh(tree, motif_7_data)

mean(AIC_dist_0$thresh < AIC(thresh_fit_alt))
mean(AIC_dist_0$ARD < AIC(Mk_ard_fit_i))


######### 
#
#
#
#
# comparing leave-k-out prediction performance
#
#
#
#
#
#
######### 

get_number_of_pruned_descendants <- function(tree, tree_){
  kept_tips <- match(tree_$tip.label, tree$tip.label)
  if (length(kept_tips) != length(tree_$tip.label)){
    warning('pruned tree contains tips not in original tree')
  }
  
  descendants_full_tree <- Descendants(tree, node=1:(Ntip(tree) + tree$Nnode),
                                       type='tips')

  # list of number of descendants each original node has in pruned tree
  n_kept_tips <- lapply(descendants_full_tree, \(x) sum (x %in% kept_tips))
  
  return(list(descendants_full_tree, n_kept_tips))
}

get_corresponding_pruned_node <- function(node, tree, tree_, descendants_full){
  
  kept_tips <- tree_$tip.label

  
  descs <- descendants_full[[node]]
  descs <- tree$tip.label[descs]
                          
  kept_descendants <- intersect(descs, kept_tips)

  node_out <- getMRCA(tree_, kept_descendants)
  
  return(node_out)
}

get_pruned_ancestor_node <- function(node, tree, tree_pruned, n_descendants_pruned){
  # go up descendents of node until one appears in tree_pruned
  # returns the index in tree_pruned, and the distance to node
  
  # for each parent, check whether it has at least 2 children in pruned tree
  ancs <- Ancestors(tree, node, type='all')
  for (anc in ancs) {
    if (n_descendants_pruned[[anc]] >= 2) {
      return(anc)
    }
  }
  warning(paste('all ancestors of node', node, 'pruned, ancs:', ancs))
  return(NA)
}

get_ARD_model_leave_k_out_CV <- function(tree, data, k){
  # performs leave-k-out cross validation for the ARD model on motif i,
  # returns vector of actual values and their probabilities under the model
  
  # data must be a named vector of binary values
  
  obs_out <- numeric(k)
  model_prob_out <- numeric(k)
  
  data_ <- sample(data, size=length(data)-k, replace=FALSE)

  # dropped species removed from tree
  dropped_species <- setdiff(tree$tip.label, names(data_))
  tree_ <- drop.tip(tree, dropped_species)
  
  res <- get_number_of_pruned_descendants(tree, tree_)
  descendants_full <- res[[1]]
  n_descendants_pruned <- res[[2]]
  
  fit_pruned <- fitMk(tree_, data_, model='ARD', pi='fitzjohn')
  fit_ancr <- ancr(fit_pruned)$ace
  
  
  r0 <- fit_pruned$rates[1]
  r1 <- fit_pruned$rates[2]
  
  fit_rate_matrix <- matrix(c(-r1, r1, r0, -r0), byrow=TRUE, 2, 2)

  dropped_inds <- match(dropped_species, tree$tip.label)
  
  i=0
  for (s in dropped_inds) {
    i <- i+1
    # ancestor node in full tree which survided in pruned tree
    full_node <- get_pruned_ancestor_node(s, tree, tree_, n_descendants_pruned)
    pruned_node <- get_corresponding_pruned_node(full_node, tree,
                                                 tree_, descendants_full)

    anc_fit_s <- fit_ancr[pruned_node - length(data) + k,]

    dist_node_to_s <- dist.nodes(tree)[full_node, s]
    
    transition_matrix <- expm(fit_rate_matrix*dist_node_to_s)
    
    prob_s <- as.numeric(anc_fit_s %*% transition_matrix)

    dat_s <- as.integer(data[[s]])
    
    obs_out[i] <- dat_s
    model_prob_out[i] <- prob_s[dat_s+1]
  }
  
  return(data.frame(observed = obs_out, model_prob = model_prob_out))
}

get_thresh_probabilities_from_pruned_tree <- function(obs_values, node, tree, tree_, avg_liab,
                                       full_node, descendants_full, full_tree_distances){
  # we assume the liabilities to be fixed at their average MCMC values
  # TODO: extend to include variance at internal node and tip, see Grimmett and Stirzaker
  
  # finding closest relative of node in full tree which exists in pruned tree,
  pruned_node_descendants <- intersect(descendants_full[[full_node]],
                                       match(tree_$tip.label, tree$tip.label))
  
  dists <- full_tree_distances[node, pruned_node_descendants]
  closest_pruned <- pruned_node_descendants[which.min(dists)]
  
  # MCRA on full tree (point at which node departs from pruned tree)
  MRCA_ <- MRCA(tree, c(node, closest_pruned))
  print(full_node)
  print(MRCA_)
  
  # computing liability distribution for MRCA under Brownian bridge assumption
  L1 <- avg_liab[full_node]
  L2 <- avg_liab[closest_pruned]
  
  t <- full_tree_distances[full_node, closest_pruned]
  s <- full_tree_distances[full_node, MRCA_]
  mu <- (1 - s/t)*L1 + (s/t)*L2
  sigma2 <- (s/t)*(s-t)
  
  # getting probabilities for observed value, assuming 0 threshold value
  tip_mu <- mu
  tip_sigma2 <- sigma2
  
  p1 <- 1 - prnom(0, mu=tip_mu, sd=sqrt(tip_sigma2))
  
  return(c(1-p1, p1))
}

get_thresh_model_leave_k_out_CV <- function(tree, data, k){
  # performs leave-k-out cross validation for the thresh model on motif i,
  # returns vector of actual values and their probabilities under the model
  # data must be a named vector of binary values
  
  dists <- dist.nodes(tree)
  obs_out <- numeric(k)
  model_prob_out <- numeric(k)
  
  data_ <- sample(data, size=length(data)-k, replace=FALSE)
  
  # dropped species removed from tree
  dropped_species <- setdiff(tree$tip.label, names(data_))
  tree_ <- drop.tip(tree, dropped_species)
  print(paste('checking names:', geiger::name.check(tree_, data_)))
  res <- get_number_of_pruned_descendants(tree, tree_)
  descendants_full <- res[[1]]
  n_descendants_pruned <- res[[2]]
  
  fit_pruned <- ancThresh(tree_, data_, ngen=150000,model='BM', control=list(print=FALSE))
  
  # observed liabilities from MCMC procedure, accounting for burnin
  fit_ancr <- fit_pruned$liab[20:151,]
  avg_liabs <- colMeans(fit_ancr)
  
  dropped_inds <- match(dropped_species, tree$tip.label)
  
  i=0
  for (s in dropped_inds) {
    i <- i+1
    ancs <- Ancestors(tree, s)
    # ancestor node in full tree which survided in pruned tree
    full_node <- get_pruned_ancestor_node(s, tree, tree_, n_descendants_pruned)
    pruned_node <- get_corresponding_pruned_node(full_node, tree,
                                                 tree_, descendants_full)
    
    # following define the mean and variance of the model posterior distribution at s
    anc_liab <- mean(fit_ancr[, pruned_node - length(data) + k])
    dist_node_to_s <- dist.nodes(tree)[full_node, s]

    dat_s <- as.integer(data[[s]])
    # the mean value of the ancestral liability is assumed to be the same
    if(dat_s == 1){
      print(dat_s)
      prob_s <- 1 - pnorm(0, mean=anc_liab, sd=sqrt(dist_node_to_s))
    } else {
      prob_s <- pnorm(0, mean=anc_liab, sd=sqrt(dist_node_to_s))
    }
    
    obs_out[i] <- dat_s
    model_prob_out[i] <- prob_s
  }
  
  return(data.frame(observed = obs_out, model_prob = model_prob_out))
}


thresh_CV <- get_thresh_model_leave_k_out_CV(tree, motif_7_data, 10)
ARD_CV <- get_ARD_model_leave_k_out_CV(tree, motif_7_data, 10)


get_CV_score <- function(i, model='thresh', num_its, k){
  # model either 'thresh' or 'ARD'
  # computes average CV score
  
  data_i <- motif_data[,i+1]

  names(data_i) <- species
  
  if (model=='thresh') {
    CV_fun <- get_thresh_model_leave_k_out_CV
  } else {
    CV_fun <- get_ARD_model_leave_k_out_CV
  }
  
  out_kl <- numeric(num_its)
  for (i in 1:num_its){
    cv_res <- CV_fun(tree, data_i, k)
    out_kl[i] <- sum(-log(cv_res$model_prob))
  }
  return(out_kl)
}


ard_kl_means <- numeric(8)
thresh_kl_means <- numeric(8)

for (i in 0:7){
  ard_kl_mean <- mean(get_CV_score(i, model='ARD', num_its = 20, k=5))
  thresh_kl_mean <- mean(get_CV_score(i, model='thresh', num_its = 20, k=5)) 
  
  ard_kl_means[i+1] <- ard_kl_mean
  thresh_kl_means[i+1] <- thresh_kl_mean
  print(paste('motfi', i))
  print(paste('ARD:', round(ard_kl_mean, 2), 'Thresh', round(thresh_kl_mean, 2)))
}


ard_kl_means_2 <- numeric(8)
thresh_OU_kl_means <- numeric(8)

for (i in 0:7){
  ard_kl_mean <- mean(get_CV_score(i, model='ARD', num_its = 10, k=5))
  thresh_kl_mean <- mean(get_CV_score(i, model='thresh', num_its = 10, k=5)) 
  
  ard_kl_means_2[i+1] <- ard_kl_mean
  thresh_OU_kl_means[i+1] <- thresh_kl_mean
  print(paste('motfi', i))
  print(paste('ARD:', round(ard_kl_mean, 2), 'Thresh', round(thresh_kl_mean, 2)))
}


