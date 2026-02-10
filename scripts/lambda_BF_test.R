##########################
#
#
# Bayes factor computation on lambda in univariate threshold model
#
#
##########################

library(ape)
library(stringr)

source('utils/data_handler.R')
data_path <- 'data/thinned_suboscine_data.csv'
tree_path <- 'data/thinned_suboscine_tree.tree'

dat <- read.csv(data_path)
tree <- read.tree(tree_path)

get_null_likelihood <- function(tree, X, motif_ind, thresh_prior_mean=0.0, thresh_prior_sd=10.0){
  # computes likelihood under null model lambda = 0
  # assumes a normal prior on threshold
  Sigma <- vcv.phylo(tree)
  dSigma <- diag(Sigma)
  
  
  res <- motif_data_to_binary(tree, X)
  
  tree <- res[[1]]
  motif_data <- res[[2]]
  
  X <- motif_data[,motif_ind+1]

  # define likelihood
  lhood_thresh <- function(dSig, thresh){
    out = 1.0
    for (i in 1:length(X)) {
      out <- out * pnorm(thresh, lower.tail = as.logical(X[i]==0),
                         mean=0, sd=sqrt(dSig[i]), log.p=FALSE)
    }
  out
  }
  
  cond_lhood <- \(z) dnorm(z,mean=thresh_prior_mean, sd=thresh_prior_sd, log=FALSE) *
                      lhood_thresh(dSigma, z)
  # integrating over prior on threshold
  lower = thresh_prior_mean - 5*thresh_prior_sd
  upper = thresh_prior_mean + 5*thresh_prior_sd
  ML <- log(integrate(cond_lhood, lower, upper, subdivisions=1e6, abs.tol=1e-4)$value)
  
  ML
}


get_ml_full_model <- function(tree, X, motif_ind, n_generations=3e6, burnin = 6e5, n_stones=50) {

  res <- motif_data_to_binary(tree, X)
  
  tree <- res[[1]]
  motif_data <- res[[2]]

  write.tree(tree, 'data/tree_tmp.tree')
  write_to_nexus(motif_data, motif_ind)
  
  # parameters
  model_type = 'threshold_alt' # available: Mk, threshold, threshold_alt
  n_chains=1
  analysis='ML'
  burnin_generations = burnin # for ML only
  tuning_interval = 1000 # for ML only
  printgen=500
  run_id = 'ML'

  ### populate with a run id, else date and time used
  
  if (is.null(run_id)) {
    run_id = date()
  }
  
  # folder of run results
  save_loc = paste('output/',
                   model_type,
                   '/',
                   run_id,
                   '/',
                   as.character(motif_ind),
                   sep=''
  )
  save_loc <- str_replace_all(save_loc, ':', '')
  
  # temporary Rev script
  tmp <- tempfile(fileext = '.Rev')
  
  con <- file(tmp, open='a')
  
  # write parameter values
  lines_0 <- c(
    paste('save_loc ="', save_loc, '"', sep=''),
    paste('n_chains =', n_chains, sep=''),
    paste('burnin_generations =', burnin_generations, sep=''),
    paste('n_generations =', n_generations, sep=''),
    paste('tuning_interval =', tuning_interval, sep=''),
    paste('printgen =', printgen, sep='')
  )
  
  writeLines(lines_0, con)
  
  # add model code
  lines_model <- readLines('rb-scripts/lines/threshold_lambda_ML.Rev')
  writeLines(lines_model, con)
  
  close(con)
  
  # run script
  system(paste('rb', tmp))
  
  unlink(tmp)
  print(paste('results at:', save_loc))
  return(save_loc)
}



get_lambda_bayes_factor <- function(tree, X, motif_inds, thresh_prior_mean=0.0,
                        thresh_prior_sd=10.0, n_generations=3e6){
  res <- list()
  for (motif_ind in motif_inds){
  ML_null <- get_null_likelihood(tree,X,motif_ind, thresh_prior_mean, thresh_prior_sd)
  res_file <- get_ml_full_model(tree,X, motif_ind, n_generations=n_generations)
  
  ML_model <- as.numeric(read_file(paste(res_file, '/ml_res.file', sep='')))
  print(paste(motif_ind, ':', ML_model - ML_null))
  res[[motif_ind+1]] <- (ML_model - ML_null)
  }
  res
}



get_null_bernoulli_likelihood <- function(tree, X, motif_ind,
                                          thresh_prior_mean=0.0, thresh_prior_sd=10.0){
  # computes marginal likelihood of null Bernoulli model, where
  # expects empirical samples from the prior on proportion
  
  average_tip_age <- mean(sapply(1:length(tree$tip.label), FUN= \(i) nodeheight(tree, i)))
  
  res <- motif_data_to_binary(tree, X)
  
  tree <- res[[1]]
  motif_data <- res[[2]]
  
  X_ <- motif_data[, motif_ind+1]
  m <- length(X_)
  n <- sum(X_)

  lhood <- \(p) ((p^n)*((1-p)^(m-n)))
  
  lhood_th <- \(th) dnorm(th, mean=thresh_prior_mean, sd=thresh_prior_sd)*(lhood(pnorm(th, mean=0, sd=sqrt(average_tip_age), lower.tail=FALSE)))

  lower = thresh_prior_mean - 5*thresh_prior_sd
  upper = thresh_prior_mean + 5*thresh_prior_sd

  ML <- log(integrate(lhood_th, lower, upper, subdivisions=1e6, abs.tol=1e-4)$value)

  ML
}

plot_induced_prop_prior <- function(tree, thresh_prior_mean=0.0, thresh_prior_sd=10.0){
  
  average_tip_age <- mean(sapply(1:length(tree$tip.label), FUN= \(i) nodeheight(tree, i)))
  
  x <- seq(thresh_prior_mean - 4*thresh_prior_sd, thresh_prior_mean + 4*thresh_prior_sd,
           0.01)
  
  plot(pnorm(x, mean=0, sd=sqrt(average_tip_age), lower.tail = TRUE), pnorm(x, mean=thresh_prior_mean, sd=thresh_prior_sd), type='l')
  
}

read_ML_from_file <- function(motif_ind, class){
  ML_result <- read_file(paste('output/threshold_alt/ML/', class, '/',
                               as.character(motif_ind),
                               '/ml_res.file',
                               sep=''))
  return(as.numeric(ML_result))
}

