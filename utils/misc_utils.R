##########################
#
#
# misc. utilities
#
#
##########################

library(RcppAlgos)
library(phytools)
library(ape)
library(cauphy)
library(matrixcalc)
library(PNC)
library(expm)

source('utils/plot.R')

get_weighted_hamming <- function(tree, X, average=TRUE){
  # for a single binary trait, computed weighted hamming distance
  # this is simply the average (or total) distance 
  # between points with differing trait values
  
  tot <- 0
  n <- length(tree$tip.label)
  D <- cophenetic.phylo(tree)
  
  # all unique pairs of elements
  pairs <- RcppAlgos::comboGeneral(1:n, m=2, repetition=FALSE)
  
  for (i in 1:nrow(pairs)){
    
    j = pairs[i, 1]
    k = pairs[i, 2]

    if (X[j] != X[k]) {
      tot = tot + D[j,k]
    }
    
  }
  
  tot / ifelse(average, n*(n-1), 1)
}


simulate_threshold_model <- function(tree, thresholds, lambdas=NULL) {
  # simulates from (uncorrelated) threshold model with number of dimensions
  # specified by number of thresholds provided
  
  if (is.null(lambdas)) {
    lambdas=rep(1, length(thresholds))
  } else {
    if (length(lambdas) != length(thresholds)) {
      warning('incompatible parameters provided')
      return(NULL)
    }
  }
  
  out_mat <- matrix(NA, nrow=Ntip(tree), ncol=length(thresholds))

  for (i in 1:length(thresholds)){
    lam_i <- lambdas[i]
    th_i <- thresholds[i]
    
    sig <- vcv.phylo(tree)
    diag_sig <- diag(diag(sig))
    sig <- lam_i*(sig-diag_sig) + diag_sig
    eps <- rnorm(Ntip(tree))
    liab <- sqrtm(sig)%*%eps  # liability via transform of std normal
    
    out_mat[,i] <- as.integer(liab>th_i)
  }
  
  return(data.frame(out_mat))
}

compute_unbounded_lambda <- function(tree, x, test=TRUE, estimate=FALSE) {
  # computed unbounded estimate of pagel's lambda via maximum likelihood,
  # if estimate=FALSE, assumes mean of phylogenetic brownian motion is 0, and variance is 1, both known
  
  
  vcv_phylo <- vcv.phylo(tree, model='Brownian')
  vcv_diag <- diag(diag(vcv_phylo))

  
  # likelihood of data given lambda
  lhood <- function(lambda) {
    # define covariance matrix
    sig <- lambda*(vcv_phylo - vcv_diag)
    sig <- sig + vcv_diag
    invSig <- solve(sig)
    n <- ncol(sig)
    
    if(estimate){ # estimate mean and variance
      mu <- as.numeric(sum(invSig%*%x)/sum(invSig))
      sig2 <- as.numeric(t(x-mu)%*%invSig%*%(x-mu)/n)
    } else {
      mu <- 0.0
      sig2 <- 1.0
    }
    
    if (!is.positive.definite(sig)){
      return(-Inf)
    }
    
    # compute likelihood
    return(- t(x - mu)%*%((1/sig2)*invSig)%*%(x-mu)/2 - determinant(sig2 * sig, logarithm = TRUE)$modulus[1]/2 - n*log(2*pi)/2)
  }
  
  H <- nodeHeights(tree)
  max_lambda <- max(H[,2]/max(H[,1]))  # maximum possible scale factor
  
  
  # optimise likelihood over lambda
  # TODO: what is the theoretical minimum lambda?
  res <- optimize(f=lhood, interval=c(-0.05,max_lambda), maximum=TRUE)
  lam_ml <- res$maximum
  
  
  # test via chi-squared test
  if (test){
    stat <- 2*(res$objective - lhood(0))
    p_val <- pchisq(stat, df=1, lower.tail=FALSE, log.p=FALSE)
    
    return(list(lambda=lam_ml,
                logLik = res$objective,
                p=p_val))
  }
  
  return(list(lambda=lam_ml,
              logLik=res$objective))
}



moran_I_permutation_test <- function(tree, X, n_samps=2500) {
  # performs a one-sided randomisation test for signifcance of Moran I
  # on phylogenetic data X assumed binary
  
  w <- 1/cophenetic(tree)
  diag(w) <- 0
  
  obs <- Moran.I(X, w)$observed
  
  null_samps <- numeric(n_samps)
  
  for (n in 1:n_samps){
    X_n <- sample(X, length(X), replace=TRUE)

    if (var(X_n) > 0){
      null_samps[n] <- Moran.I(X_n, w)$observed
    }
    else {
      null_samps[n] <- NA
    }
  }
  
  return(list(observed=obs,
    p=mean(null_samps[!is.na(null_samps)] > obs),
    null_samples=null_samps
  ))
}



plot_lambda_against_moran_I <- function(tree, lambda_range=NULL,
                                           n_points = 400,
                                           n_samps_per_point=250,
                                        quiet=TRUE){
  # compares underlying lambda value in threshold model to induced Moran I value
  
  w <- 1/cophenetic(tree)
  diag(w) <- 0
  
  if (is.null(lambda_range)){
    H <- nodeHeights(tree)
    max_lambda <- max(H[,2]/max(H[,1])) 
    lambda_range <- c(0.0, 1.0)
  }
  lamb_grid <- seq(lambda_range[1],lambda_range[2],
                   (lambda_range[2]-lambda_range[1])/n_points)
  
  out <- numeric(length(lamb_grid))
  
  for (n in 1:length(out)){
    if (!quiet & n%%10==0){
      print(paste('progress:', 100*n/length(out), '/100'))
    }
    
    lamb_n <- lamb_grid[n]
    runs_n <- numeric(n_samps_per_point)
    for (m in 1:n_samps_per_point){
      sig <- vcv.phylo(tree)
      diag_sig <- diag(diag(sig))
      sig <- lamb_n*(sig-diag_sig) + diag_sig
      eps <- rnorm(Ntip(tree))
      liab <- sqrtm(sig)%*%eps  # liability via transform of std normal

      runs_n[m] <- Moran.I(as.integer(liab > 0), w)$observed
    }
    out[n] <- mean(runs_n)
  }
  
  plot(lamb_grid, out, type='l', xlab='Lambda', ylab='I', main="Moran's I, threshold model")
  
  return(list(MoranI=out,
              grid=lamb_grid))
  }



read_mv_ancestral_values <- function(tree, df, motif_inds, raw=FALSE, burnin=0.2,
                                     match_rb_order_to_ape = TRUE){
  # reads multivariate ancestral values
  # returns list of data frames for each motif, if raw: returns continuous values
  
  n <- length(motif_inds)
  m <- ncol(df)
  l <- nrow(df)
  n_taxa <- Ntip(tree)
  
  df <- df[floor(burnin*l):l,]
  
  thresholds <- df[, (m-n+1):m]
  
  df <- df[,5:m]
  
  perm <- match_ape_rb_nodes(tree)
  
  out <- vector('list', n)
  for (i in 1:n) {
    inds <- seq(n*n_taxa + i, n*(2*n_taxa - 3) + i, n)
    
    
    
    ace_i <- df[, inds]
    
    # add root node, which is not recorded and fixed at 0
    ace_i[paste('liability.', 2*n_taxa - 1, '..',i,'.', sep='')] <- rep(0.0, nrow(ace_i))
    
    if (!raw) {
      for (j in 1:nrow(ace_i)) {
        ace_i[j,] <- as.integer(ace_i[j,] > thresholds[j,i])
      
      }
    }
    
    if (match_rb_order_to_ape) {  # reorder to align 
      ace_i <- ace_i[,perm]
    }
    
    out[[i]] <- ace_i
  }
  return(out)
}












