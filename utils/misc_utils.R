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


simulate_threshold_model <- function(tree, n_samps, lambdas) {
  ### todo
  warning('Not implemented')
  return(NULL)
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
  res <- optimize(f=lhood, interval=c(-0.2,max_lambda), maximum=TRUE)
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


