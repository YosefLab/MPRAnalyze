#' likelihood of dna observations
#' 
#' @name ll.dna
#' @rdname ll.dna
#' 
#' @aliases 
#' ll.dna.gamma.pois
#' ll.dna.ln.nb
#'
#' @param theta the vector of rna model parameters to evaluate likelihood for
#' (numeric, rna parameters)
#' @param dcounts the observed DNA counts (integer, enhancers x samples)
#' @param log.ddepth dna library size correction vector, log scale (numeric, samples)
#' @param ddesign.mat the dna model design matrix (logical, samples x dna parameters)
#'
#' @return negative log likelihood of dna  observations under the specified model
NULL

#' @rdname ll.dna
ll.dna.gamma.pois <- function(theta,
                              dcounts, log.ddepth,
                              ddesign.mat) {
    ## limit theta to avoid model shrinkage\explosion
    # theta <- pmax(pmin(theta, 23), -23)
    
    # (enhancers x parameters) x (parameters x  samples) = (enhancers x samples)
    beta.inv <- theta[,1] + theta[,-1] %*% t(ddesign.mat) + log.ddepth
    
    ll <- sum(dgamma(x = dcounts,
                     shape = exp(theta[,1]), # alpha
                     rate = exp(-beta.inv),
                     log = TRUE))
    
    return(-ll)
}

#' @rdname ll.dna
ll.dna.ln.nb <- function(theta,
                         dcounts, log.ddepth,
                         ddesign.mat) {
    ## limit theta to avoid model shrinkage\explosion
    # theta <- pmax(pmin(theta, 23), -23)
    
    # (enhancers x parameters) x (parameters x  samples) = (enhancers x samples)
    d.est <- theta[,-1] %*% t(ddesign.mat) + log.ddepth
    
    # the sum of log-likelihoods
    ll <- sum(dlnorm(x = dcounts,
                     meanlog = d.est,
                     sdlog = theta[,1],
                     log = TRUE))
    
    return(-ll)
}

#' likelihood of rna observations
#' 
#' @name ll.rna
#' @rdname ll.rna
#' 
#' @aliases 
#' ll.rna.gamma.pois
#' ll.rna.ln.nb
#'
#' @param theta the vector of rna model parameters to evaluate likelihood for
#' (numeric, rna parameters)
#' @param theta.d dna model parameters to condition likelihood on
#' @param rcounts the observed RNA counts (integer, enhancers x samples)
#' @param log.rdepth rna library size correction vector, log scale (numeric, samples)
#' @param ddesign.mat the dna model design matrix (logical, samples x dna parameters)
#' @param rdesign.mat the rna model design matrix (logical, samples x rna parameters)
#'
#' @return negative log likelihood of rna observations under the specified model
NULL

#' @rdname ll.rna
ll.rna.gamma.pois <- function(theta, theta.d, 
                              rcounts, log.rdepth,
                              ddesign.mat, rdesign.mat) {
    ## limit theta to avoid model shrinkage\explosion
    # theta <- pmax(pmin(theta, 23), -23)
    
    # (enhancers x parameters) x (parameters x  samples) = (enhancers x samples)
    log.d.est <- log(theta.d[,1] + theta.d[,-1] %*% t(ddesign.mat))
    log.r.est <- theta[,-1] %*% t(rdesign.mat) + log.rdepth + log.d.est
    
    ## compute likelihood
    ll <- sum(dnbinom(x = rcounts,
                      size = exp(theta[,1]),
                      mu = exp(log.r.est),
                      log = TRUE))
    
    return(-ll)
}

#' @rdname ll.rna
ll.rna.ln.nb <- function(theta, theta.d, 
                         rcounts, log.rdepth,
                         ddesign.mat, rdesign.mat) {
    ## limit theta to avoid model shrinkage\explosion
    # theta <- pmax(pmin(theta, 23), -23)
    
    # (enhancers x parameters) x (parameters x  samples) = (enhancers x samples)
    r.est <- theta.d[,-1] %*% t(ddesign.mat) + theta[,-1] %*% t(rdesign.mat) + log.rdepth
    
    ## compute likelihood
    ll <- sum(dnbinom(x = rcounts,
                      size = exp(theta[,1]),
                      mu = exp(r.est),
                      log = TRUE))
    
    return(-ll)
}