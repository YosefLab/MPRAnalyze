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
    theta <- pmax(pmin(theta, 23), -23)
    
    d.est <- theta[1] + ddesign.mat %*% theta[-1] + log.ddepth
    
    ll <- sum(dgamma(x = dcounts,
                     shape = exp(theta[1]), # alpha
                     rate = exp(d.est),
                     log = TRUE))
    
    return(-ll)
}

#' @rdname ll.dna
ll.dna.ln.nb <- function(theta,
                         dcounts, log.ddepth,
                         ddesign.mat) {
    ## limit theta to avoid model shrinkage\explosion
    # theta <- pmax(pmin(theta, 23), -23)
    
    d.est <- ddesign.mat %*% theta[-1] + log.ddepth
    
    # the sum of log-likelihoods
    ll <- sum(dlnorm(x = dcounts,
                     meanlog = exp(d.est),
                     sdlog = exp(theta[1]),
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
    theta <- pmax(pmin(theta, 23), -23)
    
    log.d.est <- theta.d[,1] + theta.d[,-1] %*% t(ddesign.mat)
    if(NROW(theta.d)>1) {
        print(dim(theta.d)) 
        print(dim(t(rdesign.mat)))
    }
    log.r.est <- theta %*% t(rdesign.mat) + log.rdepth + log.d.est
    
    ## compute likelihood
    ll <- sum(dnbinom(x = rcounts,
                      size = exp(theta[1]),
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
    
    r.est <- ddesign.mat %*% theta.d[,-1] + rdesign.mat %*% theta + log.rdepth
    
    ## compute likelihood
    ll <- sum(dnbinom(x = rcounts,
                      size = exp(theta[1]),
                      mu = exp(r.est),
                      log = TRUE))
    
    return(-ll)
}