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
#' 
#' @details matrix multiplications: object dimensions
#' theta: (dna param)
#' ddesign.mat: (samples x dna param-1)
#' log.ddepth, dcounts: (samples)
#' @details matrix multiplications: d.est
#' d.est <- theta[1] + ddesign.mat %*% theta[-1] + log.ddepth (gamma.pois)
#' d.est <- ddesign.mat %*% theta[-1] + log.ddepth (ln.nb)
#' dimensions( ddesign.mat %*% theta[-1] )
#' =(samples x dna param-1) x (dna param -1)
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
    theta <- pmax(pmin(theta, 23), -23)
    
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
#' 
#' @details matrix multiplications: object dimensions
#' theta: (rna param)
#' theta.d: (enhancers x dna param)
#' ddesign.mat: (samples x dna param-1)
#' rdesign.mat: (samples x rna param-1)
#' dcounts, rcounts, log.ddepth, log.rdepth: (samples)
#' @details matrix multiplications: log.d.est
#' log.d.est <- theta.d[,1] + theta.d[,-1] %*% t(ddesign.mat) (gamma.pois)
#' log.d.est <- theta.d[,-1] %*% t(ddesign.mat) (ln.nb)
#' dimensions( theta.d[,-1] %*% t(ddesign.mat) )
#' =(enhancers x dna param -1) x (samples x dna param-1)^T
#' =(enhancers x dna param -1) x (dna param-1 x samples)
#' =(enhancers x samples)
#' @details matrix multiplication: log.r.est
#' log.r.est <- log.d.est + (theta %*% t(rdesign.mat))[1,] + log.rdepth (gamma.pois)
#' log.r.est <- log.d.est + (theta %*% t(rdesign.mat))[1,] + log.rdepth (ln.nb)
#' dimensions( (theta %*% t(rdesign.mat))[1,] )
#' =((1 x rna param) x (samples x rna param)^T)[1,]
#' =(1 x samples)[1,] = (samples)
NULL

#' @rdname ll.rna
ll.rna.gamma.pois <- function(theta, theta.d, 
                              rcounts, log.rdepth,
                              ddesign.mat, rdesign.mat) {
    ## limit theta to avoid model shrinkage\explosion
    theta <- pmax(pmin(theta, 23), -23)
    
    log.d.est <- theta.d[,1] + theta.d[,-1] %*% t(ddesign.mat)
    #print(dim(rcounts))
    #print(length(log.rdepth))
    #print(dim(log.d.est))
    #print(dim(theta %*% t(rdesign.mat)))
    log.r.est <- log.d.est + 
        matrix((theta %*% t(rdesign.mat))[1,] + log.rdepth,
               nrow=NROW(rcounts), ncol=NCOL(rcounts), byrow=TRUE)
    
    ## compute likelihood
    ll <- sum(dnbinom(x = rcounts,
                      size = exp(theta.d[,1]),
                      mu = exp(log.r.est),
                      log = TRUE))
    
    return(-ll)
}

#' @rdname ll.rna
ll.rna.ln.nb <- function(theta, theta.d, 
                         rcounts, log.rdepth,
                         ddesign.mat, rdesign.mat) {
    ## limit theta to avoid model shrinkage\explosion
    theta <- pmax(pmin(theta, 23), -23)
    
    log.d.est <- theta.d[,-1] %*% t(ddesign.mat)
    log.r.est <- log.d.est + (theta[-1] %*% t(rdesign.mat))[1,] + log.rdepth
    
    ## compute likelihood
    ll <- sum(dnbinom(x = rcounts,
                      size = exp(theta[1]),
                      mu = matrix(exp(log.r.est) , nrow=NROW(rcounts), 
                                  ncol=length(log.r.est), byrow=TRUE),
                      log = TRUE))
    
    return(-ll)
}