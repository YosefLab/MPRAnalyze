# new likelihood functions
ll.dna.gamma.pois <- function(theta,
                              dcounts, logdepth,
                              ddesign.mat) {
    ## limit theta to avoid model shrinkage\explosion
    # theta <- pmax(pmin(theta, 23), -23)
    
    beta.inv <- as.numeric((ddesign.mat %*% theta[,-1]) + logdepth)
    
    ll <- sum(dgamma(x = dcounts,
                     shape = exp(theta[,1]),
                     rate = exp(-beta.inv),
                     log = TRUE))
    
    return(-ll)
}

ll.rna.gamma.pois <- function(theta, theta.d, 
                              rcounts, logdepth,
                              ddesign.mat, rdesign.mat) {
    ## limit theta to avoid model shrinkage\explosion
    # theta <- pmax(pmin(theta, 23), -23)
    
    ## get sample-specific estimators
    log.d.est <- log(theta.d[,1] + ddesign.mat %*% theta.d[,-1])
    r.est <- as.numeric((rdesign.mat %*% theta[,-1]) + logdepth + log.d.est)
    
    ## compute total likelihood
    ll <- sum(dnbinom(x = rcounts,
                      size = exp(theta[,1]),
                      mu = exp(r.est),
                      log = TRUE))
    
    return(-ll)
}

ll.dna.ln.nb <- function(theta,
                         dcounts, logdepth,
                         ddesign.mat) {
    ## limit theta to avoid model shrinkage\explosion
    # theta <- pmax(pmin(theta, 23), -23)
    
    d.est <- ddesign.mat %*% theta[,-1] + logdepth
    
    # the sum of log-likelihoods
    ll <- sum(dlnorm(x = dcounts,
                     meanlog = d.est,
                     sdlog = theta[,1],
                     log = TRUE))
    
    return(-ll)
}


ll.rna.ln.nb <- function(theta, theta.d, 
                         rcounts, logdepth,
                         ddesign.mat, rdesign.mat) {
    ## limit theta to avoid model shrinkage\explosion
    # theta <- pmax(pmin(theta, 23), -23)
    
    r.est <- ddesign.mat %*% theta.d[,-1] + rdesign.mat %*% theta[,-1] + logdepth
        
    ## compute total likelihood
    ll <- sum(dnbinom(x = rcounts,
                      size = exp(theta[,1]),
                      mu = exp(r.est),
                      log = TRUE))
    
    return(-ll)
}


#' likelihood objective function for lognormal point estimator for DNA counts
#'
#' @param theta the vector of parameters to fit (numeric, K + 1). First parameter
#' in the vector is the dispersion of the dstribution, the rest are ordered
#' in correspondence with the columns in the design matrix
#' @param dcounts the observed DNA counts (integer, N)
#' @param logdepth library size correction vector, log scale (numeric, N)
#' @param design.mat the design matrix (logical, N x K)
#'
#' @return the minus log likeihood under the specified model
ll.lnDNA <- function(theta, dcounts, logdepth, design.mat) {
    ## this increases speed and prevent warnings, but crashes Hessian computation
    # theta <- pmax(pmin(theta, 23), -23)
    
    # the estimated mean of the distribution
    est <- as.numeric((design.mat %*% theta[-1]) + logdepth)
    
    # the sum of log-likelihoods
    ll <- sum(dlnorm(x = dcounts,
                     meanlog = est,
                     sdlog = theta[1],
                     log = TRUE))
    return(-ll)
}

#' Likelihood objective function for gamma distribution for the DNA counts
#' 
#' K: number of parameter of enhancer-wise DNA model
#' N: number of sampkes
#' C: number of control enhancers fit (C+1: number of all enhancers fit)
#'
#' @param theta a vector of parameters to be optimized. First is the shape
#' parameter (alpha), the rest correspond to columns of the design matrix.
#' Note that the first column is the inverse of beta, since that makes computing
#' the estimate easier (numeric, K + 1). If control enhancers are co-fit,
#' theta, is a matrix of (C+1) x K
#' @param dcounts the DNA counts (integer, (C+1) x N)
#' @param logdepth library size correction factors (numeric, (C+1) x N)
#' @param design.mat the design matrix (logical, N x K)
#' Note that the design matrix is shared among all enhancers
#'
#' @return the minus log likelihood of the specified model
ll.gammaDNA <- function(theta, dcounts, logdepth, design.mat) {
    # theta <- pmax(pmin(theta, 23), -23)
    
    alpha <- theta[,1]
    
    beta.inv <- as.numeric((design.mat %*% theta[,-1]) + logdepth)
    
    ll <- sum(dgamma(x = dcounts,
                     shape = exp(alpha),
                     rate = exp(-beta.inv),
                     log = TRUE))
    
    return(-ll)
}

#' likelihood objective function for negative binomial RNA counts given a DNA
#' count estimator
#'
#' @param theta the parameter vector being optimized. First value is the
#' dispersion parameter, afterhich come factors ordered in correspondence
#' with the columns in the design.mat matrix (numeric, K + 1)
#' @param rcounts the RNA count vector (integer, N)
#' @param depth library size correction factors (numeric, N)
#' @param d.est the DNA point estimator to use, in log space (numeric, N)
#' @param design.mat the design matrix (logical, N x K)
#'
#' @return the minus log likelihood under the specified model
ll.nbRNA <- function(theta, rcounts, logdepth, log.dest, design.mat) {
    ## limit theta to avoid model shrinkage\explosion
    # theta <- pmax(pmin(theta, 23), -23)
    
    ## get sample-specific estimators
    est <- as.numeric((design.mat %*% theta[,-1]) + logdepth + log.dest)
    
    ## compute total likelihood
    ll <- sum(dnbinom(x = rcounts,
                      size = exp(theta[,1]),
                      mu = exp(est),
                      log = TRUE))
    
    return(-ll)
}

#' likelihood objective function for joint gamma-poisson mixture model
#'
#' @param theta the parameter vector being optimized (numeric, K + M + 1)
#' @param dcounts the DNA count vector (integer, N)
#' @param rcounts the RNA count vector (integer, N)
#' @param ddepth the DNA library size correction factors (numeric, N)
#' @param rdepth the RNA library size correction factors (numeric, N)
#' @param ddesign.mat the DNA design matrix (logical, N x K)
#' @param rdesign.mat the RNA design matrix (logical, N x M)
ll.mixture.gammaPoisson <- function(theta, dcounts, rcounts,
                                    log.ddepth, log.rdepth,
                                    ddesign.mat, rdesign.mat) {
    alpha <- theta[1]
    d.param.idx <- 2:(NCOL(ddesign.mat) + 1)
    theta.d <- t(matrix(alpha, theta[d.param.idx]))
    theta.r <- theta[-d.param.idx]
    
    
    d.ll <- ll.gammaDNA(theta = theta.d,
                        dcounts = dcounts,
                        logdepth = log.ddepth,
                        design.mat = ddesign.mat)
    
    ## ddepth is not accounted for in this estimate since it is a technical
    ## artifact that does not affect RNA counts
    log.dest <- alpha + (ddesign.mat %*% theta.d[-1])
    
    r.ll <- ll.nbRNA(theta = theta.r,
                     rcounts = rcounts,
                     logdepth = log.rdepth,
                     log.dest = log.dest,
                     design.mat = rdesign.mat)
    
    return(d.ll + r.ll)
}

#' likelihood objective function for joint gamma-poisson mixture model
#' on case and control enhancers
#' 
#' TODO: David - I wrote this so that this function can evalute likelihood
#' for one enhancer or one case together with controls with and without
#' DNA ctrl models prefit. Check whether this is ok, would be good to just
#' have one function in terms of code complexity if readability is not too bad?
#' Everything framed as matrix operations so I think it s still readable ish.
#' 
#' K: number of parameter of enhancer-wise DNA model
#' M: number of parameters of RNA model
#' Mctrl: number of extra parameters of RNA model for control enhancers
#' N: number of sampkes
#' C: number of control enhancers fit (C+1: number of all enhancers fit)
#'
#' @param theta the parameter vector that is optimized 
#' (numeric, (C+1) x K + M + Mctrl + 1)
#' @param dcounts the DNA count matrix (integer, (C+1) x N)
#' If control enhancer DNA models were prefit, this matrix only contains
#' the case enhancer (1 x N) as the LL terms with contribution of 
#' the control enhancer DNA observations are constant with respect to 
#' the parameters optimised in this scenario.
#' @param rcounts the RNA count matrix (integer, (C+1) x N)
#' @param ddepth the DNA library size correction factors (numeric, N)
#' @param rdepth the RNA library size correction factors (numeric, N)
#' @param ddesign.mat the DNA design matrix (logical, N x K)
#' @param rdesign.mat the RNA design matrix (logical, N x M)
#' @param rdesign.mat.correction design matrix for prefitconfounding effects
#' @param theta.r.correction prefit parameters to correct for confounding effects
ll.dnarna.gammaPoisson.one <- function(theta, dcounts, rcounts,
                                       log.ddepth, log.rdepth, rctrlscale,
                                       ddesign.mat, rdesign.mat, rctrldesign.mat) {
    nsamples <- length(log.ddepth)
    npar.d.beta <- NCOL(ddesign.mat)
    npar.r <- NCOL(rdesign.mat)
    
    # extract parameter vectors by model part
    alpha <- theta[npar.d.alpha]
    theta.d <- cbind(
        alpha, # alpha
        matrix(theta[seq(2, npar.d.alpha+npar.d.beta, by=1)],
               nrow=1, ncol=npar.d.beta, byrow=TRUE))
    theta.r <- theta[seq(npar.d.alpha+npar.d.beta+1, 
                         npar.d.alpha+npar.d.beta+npar.r), by=1)]
    
    d.ll <- ll.gammaDNA(theta = theta.d,
                        dcounts = dcounts,
                        logdepth = log.ddepth,
                        design.mat = ddesign.mat)
    
    ## ddepth is not accounted for in this estimate since it is a technical
    ## artifact that does not affect RNA counts
    log.dest <- alpha + ddesign.mat %*% theta.d[,-1]
    
    # estimate ll on case enhancer
    r.ll <- ll.nbRNA(theta = c(theta.r, rctrlscale),
                     rcounts = rcounts[1,,drop=FALSE],
                     logdepth = log.rdepth,
                     log.dest = log.dest[1,,drop=FALSE],
                     design.mat = rbind(rdesign.mat, rctrldesign.mat) )
    
    return(d.ll + r.ll)
}

ll.dna.gammaPoisson.one <- function(theta, dcounts, rcounts,
                                    log.ddepth, log.rdepth, rctrlscale,
                                    ddesign.mat, rdesign.mat, rctrldesign.mat) {
    nsamples <- length(log.ddepth)
    npar.d.beta <- NCOL(ddesign.mat)
    
    # extract parameter vectors by model part
    alpha <- theta[npar.d.alpha]
    theta.d <- cbind(
        alpha, # alpha
        matrix(theta[seq(2, npar.d.alpha+npar.d.beta, by=1)],
               nrow=1, ncol=npar.d.beta, byrow=TRUE))
    theta.r <- theta[seq(npar.d.alpha+npar.d.beta+1, 
                         npar.d.alpha+npar.d.beta+npar.r), by=1)]
    
    d.ll <- ll.gammaDNA(theta = theta.d,
                        dcounts = dcounts,
                        logdepth = log.ddepth,
                        design.mat = ddesign.mat)
    
    ## ddepth is not accounted for in this estimate since it is a technical
    ## artifact that does not affect RNA counts
    log.dest <- alpha + ddesign.mat %*% theta.d[,-1]
    
    # estimate ll on case enhancer
    r.ll <- ll.nbRNA(theta = c(theta.r, rctrlscale),
                     rcounts = rcounts[1,,drop=FALSE],
                     logdepth = log.rdepth,
                     log.dest = log.dest[1,,drop=FALSE],
                     design.mat = rbind(rdesign.mat, rctrldesign.mat) )
    
    return(d.ll + r.ll)
}

#' dcounts matrix 1 x sample: only case enhancer, assume controls are prefit
#' and given via theta.d.ctrl.prefit
ll.dnarna.gammaPoisson.wctrl <- function(theta, theta.d.ctrl.prefit,
                                         dcounts, rcounts,
                                         log.ddepth, log.rdepth, 
                                         ddesign.mat, rdesign.mat, rdesign.ctrl.mat) {
    
    nenhancers.obs <- NROW(rcounts) 
    nsamples <- length(log.rdepth)
    npar.d.beta <- NCOL(ddesign.mat)
    npar.r <- NCOL(rdesign.mat)
    npar.r.ctrl <- NCOL(rdesign.ctrl.mat)
    
    # extract parameter vectors by model part
    theta.d <- theta[seq(1, 1+npar.d.beta, by=1)] # alpha and beta
    theta.r <- theta[seq(1+npar.d.beta+1, 
                         1+npar.d.beta+npar.r), by=1)]
    theta.r.ctrl <- theta[
        seq(1+npar.d.beta+npar.r+1, 
            1+npar.d.beta+npar.r+npar.r.ctrl, by=1)]
    
    d.ll <- ll.gammaDNA(theta = theta.d,
                        dcounts = dcounts,
                        logdepth = log.ddepth, 
                        design.mat = ddesign.mat)
    
    ## ddepth is not accounted for in this estimate since it is a technical
    ## artifact that does not affect RNA counts
    # dimensions: (C+1) x N + ((C+1) x K) %*% t(N x K) = (C+1) x N
    # note that irrespective of whether or not control enhancer DNA models were
    # prefit, this matrix is (C+1) x N
    log.dest <- matrix(c(theta.d[1], theta.d.ctrl.prefit[,1,drop=FALSE]), 
                       nrow=nenhancers.obs, ncol=nsamples, byrow=FALSE) +
        (rbind(theta.d[,-1], theta.d.ctrl.prefit[,-1]) %*% t(ddesign.mat))
    
    # compute ll on case enhancer
    r.ll <- ll.nbRNA(
        theta = theta.r,
        rcounts = rcounts[1,,drop=FALSE],
        logdepth = log.rdepth,
        log.dest = log.dest[1,,drop=FALSE],
        design.mat = rdesign.mat )
    
    # compute ll on control enhancers
    r.ll <- r.ll + ll.nbRNA(
        theta = c(theta.r, theta.r.ctrl),
        rcounts = rcounts[-1,,drop=FALSE],
        logdepth = matrix(log.rdepth, nrow=nenhancers.obs-1, 
                          ncol=nsamples, byrow=TRUE),
        log.dest = log.dest[-1,,drop=FALSE],
        design.mat = cbind(rdesign.mat, rdesign.ctrl.mat) )
    
    return(d.ll + r.ll)
}

#' dcounts matrix 1 x sample: only case enhancer, assume controls are prefit
#' and given via theta.d.ctrl.prefit
ll.dna.wctrl <- function(theta, theta.r,
                         costfnDNA, costfnRNA,
                         dcounts, rcounts,
                         log.ddepth, log.rdepth,
                         ddesign.mat, rdesign.mat, rdesign.ctrl.mat) {
    
    # extract parameter vectors by model part
    theta.d <- theta[seq(1, 1+NCOL(ddesign.mat), by=1)] # alpha and beta
    
    d.ll <- costfnDNA(theta = theta.d,
                      dcounts = dcounts,
                      logdepth = log.ddepth, 
                      design.mat = ddesign.mat)
    
    # note that only the LL terms with observations of the case
    # enhancers have non-zero derivatives wrt to the case enhancer
    # dna model parameters
    log.dest <- theta.d[1] + ddesign.mat %*% theta.d[,-1]
    
    r.ll <- costfnRNA(
        theta = theta.r,
        rcounts = rcounts,
        logdepth = log.rdepth,
        log.dest = log.dest,
        design.mat = rdesign.mat )
    
    return(d.ll + r.ll)
}

#' dcounts matrix 1 x sample: only case enhancer, assume controls are prefit
#' and given via theta.d.ctrl.prefit
ll.rna.gammaPoisson.wctrl <- function(theta, theta.d, theta.d.ctrl.prefit,
                                      rcounts,
                                      log.rdepth, 
                                      ddesign.mat, rdesign.mat, rdesign.ctrl.mat) {
    
    nenhancers.obs <- NROW(rcounts) 
    nsamples <- length(log.rdepth)
    npar.r <- NCOL(rdesign.mat)
    npar.r.ctrl <- NCOL(rdesign.ctrl.mat)
    
    # extract parameter vectors by model part
    theta.r <- theta[seq(1, npar.r), by=1)]
    theta.r.ctrl <- theta[
        seq(npar.r+1, 
            npar.r+npar.r.ctrl, by=1)]
    
    # LL terms with DNA observation have zero derivative wrt to RNA model
    # and do not need to be computed here.
    
    ## ddepth is not accounted for in this estimate since it is a technical
    ## artifact that does not affect RNA counts
    # dimensions: (C+1) x N + ((C+1) x K) %*% t(N x K) = (C+1) x N
    # note that irrespective of whether or not control enhancer DNA models were
    # prefit, this matrix is (C+1) x N
    log.dest <- matrix(c(theta.d[1], theta.d.ctrl.prefit[,1,drop=FALSE]), 
                       nrow=nenhancers.obs, ncol=nsamples, byrow=FALSE) +
        (rbind(theta.d[,-1], theta.d.ctrl.prefit[,-1]) %*% t(ddesign.mat))
    
    # compute ll on case enhancer
    r.ll <- ll.nbRNA(
        theta = theta.r,
        rcounts = rcounts[1,,drop=FALSE],
        logdepth = log.rdepth,
        log.dest = log.dest[1,,drop=FALSE],
        design.mat = rdesign.mat )
    
    # compute ll on control enhancers
    r.ll <- r.ll + ll.nbRNA(
        theta = c(theta.r, theta.r.ctrl),
        rcounts = rcounts[-1,,drop=FALSE],
        logdepth = matrix(log.rdepth, nrow=nenhancers.obs-1, 
                          ncol=nsamples, byrow=TRUE),
        log.dest = log.dest[-1,,drop=FALSE],
        design.mat = cbind(rdesign.mat, rdesign.ctrl.mat) )
    
    return(d.ll + r.ll)
}

