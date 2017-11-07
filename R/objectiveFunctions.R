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
#' @param theta a vector of parameters to be optimized. First is the shape
#' parameter (alpha), the rest correspond to columns of the design matrix.
#' Note that the first column is the inverse of beta, since that makes computing
#' the estimate easier (numeric, K + 1)
#' @param dcounts the DNA counts (integer, N)
#' @param logdepth library size correction factors (numeric, N)
#' @param design.mat the design matrix (logical, N x K)
#'
#' @return the minus log likelihood of the specified model
ll.gammaDNA <- function(theta, dcounts, logdepth, design.mat) {
    # theta <- pmax(pmin(theta, 23), -23)

    alpha <- theta[1]

    beta.inv <- as.numeric((design.mat %*% theta[-1]) + logdepth)

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
    est <- as.numeric((design.mat %*% theta[-1]) + logdepth + log.dest)

    ## compute total likelihood
    ll <- sum(dnbinom(x = rcounts,
                      size = exp(theta[1]),
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
    theta.d <- c(alpha, theta[d.param.idx])
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