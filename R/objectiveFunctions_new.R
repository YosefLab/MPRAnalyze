### Low level likelihoods for DNA/RNA observations

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

### High level likelihoods for optimisation

#' likelihood wrapper to compute terms of likelihood with non-zero derivative
#' with respect to case dna model and rna model 
#' if control enhancers are not used
#' 
#' @param theta the vector of rna model parameters to evaluate likelihood for
#' (numeric, rna parameters)
#' @param llfnDNA cost function to compute dna likelihood terms on (function)
#' @param llfnRNA cost function to compute rna likelihood terms on (function)
#' @param dcounts the observed case DNA counts (integer, enhancers x samples)
#' @param rcounts the observed RNA counts (integer, enhancers x samples)
#' @param log.ddepth dna library size correction vector, log scale (numeric, samples)
#' @param log.rdepth rna library size correction vector, log scale (numeric, samples)
#' @param rctrlscale vector of scaling parameters for rna model that correct for 
#' variation between conditions pre-fit on control enhancers 
#' (numeric, ctrl rna parameters)
#' @param ddesign.mat the dna model design matrix (logical, samples x dna parameters)
#' @param rdesign.mat the rna model design matrix (logical, samples x rna parameters)
#' @param rdesign.ctrl.mat the control rna model design matrix 
#' (logical, samples x rna parameters)
#' 
#' @return negative sum of log likelihood terms with non-zero derivative
#' with respect to case dna model and rna model 
#' if control enhancers are not used
cost.dnarna <- function(theta, dcounts, rcounts,
                        llfnDNA, llfnRNA,
                        log.ddepth, log.rdepth, rctrlscale,
                        ddesign.mat, rdesign.mat, rctrldesign.mat) {
    
    ## extract parameter vectors by model part
    # first parameter of DNA model is the variance(-link) parameter
    theta.d <- theta[seq(1, 1+NCOL(ddesign.mat), by=1)]
    theta.r <- theta[seq(1+NCOL(ddesign.mat)+1, 
                         1+NCOL(ddesign.mat)+NCOL(rdesign.mat), by=1)]
    
    ## compute likelihood
    # likelihood of case dna observations
    d.ll <- llfnDNA(theta = theta.d,
                    dcounts = dcounts,
                    log.ddepth = log.ddepth,
                    ddesign.mat = ddesign.mat)
    # likelihood of case rna observations
    r.ll <- llfnRNA(theta = c(theta.r, rctrlscale),
                    theta.d = theta.d,
                    rcounts = rcounts,
                    log.rdepth = log.rdepth,
                    ddesign.mat = ddesign.mat,
                    rdesign.mat = rbind(rdesign.mat, rctrldesign.mat) )
    
    return(d.ll + r.ll)
}

# this is for pre-fitting
cost.dna <- function(theta, theta.r,
                     llfnDNA, llfnRNA,
                     dcounts, rcounts,
                     log.ddepth, log.rdepth,
                     ddesign.mat, rdesign.mat) {
    
    ## compute likelihood
    d.ll <- llfnDNA(theta = theta,
                    dcounts = dcounts,
                    log.ddepth = log.ddepth,
                    ddesign.mat = ddesign.mat)
    r.ll <- llfnRNA(theta = theta.r,
                    theta.d = theta,
                    rcounts = rcounts,
                    log.rdepth = log.rdepth,
                    ddesign.mat = ddesign.mat,
                    rdesign.mat = rdesign.mat )
    
    return(d.ll + r.ll)
}
# this is for pre-fitting
cost.rna <- function(theta, theta.d,
                     llfnRNA,
                     rcounts,
                     log.rdepth,
                     ddesign.mat, rdesign.mat) {
    
    ## compute likelihood
    r.ll <- llfnRNA(theta = theta,
                    theta.d = theta,
                    rcounts = rcounts,
                    log.rdepth = log.rdepth,
                    ddesign.mat = ddesign.mat,
                    rdesign.mat = rdesign.mat )
    
    return(r.ll)
}


#' likelihood wrapper to compute terms of likelihood with non-zero derivative
#' with respect to case dna model and rna model if control enhancers are used
#' 
#' @param theta the vector of rna model parameters to evaluate likelihood for
#' (numeric, rna parameters)
#' @param theta.d.ctrl.prefit ctrl dna model parameters to condition likelihood on
#' (numeric, control enhancers x dna parameters)
#' @param llfnDNA cost function to compute dna likelihood terms on (function)
#' @param llfnRNA cost function to compute rna likelihood terms on (function)
#' @param dcounts the observed case DNA counts (integer, 1 x samples)
#' @param rcounts the observed RNA counts (integer, enhancers x samples)
#' @param log.ddepth dna library size correction vector, log scale (numeric, samples)
#' @param log.rdepth rna library size correction vector, log scale (numeric, samples)
#' @param ddesign.mat the dna model design matrix (logical, samples x dna parameters)
#' @param rdesign.mat the rna model design matrix (logical, samples x rna parameters)
#' @param rdesign.ctrl.mat the control rna model design matrix 
#' (logical, samples x rna parameters)
#' 
#' @return negative sum of log likelihood terms with non-zero derivative
#' with respect to case dna model and rna model if control enhancers are used
cost.dnarna.wctrl <- function(theta, theta.d.ctrl.prefit,
                              llfnDNA, llfnRNA,
                              dcounts, rcounts,
                              log.ddepth, log.rdepth, 
                              ddesign.mat, rdesign.mat, rdesign.ctrl.mat) {
    
    ## extract parameter vectors by model part
    # first parameter of DNA model is the variance(-link) parameter
    theta.d <- theta[seq(1, 1+NCOL(ddesign.mat), by=1)]
    theta.r <- theta[
        seq(1+NCOL(ddesign.mat)+1, 
            1+NCOL(ddesign.mat)+NCOL(rdesign.mat), by=1)]
    theta.r.ctrl <- theta[
        seq(1+NCOL(ddesign.mat)+NCOL(rdesign.mat)+1, 
            1+NCOL(ddesign.mat)+NCOL(rdesign.mat)+NCOL(rdesign.ctrl.mat), by=1)]
    
    ## compute likelihood
    # likelihood of case dna observations
    d.ll <- llfnDNA(theta = theta.d,
                    dcounts = dcounts,
                    log.ddepth = log.ddepth, 
                    ddesign.mat = ddesign.mat)
    # likelihood of case rna observations
    r.ll.case <- llfnRNA(
        theta = theta.r,
        theta.d = theta.d,
        rcounts = rcounts[1,,drop=FALSE],
        log.rdepth = log.rdepth,
        ddesign.mat = ddesign.mat,
        rdesign.mat = rdesign.mat )
    # likelihood of ctrl rna observations
    r.ll.ctrl <- llfnRNA(
        theta = c(theta.r, theta.r.ctrl),
        theta.d = theta.d.ctrl.prefit,
        rcounts = rcounts[-1,,drop=FALSE],
        log.rdepth = log.rdepth,
        ddesign.mat = ddesign.mat,
        rdesign.mat = cbind(rdesign.mat, rdesign.ctrl.mat) )
    
    return(d.ll + r.ll.case + r.ll.ctrl)
}

#' likelihood wrapper to compute terms of likelihood with non-zero derivative
#' with respect to case dna model if control enhancers are used
#' 
#' This is conditioned on RNA model.
#' 
#' @param theta the vector of rna model parameters to evaluate likelihood for
#' (numeric, rna parameters)
#' @param theta.r rna model parameters to condition likelihood on 
#' (numeric, rna parameters)
#' @param llfnDNA cost function to compute dna likelihood terms on (function)
#' @param llfnRNA cost function to compute rna likelihood terms on (function)
#' @param dcounts the observed case DNA counts (integer, 1 x samples)
#' @param rcounts the observed case RNA counts (integer, 1 x samples)
#' @param log.ddepth dna library size correction vector, log scale (numeric, samples)
#' @param log.rdepth rna library size correction vector, log scale (numeric, samples)
#' @param ddesign.mat the dna model design matrix (logical, samples x dna parameters)
#' @param rdesign.mat the rna model design matrix (logical, samples x rna parameters)
#' @param rdesign.ctrl.mat the control rna model design matrix 
#' (logical, samples x rna parameters)
#' 
#' @return negative sum of log likelihood terms with non-zero derivative
#' with respect to case dna model if control enhancers are used
cost.dna.wctrl <- function(theta, theta.r,
                           llfnDNA, llfnRNA,
                           dcounts, rcounts,
                           log.ddepth, log.rdepth,
                           ddesign.mat, rdesign.mat, rdesign.ctrl.mat) {
    
    ## extract parameter vectors by model part
    # first parameter of DNA model is the variance(-link) parameter
    theta.d <- theta[seq(1, 1+NCOL(ddesign.mat), by=1)]
    
    ## compute likelihood
    # likelihood of case dna observations
    d.ll <- llfnDNA(theta = theta.d,
                    dcounts = dcounts,
                    log.ddepth = log.ddepth, 
                    ddesign.mat = ddesign.mat)
    # likelihood of case rna observations
    r.ll <- llfnRNA(theta = theta.r,
                    theta.d = theta.d,
                    rcounts = rcounts,
                    log.rdepth = log.rdepth,
                    ddesign.mat = ddesign.mat,
                    rdesign.mat = rdesign.mat)
    
    return(d.ll + r.ll)
}

#' likelihood wrapper to compute terms of likelihood with non-zero derivative
#' with respect to rna model if control enhancers are used
#' 
#' This is conditioned on DNA models.
#' @param theta the vector of rna model parameters to evaluate likelihood for
#' (numeric, rna parameters)
#' @param theta.d dna model parameters to condition likelihood on 
#' (numeric, dna parameters)
#' @param theta.d.ctrl.prefit ctrl dna model parameters to condition likelihood on
#' (numeric, control enhancers x dna parameters)
#' @param llfnRNA cost function to compute likelihood on (function)
#' @param rcounts the observed RNA counts (integer, enhancers x samples)
#' @param log.rdepth rna library size correction vector, log scale (numeric, samples)
#' @param ddesign.mat the dna model design matrix (logical, samples x dna parameters)
#' @param rdesign.mat the rna model design matrix (logical, samples x rna parameters)
#' @param rdesign.ctrl.mat the control rna model design matrix 
#' (logical, samples x rna parameters)
#' 
#' @return negative sum of log likelihood terms with non-zero derivative
#' with respect to rna model if control enhancers are used
cost.rna.wctrl <- function(theta, theta.d, theta.d.ctrl.prefit,
                           llfnRNA,
                           rcounts,
                           log.rdepth, 
                           ddesign.mat, rdesign.mat, rdesign.ctrl.mat) {
    
    ## extract parameter vectors by model part
    theta.r <- theta[seq(1, NCOL(rdesign.mat)), by=1)]
    theta.r.ctrl <- theta[seq(NCOL(rdesign.mat)+1, 
                              NCOL(rdesign.mat)+NCOL(rdesign.ctrl.mat), by=1)]
    
    ## compute liklihood
    # likelihood of case rna observations
    r.ll.case <- llfnRNA(
        theta = theta.r,
        theta.d = theta.d,
        rcounts = rcounts[1,,drop=FALSE],
        log.rdepth = log.rdepth,
        ddesign.mat = ddesign.mat,
        rdesign.mat = rdesign.mat )
    # likelihood of control rna observations
    r.ll.ctrl <- llfnRNA(
        theta = c(theta.r, theta.r.ctrl),
        theta.d = theta.d.ctrl.prefit,
        rcounts = rcounts[-1,,drop=FALSE],
        log.rdepth = log.rdepth,
        ddesign.mat = ddesign.mat,
        rdesign.mat = cbind(rdesign.mat, rdesign.ctrl.mat) )
    
    return(r.ll.case + r.ll.ctrl)
}

