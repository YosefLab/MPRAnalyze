#' Get appropriate likelihood function accourding to the selected model
#' 
#' @param model the distributional model selected
#' 
#' @return a list with two items: dna, rna. Each item is the corresponding 
#' likelihod function
get.ll.functions <- function(model) {
    if(model == "gamma.pois"){
        llfnDNA <- ll.dna.gamma.pois
        llfnRNA <- ll.rna.gamma.pois
    } else if(model == "ln.nb"){
        llfnDNA <- ll.dna.ln
        llfnRNA <- ll.rna.ln.nb
    } else if(model == "ln.ln"){
        llfnDNA <- ll.dna.ln
        llfnRNA <- ll.rna.ln.ln
    } else {
        stop("model ", model, " not supported")
    }
    return(list(dna=llfnDNA, rna=llfnRNA))
}

#' likelihood of dna observations
#' 
#' @name ll.dna
#' @rdname ll.dna
#' 
#' @aliases 
#' ll.dna.gamma.pois
#' ll.dna.ln.nb
#' ll.dna.ln.ln
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
ll.dna.gamma.pois <- function(theta, dcounts, log.ddepth, ddesign.mat) {
    
    log.rate.est <- matrix(rep(-(ddesign.mat %*% theta[-1] + log.ddepth), 
                            NCOL(dcounts)), ncol = NCOL(dcounts))
    ## compute likelihood
    # shape_gamma = alpha
    # rate_gamma = beta = 1/model_dna
    ll <- sum(dgamma(x = dcounts,
                    shape = exp(theta[1]), # alpha
                    rate = exp(log.rate.est),
                    log = TRUE))
    return(-ll)
}

#' @rdname ll.dna
ll.dna.ln <- function(theta, dcounts, log.ddepth, ddesign.mat) {
    
    log.d.est <- matrix(rep(ddesign.mat %*% theta[-1] + log.ddepth, 
                            NCOL(dcounts)), ncol=NCOL(dcounts))
    
    ## compute likelihood
    ll <- sum(dlnorm(x = dcounts,
                    meanlog = log.d.est,
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
#' ll.rna.ln.ln
#'
#' @param theta the vector of rna model parameters to evaluate likelihood for
#' (numeric, rna parameters)
#' @param theta.d dna model parameters to condition likelihood on
#' @param rcounts the observed RNA counts (integer, samples x enhancers)
#' @param log.rdepth rna library size correction vector, log scale (numeric, samples)
#' @param d2rdesign.mat the trasitional matrix to distribute and match DNA 
#' estimates to RNA observations (logical, rna samples x dna parameters)
#' @param rdesign.mat the rna model design matrix (logical, rna samples x rna parameters)
#'
#' @return negative log likelihood of rna observations under the specified model
NULL

#' @rdname ll.rna
ll.rna.gamma.pois <- function(theta, theta.d, 
                            rcounts, log.rdepth,
                            d2rdesign.mat, rdesign.mat) {
    
    alpha.mat <- matrix(rep(theta.d[1,], each=NROW(d2rdesign.mat)), 
                        nrow=NROW(d2rdesign.mat))
    log.d.est <- (d2rdesign.mat %*% theta.d[-1,,drop=FALSE]) + alpha.mat
    log.r.est <- log.d.est + rep((rdesign.mat %*% theta[-1]) + log.rdepth, 
                                NCOL(log.d.est))
    
    ## compute likelihood
    # mu_NB = alpha / beta * rna_model
    #       = alpha * dna_model * rna_model
    # size_alpha = alpha
    ll <- sum(dnbinom(x = rcounts,
                    size = exp(alpha.mat),
                    mu = exp(log.r.est),
                    log = TRUE))
    return(-ll)
}

#' @rdname ll.rna
ll.rna.ln.nb <- function(theta, theta.d, 
                        rcounts, log.rdepth,
                        d2rdesign.mat, rdesign.mat) {
    
    log.d.est <- d2rdesign.mat %*% theta.d[-1,]
    log.r.est <- log.d.est + rep(((rdesign.mat %*% theta[-1]) + log.rdepth), 
                                NCOL(log.d.est))
    
    ## compute likelihood
    ll <- sum(dnbinom(x = rcounts,
                    size = exp(exp(theta[1])),
                    mu = exp(log.r.est),
                    log = TRUE))
    
    return(-ll)
}

#' @rdname ll.rna
ll.rna.ln.ln <- function(theta, theta.d, 
                        rcounts, log.rdepth,
                        d2rdesign.mat, rdesign.mat) {
    
    log.d.est <- d2rdesign.mat %*% theta.d[-1,]
    log.r.est <- log.d.est + rep(((rdesign.mat %*% theta[-1]) + log.rdepth), 
                                NCOL(log.d.est))
    
    ## compute likelihood
    ll <- sum(dlnorm(x = rcounts,
                    meanlog = log.r.est,
                    sdlog = exp(theta[1]),
                    log = TRUE))
    
    return(-ll)
}