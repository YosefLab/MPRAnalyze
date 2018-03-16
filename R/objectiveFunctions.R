#' likelihood wrapper to compute terms of likelihood with non-zero derivative
#' with respect to case dna model and rna model 
#' 
#' nomenclature cost.[model components to be estimated].[if ctrl obs are included in fitting]
#' 
#' @name cost.model
#' @rdname cost.model
#' 
#' @param theta vector of model parameters to evaluate likelihood for
#' (numeric, parameters)
#' @param theta.d vector of dna model parameters to condition likelihood on
#' (numeric, dna parameters)
#' @param theta.r vector of rna model parameters to condition likelihood on
#' (numeric, rna parameters)
#' @param llfnDNA cost function to compute dna likelihood terms on (function)
#' @param llfnRNA cost function to compute rna likelihood terms on (function)
#' @param dcounts the observed case DNA counts (integer, enhancers x dna samples)
#' @param rcounts the observed RNA counts (integer, enhancers x rna samples)
#' @param log.ddepth dna library size correction vector, log scale (numeric, dna samples)
#' @param log.rdepth rna library size correction vector, log scale (numeric, rna samples)
#' @param ddesign.mat the dna model design matrix (logical, dna samples x dna parameters)
#' @param rdesign.mat the rna model design matrix (logical, rna samples x rna parameters)
#' @param d2rdesign.mat the transitional design matrix, distributiong the DNA estimates to the RNA observations 
#' @param rdesign.ctrl.mat the control rna model design matrix 
#' (logical, samples x rna parameters)
NULL

#' likelihood wrapper to compute terms of likelihood with non-zero derivative
#' with respect to case dna model and rna model without control enhancer observations
#' 
#' nomenclature cost.[model components to be estimated]
#' 
#' Each of these function can only be called to optimise a single (enhancer) model:
#' cost.dnarna and cost.dna can only be used on observations from one enhancer
#' cost.rna can be used on observation from multiple enhancers under the same rna model
#' 
#' @name cost.model.noctrl
#' @rdname cost.model.noctrl
#' 
#' @aliases 
#' cost.dnarna
#' cost.dna
#' cost.rna
#' 
#' @inheritParams cost.model
#' @param rctrlscale vector of scaling parameters for rna model that correct for 
#' variation between conditions pre-fit on control enhancers 
#' (numeric, ctrl rna parameters)
#' 
#' @return negative sum of log likelihood terms with non-zero derivative
#' with respect to case model 
NULL

#' @rdname cost.model.noctrl
cost.dnarna <- function(theta, dcounts, rcounts,
                        llfnDNA, llfnRNA,
                        log.ddepth, log.rdepth, rctrlscale=NULL,
                        ddesign.mat, rdesign.mat, d2rdesign.mat,
                        rdesign.ctrl.mat=NULL) {
    ## extract parameter vectors by model part
    # first parameter of DNA model is the variance(-link) parameter
    theta.d <- theta[seq(1, 1+NCOL(ddesign.mat), by=1)]
    theta.r <- theta[c(1, seq(1+NCOL(ddesign.mat)+1, 
                            1+NCOL(ddesign.mat)+NCOL(rdesign.mat), by=1))]
    ## compute likelihood
    # likelihood of case dna observations
    d.ll <- llfnDNA(theta = theta.d,
                    dcounts = dcounts,
                    log.ddepth = log.ddepth,
                    ddesign.mat = ddesign.mat)
    # likelihood of case rna observations
    r.ll <- llfnRNA(theta = c(theta.r, rctrlscale),
                    theta.d = matrix(theta.d),
                    rcounts = rcounts,
                    log.rdepth = log.rdepth,
                    d2rdesign.mat = d2rdesign.mat,
                    rdesign.mat = cbind(rdesign.mat, rdesign.ctrl.mat)) 

    return(d.ll + r.ll)
}

#' @rdname cost.model.noctrl
cost.dna <- function(theta, theta.r,
                     llfnDNA, llfnRNA,
                     dcounts, rcounts,
                     log.ddepth, log.rdepth, rctrlscale=NULL,
                     ddesign.mat, rdesign.mat, d2rdesign.mat, 
                     rdesign.ctrl.mat=NULL) {
    
    # theta <- pmax(pmin(theta, 23), -23)
    ## compute likelihood
    d.ll <- llfnDNA(theta = theta,
                    dcounts = dcounts,
                    log.ddepth = log.ddepth,
                    ddesign.mat = ddesign.mat)
    
    r.ll <- llfnRNA(theta = c(theta.r, rctrlscale),
                    theta.d = matrix(theta),
                    rcounts = rcounts,
                    log.rdepth = log.rdepth,
                    d2rdesign.mat = d2rdesign.mat,
                    rdesign.mat = cbind(rdesign.mat, rdesign.ctrl.mat) )
    return(d.ll + r.ll)
}

#' @rdname cost.model.noctrl
cost.rna <- function(theta, theta.d, llfnRNA, rcounts,
                     log.rdepth, rctrlscale=NULL,
                     d2rdesign.mat, rdesign.mat,
                     rdesign.ctrl.mat=NULL) {
    
    # theta <- pmax(pmin(theta, 23), -23)
    ## compute likelihood
    r.ll <- llfnRNA(theta = c(theta, rctrlscale),
                    theta.d = theta.d,
                    rcounts = rcounts,
                    log.rdepth = log.rdepth,
                    d2rdesign.mat = d2rdesign.mat,
                    rdesign.mat = cbind(rdesign.mat, rdesign.ctrl.mat))
    return(r.ll)
}


#' likelihood wrapper to compute terms of likelihood with non-zero derivative
#' with respect to case dna model and rna model with control enhancer observations
#' 
#' nomenclature cost.[model components to be estimated].wctrl
#' 
#' @name cost.model.wctrl
#' @rdname cost.model.wctrl
#' 
#' @aliases 
#' cost.dna.wctrl
#' cost.rna.wctrl
#' 
#' @inheritParams cost.model
#' @param theta.d.ctrl.prefit prefit control DNA model parameters
#' (numeric, control enhancers x dna model parameters)
#' @param d2rdesign.ctrl.mat model design to distribution DNA pre-fit estimates
#' to the RNA observations.
#' (logical, rna samples x dna parameters)
#' 
#' @return negative sum of log likelihood terms with non-zero derivative
#' with respect to case model 
NULL

#' @rdname cost.model.wctrl
cost.dna.wctrl <- function(theta, theta.r, llfnDNA, llfnRNA,
                           dcounts, rcounts, log.ddepth, log.rdepth,
                           ddesign.mat, rdesign.mat, d2rdesign.mat) {
    
    # theta <- pmax(pmin(theta, 23), -23)
    ## likelihood of case dna observations
    d.ll <- llfnDNA(theta = theta,
                    dcounts = dcounts,
                    log.ddepth = log.ddepth, 
                    ddesign.mat = ddesign.mat)
    
    ## likelihood of case rna observations
    r.ll <- llfnRNA(theta = theta.r,
                    theta.d = matrix(theta),
                    rcounts = rcounts[1,,drop=FALSE],
                    log.rdepth = log.rdepth,
                    d2rdesign.mat = d2rdesign.mat,
                    rdesign.mat = rdesign.mat)
    
    return(d.ll + r.ll)
}

#' @rdname cost.model.wctrl
cost.rna.wctrl <- function(theta, theta.d, theta.d.ctrl.prefit,
                           llfnRNA, 
                           rcounts,
                           log.ddepth, log.rdepth, 
                           d2rdesign.mat, rdesign.mat,
                           d2rdesign.ctrl.mat, rdesign.ctrl.mat) {
    
    # theta <- pmax(pmin(theta, 23), -23)
    ## extract parameter vectors by model part
    theta.r <- theta[seq(1, 1+NCOL(rdesign.mat), by=1)]
    if(!is.null(rdesign.ctrl.mat)) {
        theta.r.ctrl <- theta[seq(1+NCOL(rdesign.mat)+1, 
                                  1+NCOL(rdesign.mat)+NCOL(rdesign.ctrl.mat), by=1),
                              drop=FALSE]
    } else {
        theta.r.ctrl <- NULL
    }
    
    ## compute liklihood
    # likelihood of case rna observations
    r.ll.case <- llfnRNA(
        theta = theta.r,
        theta.d = matrix(theta.d),
        rcounts = rcounts[1,,drop=FALSE],
        log.rdepth = log.rdepth,
        d2rdesign.mat = d2rdesign.mat,
        rdesign.mat = rdesign.mat)
    
    # likelihood of control rna observations
    r.ll.ctrl <- llfnRNA(
        theta = c(theta.r, theta.r.ctrl),
        theta.d = theta.d.ctrl.prefit,
        rcounts = rcounts[-1,,drop=FALSE],
        log.rdepth = log.rdepth,
        d2rdesign.mat = d2rdesign.ctrl.mat,
        rdesign.mat = cbind(rdesign.mat, rdesign.ctrl.mat) )
    
    return(r.ll.case + r.ll.ctrl)
}

