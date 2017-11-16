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
#' @param dcounts the observed case DNA counts (integer, enhancers x samples)
#' @param rcounts the observed RNA counts (integer, enhancers x samples)
#' @param log.ddepth dna library size correction vector, log scale (numeric, samples)
#' @param log.rdepth rna library size correction vector, log scale (numeric, samples)
#' @param ddesign.mat the dna model design matrix (logical, samples x dna parameters)
#' @param rdesign.mat the rna model design matrix (logical, samples x rna parameters)
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
cost.dnarna <- function(theta, theta.d=NULL, theta.r=NULL,
                        dcounts, rcounts,
                        llfnDNA, llfnRNA,
                        log.ddepth, log.rdepth, rctrlscale=NULL,
                        ddesign.mat, rdesign.mat, rctrldesign.mat=NULL) {
    
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
                    theta.d = t(matrix(theta.d)),
                    rcounts = rcounts,
                    log.rdepth = log.rdepth,
                    ddesign.mat = ddesign.mat,
                    rdesign.mat = rbind(rdesign.mat, rctrldesign.mat) )
    
    return(d.ll + r.ll)
}

#' @rdname cost.model.noctrl
cost.dna <- function(theta, theta.d=NULL, theta.r,
                     llfnDNA, llfnRNA,
                     dcounts, rcounts,
                     log.ddepth, log.rdepth, rctrlscale=NULL,
                     ddesign.mat, rdesign.mat, rdesign.ctrl.mat=NULL) {
    
    ## compute likelihood
    d.ll <- llfnDNA(theta = theta,
                    dcounts = dcounts,
                    log.ddepth = log.ddepth,
                    ddesign.mat = ddesign.mat)
    r.ll <- llfnRNA(theta = c(theta.r, rctrlscale),
                    theta.d = t(matrix(theta)),
                    rcounts = rcounts,
                    log.rdepth = log.rdepth,
                    ddesign.mat = ddesign.mat,
                    rdesign.mat = cbind(rdesign.mat, rdesign.ctrl.mat) )
    
    return(d.ll + r.ll)
}

#' @rdname cost.model.noctrl
cost.rna <- function(theta, theta.d, theta.r=NULL, 
                     llfnDNA=NULL, llfnRNA,
                     dcounts, rcounts,
                     log.ddepth=NULL, log.rdepth, rctrlscale=NULL,
                     ddesign.mat, rdesign.mat, rdesign.ctrl.mat=NULL) {
    
    ## compute likelihood
    r.ll <- llfnRNA(theta = c(theta, rctrlscale),
                    theta.d = theta.d,
                    rcounts = rcounts,
                    log.rdepth = log.rdepth,
                    ddesign.mat = ddesign.mat,
                    rdesign.mat = cbind(rdesign.mat, rdesign.ctrl.mat) )
    
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
#' cost.dnarna
#' cost.dna
#' cost.rna
#' 
#' @inheritParams cost.model
#' @param theta.d.ctrl.prefit prefit control DNA model parameters
#' (numeric, control enhancers x dna model parameters)
#' 
#' @return negative sum of log likelihood terms with non-zero derivative
#' with respect to case model 
NULL

#' @rdname cost.model.wctrl
cost.dnarna.wctrl <- function(theta, theta.d=NULL, theta.r=NULL, theta.d.ctrl.prefit,
                              llfnDNA, llfnRNA,
                              dcounts, rcounts,
                              log.ddepth, log.rdepth, 
                              ddesign.mat, rdesign.mat, rdesign.ctrl.mat) {
    
    ## extract parameter vectors by model part
    # first parameter of DNA model is the variance(-link) parameter
    theta.d <- theta[seq(1, 1+NCOL(ddesign.mat), by=1)]
    theta.r <- theta[
        1,seq(1+NCOL(ddesign.mat)+1, 
              1+NCOL(ddesign.mat)+NCOL(rdesign.mat), by=1)]
    theta.r.ctrl <- theta[
        1,seq(1+NCOL(ddesign.mat)+NCOL(rdesign.mat)+1, 
              1+NCOL(ddesign.mat)+NCOL(rdesign.mat)+NCOL(rdesign.ctrl.mat), by=1)]
    
    ## compute likelihood
    # likelihood of case dna observations
    d.ll <- llfnDNA(theta = theta.d,
                    dcounts = dcounts,
                    log.ddepth = log.ddepth, 
                    ddesign.mat = ddesign.mat)
    # likelihood of case rna observations
    r.ll.case <- llfnRNA(theta = theta.r,
                         theta.d = theta.d,
                         rcounts = rcounts[1,],
                         log.rdepth = log.rdepth,
                         ddesign.mat = ddesign.mat,
                         rdesign.mat = rdesign.mat )
    # likelihood of ctrl rna observations
    r.ll.ctrl <- llfnRNA(theta = c(theta.r, theta.r.ctrl),
                         theta.d = theta.d.ctrl.prefit,
                         rcounts = rcounts[-1,],
                         log.rdepth = log.rdepth,
                         ddesign.mat = ddesign.mat,
                         rdesign.mat = cbind(rdesign.mat, rdesign.ctrl.mat) )
    
    return(d.ll + r.ll.case + r.ll.ctrl)
}

#' @rdname cost.model.wctrl
cost.dna.wctrl <- function(theta, theta.d=NULL, theta.r,  theta.d.ctrl.prefit=NULL,
                           llfnDNA, llfnRNA,
                           dcounts, rcounts,
                           log.ddepth, log.rdepth,
                           ddesign.mat, rdesign.mat, rdesign.ctrl.mat) {
    
    ## extract parameter vectors by model part
    # first parameter of DNA model is the variance(-link) parameter
    theta.d <- t(as.matrix(theta))#[seq(1, 1+NCOL(ddesign.mat), by=1)]
    
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

#' @rdname cost.model.wctrl
cost.rna.wctrl <- function(theta, theta.d, theta.r=NULL, theta.d.ctrl.prefit,
                           llfnDNA=NULL, llfnRNA,
                           dcounts=NULL, rcounts,
                           log.ddepth, log.rdepth, 
                           ddesign.mat, rdesign.mat, rdesign.ctrl.mat) {
    
    ## extract parameter vectors by model part
    theta.r <- theta[seq(1, NCOL(rdesign.mat), by=1)]
    theta.r.ctrl <- theta[seq(NCOL(rdesign.mat)+1, 
                              NCOL(rdesign.mat)+NCOL(rdesign.ctrl.mat), by=1),
                          drop=FALSE]
    
    ## compute liklihood
    # likelihood of case rna observations
    r.ll.case <- llfnRNA(
        theta = theta.r,
        theta.d = theta.d,
        rcounts = rcounts[1,],
        log.rdepth = log.rdepth,
        ddesign.mat = ddesign.mat,
        rdesign.mat = rdesign.mat )
    # likelihood of control rna observations
    r.ll.ctrl <- llfnRNA(
        theta = c(theta.r, theta.r.ctrl),
        theta.d = theta.d.ctrl.prefit,
        rcounts = rcounts[-1,],
        log.rdepth = log.rdepth,
        ddesign.mat = ddesign.mat,
        rdesign.mat = cbind(rdesign.mat, rdesign.ctrl.mat) )
    
    return(r.ll.case + r.ll.ctrl)
}

