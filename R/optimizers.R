#' Fit dna and rna model to a given enhancer
#' 
#' Optim wrapper that performs numerical optimisation and extracts results.
#' Depending on the function chosen, a single optimisation or iterative 
#' estimation are performed.
#' 
#' @name fit.dnarna
#' @rdname fit.dnarna
#' 
#' @import stats
#' 
#' @aliases 
#' fit.dnarna.noctrlobs
#' fit.dnarna.wctrlobs.iter
#' 
#' @details dna model
#' The dna model is encoded as a vector where the first parameter
#' is the variance-link (ln: stdv, nb: dispersion) and the remaining parameters
#' are the coefficient of the linear model that consitutes the mean parameter,
#' as encoded in the model matrix.
#' 
#' @details rna model
#' The case rna model and the control rna interaction terms are encoded as 
#' seperate linear models so that they can be added up by design matrix 
#' concatenation (cbind) if a control enhancer is considered. Note that the 
#' first parameter may be a variance link parameter (ln.nb: nb: dispersion) 
#' which does not contribute to the linear model that constitutes the mean 
#' parameter. As all parameters only act on the mean parameter in the case of 
#' gamma.pois, the rna model matrix is augmented by a first column containing 
#' only padding 0s if such a variance parameter exists, so that the mean 
#' parameter can be computed via the same matrix multiplication irrespective of
#' whether the rna model contains such parameter that does not contribute to the
#' mean model. The code that throws away non-used coefficients before fitting
#' is blocked for this parameter. Note that the ctrl rna model is not allowed to
#' have and additional variance coefficient so that the testing is only on the 
#' first moment. This is guaranteed by the default that this non-use cofficient
#' is discarded. 
#' 
#' @param model noise model
#' @param dcounts the DNA counts
#' @param rcounts the RNA counts
#' @param ddepth dna library size correction vector (numeric, dna samples)
#' @param rdepth rna library size correction vector (numeric, rna samples)
#' @param ddesign.mat the dna model design matrix 
#' (logical, rna samples x dna parameters)
#' @param rdesign.mat the rna model design matrix 
#' (logical, rna samples x rna parameters)
#' @param d2rdesign.mat the transition matrix relating DNA estimates to RNA 
#' observations (logical, rna sample x dna parameters)
#' @param rctrlscale control-based correction scalers
#' @param rdesign.ctrl.mat the control rna model design matrix 
#' (logical, samples x rna parameters)
#' @param theta.d.ctrl.prefit ctrl dna model parameters to condition likelihood
#' on (numeric, control enhancers x dna parameters)
#' @param compute.hessian if TRUE (default), compute the Hessian matrix of the
#' coefficients to facilitate coefficient-based hypothesis testing
#'
#' @return a list with components:
#' \itemize{
#'     \item d.fitval: the fitted values of the DNA counts
#'     \item d.est: the estimated true DNA levels (corrected for library size)
#'     \item d.coef: the fitted mode parameters for the DNA counts
#'     \item d.se: the standard errors of the estimates of the DNA counts
#'     \item r.est: the fitted values of the RNA counts
#'     \item r.fitval: the fitted mode parameters for the RNA counts
#'     \item r.est: the estimated true RNA levels (corrected for library size)
#'     \item r.se: the standard errors of the estimates of the RNA counts
#'     \item ll: the log likelihood of the model
#' }
#' @noRd
NULL

#' @rdname fit.dnarna
#' @noRd
fit.dnarna.noctrlobs <- function(model, dcounts, rcounts,
                                 ddepth, rdepth, rctrlscale,
                                 ddesign.mat, rdesign.mat, d2rdesign.mat,
                                 rdesign.ctrl.mat, theta.d.ctrl.prefit,
                                 compute.hessian) {
    ## get likelihood function
    ll.funs <- get.ll.functions(model)
    
    ## filter invalid counts (NAs) from data and design
    valid.c.d <- (dcounts > 0) & !is.na(dcounts)
    valid.c.r <- (rcounts > 0) & !is.na(rcounts)
    
    dcounts.valid <- dcounts[valid.c.d]
    rcounts.valid <- rcounts[1,valid.c.r,drop=FALSE]
    log.ddepth.valid <- log(ddepth[valid.c.d])
    log.rdepth.valid <- log(rdepth[valid.c.r])
    
    ## clean design matrix from unused factors
    # valid.df <- apply(ddesign.mat[valid.c.d,,drop=FALSE], 2, 
    #                   function(x) !all(x==0))
    # valid.rf <- apply(rdesign.mat[valid.c.r,,drop=FALSE], 2,
    #                   function(x) !all(x==0))
    valid.df <- apply(ddesign.mat[valid.c.d,,drop=FALSE], 2,
                      function(x) sum(x!=0) > 1)
    if(!any(valid.df)) {
        valid.df <- apply(ddesign.mat[valid.c.d,,drop=FALSE], 2,
                          function(x) !all(x==0))
    }
    valid.rf <- apply(rdesign.mat[valid.c.r,,drop=FALSE], 2,
                      function(x) sum(x!=0) > 1)
    if(!any(valid.rf)) {
        valid.rf <- apply(rdesign.mat[valid.c.d,,drop=FALSE], 2,
                          function(x) !all(x==0))
    }
    ddmat.valid <- ddesign.mat[valid.c.d,valid.df,drop=FALSE]
    rdmat.valid <- rdesign.mat[valid.c.r,valid.rf,drop=FALSE]
    d2rmat.valid <- d2rdesign.mat[valid.c.r,valid.df,drop=FALSE]
    rdmat.ctrl.valid <- rdesign.ctrl.mat[valid.c.r,,drop=FALSE]
    
    ## Initialize parameter vector with a guess
    guess <- rep(0, 1 + NCOL(ddmat.valid) + NCOL(rdmat.valid))
    if(length(dcounts.valid) > 1) {
        guess[1] <- log(sd(dcounts.valid))
    }
    means <- as.vector(log((dcounts.valid %*% ddmat.valid) / colSums(ddmat.valid)))
    guess[2:(1 + NCOL(ddmat.valid))] <- (means - means[1])
    guess[2] <- means[1]
    
    ## optimize
    suppressWarnings(
        fit <- optim(
            par = guess,
            fn = cost.dnarna, 
            llfnDNA = ll.funs$dna, 
            llfnRNA = ll.funs$rna,
            dcounts = dcounts.valid, 
            rcounts = rcounts.valid,
            log.ddepth = log.ddepth.valid, 
            log.rdepth = log.rdepth.valid,
            rctrlscale = rctrlscale,
            ddesign.mat = ddmat.valid, 
            rdesign.mat = rdmat.valid,
            d2rdesign.mat = d2rmat.valid,
            rdesign.ctrl.mat = rdmat.ctrl.valid,
            hessian = compute.hessian,
            method = "L-BFGS-B", control = list(maxit=1000), 
            lower=-23, upper=23)
    )
    
    ## split parameters to the two parts of the model
    fit$par <- pmax(pmin(fit$par, 23), -23)
    d.par <- fit$par[seq(1, 1+NCOL(ddmat.valid))]
    r.par <- fit$par[c(1, seq(1+NCOL(ddmat.valid)+1,
                              1+NCOL(ddmat.valid)+NCOL(rdmat.valid)))]
    
    d.coef <- c(d.par[1], rep(NA, NCOL(ddesign.mat)))
    d.coef[1 + which(valid.df)] <- d.par[-1]
    d.df <- length(d.par)
    
    r.coef <- c(r.par[1], rep(NA, NCOL(rdesign.mat)))
    r.coef[1 + which(valid.rf)] <- r.par[-1]
    r.df <- length(r.par)
    
    ## standard error of the estimates
    se <- NULL
    d.se <- NULL
    r.se <- NULL
    if (compute.hessian) {
        se.comp <- tryCatch({
            se <- sqrt(diag(solve(fit$hessian)))
        }, error = function(e) return(e))
        
        if(!inherits(se.comp, "error")) {
            d.se <- c(se[1], rep(NA, NCOL(ddesign.mat)))
            d.se[1 + which(valid.df)] <- se[seq(2, 1+NCOL(ddmat.valid))]
            r.se <- c(se[1], rep(NA, NCOL(rdesign.mat)))
            r.se[1 + which(valid.rf)] <- se[seq(1+NCOL(ddmat.valid)+1,
                                        1+NCOL(ddmat.valid)+NCOL(rdmat.valid))]
        } else {
            message(se.comp$message)
        }
    }
    
    return(list(
        d.coef = d.coef, d.se = d.se, d.df = d.df,
        r.coef = r.coef, r.se = r.se, r.df = r.df,
        r.ctrl.coef = NULL, r.ctrl.se = NULL, r.ctrl.df = 0,
        converged = fit$convergence,
        ll = -fit$value))
}

#' fit the control-based joint model. Seperate DNA model per enhancer, joint RNA
#' model.
#' 
#' @inheritParams fit.dnarna
#' @param print.progress print a progress report on the fitting (default:TRUE),
#' recommended since the joint model fitting is computationally intensive
#' @param BPPARAM the parallelization backend to use
#' @return a list with components:
#' \itemize{
#'     \item d.fitval: the fitted values of the DNA counts
#'     \item d.est: the estimated true DNA levels (corrected for library size)
#'     \item d.coef: the fitted mode parameters for the DNA counts
#'     \item d.se: the standard errors of the estimates of the DNA counts
#'     \item r.est: the fitted values of the RNA counts
#'     \item r.fitval: the fitted mode parameters for the RNA counts
#'     \item r.est: the estimated true RNA levels (corrected for library size)
#'     \item r.se: the standard errors of the estimates of the RNA counts
#'     \item ll: the log likelihood of the model
#' }
#' @noRd
fit.dnarna.onlyctrl.iter <- function(model, dcounts, rcounts,
                                     ddepth, rdepth,
                                     ddesign.mat, rdesign.mat, d2rdesign.mat,
                                     BPPARAM, print.progress=TRUE) {
    
    ## get cost function
    ll.funs <- get.ll.functions(model)
    
    ## filter invalid counts (0, NAs) from data and design
    valid.c.d <- (dcounts > 0) & !is.na(dcounts)
    valid.c.r.agg <- apply((rcounts > 0) & !is.na(rcounts), 2, any) 
    log.ddepth <- log(ddepth)
    log.rdepth <- log(rdepth)
    
    ## clean design matrix from unused factors
    valid.rf <- apply(rdesign.mat[valid.c.r.agg,,drop=FALSE], 2, 
                      function(x) sum(x!=0) > 1)
    if(!any(valid.rf)) {
        valid.rf <- apply(rdesign.mat[valid.c.r.agg,,drop=FALSE], 2,
                          function(x) !all(x==0))
    }
    
    ## Iterative parameter estimation: coordinate ascent
    # Iterate DNA and RNA model estimation
    # Initialize DNA model parameter vector with a guess
    d.par <- matrix(0, nrow=NROW(dcounts), ncol=1+NCOL(ddesign.mat))
    r.par <- rep(0, 1+NCOL(rdesign.mat))
    
    llold <- 0
    llnew <- 0
    iter <- 1
    converged <- TRUE
    RELTOL <- 10^(-6)
    MAXITER <- 1000
    
    
    while((llnew > llold-llold*RELTOL | iter <= 2) & iter < MAXITER) {
        ## estimate dna model for each control enhancer
        dfits <- bplapply(seq_len(NROW(dcounts)), function(i) {
            valid.df <- apply(ddesign.mat[valid.c.d[i,],,drop=FALSE], 2, 
                              function(x) !all(x==0))
            d.par <- rep(0, sum(valid.df) + 1)
            suppressWarnings(
            fit <- optim(
                par = d.par, 
                fn = cost.dna, 
                theta.r = r.par,
                llfnDNA = ll.funs$dna, 
                llfnRNA = ll.funs$rna,
                dcounts = t(dcounts[i,valid.c.d[i,],drop=FALSE]),
                rcounts = t(rcounts[i,valid.c.r.agg,drop=FALSE]),
                log.ddepth = log.ddepth[valid.c.d[i,]], 
                log.rdepth = log.rdepth[valid.c.r.agg],
                ddesign.mat = ddesign.mat[valid.c.d[i,],valid.df,drop=FALSE], 
                rdesign.mat = rdesign.mat[valid.c.r.agg,valid.rf,drop=FALSE], 
                d2rdesign.mat = d2rdesign.mat[valid.c.r.agg,
                                              valid.df,drop=FALSE],
                hessian = FALSE, 
                method = "L-BFGS-B", control = list(maxit=1000), 
                lower=-23, upper=23)
            )
            fit$par <- pmax(pmin(fit$par, 23), -23)
            return(fit)
        }, BPPARAM = BPPARAM)
        
        d.par <- matrix(0, nrow=NROW(dcounts), ncol=NCOL(ddesign.mat)+1)
        d.par[,1] <- vapply(dfits, function(x) x$par[1], 0.0)
        for(i in seq_len(NROW(dcounts))){
            valid.df <- apply(ddesign.mat[valid.c.d[i,],,drop=FALSE], 2, 
                              function(x) !all(x==0))
            d.par[i, 1 + which(valid.df)] <- dfits[[i]]$par[-1]
        }
        
        ## estimate rna model conditioned on dna model
        suppressWarnings(
            rfit <- optim(
                par = r.par, 
                fn = cost.rna, 
                theta.d = t(d.par),
                llfnRNA = ll.funs$rna, 
                rcounts = t(rcounts[,valid.c.r.agg,drop=FALSE]),
                log.rdepth = log.rdepth[valid.c.r.agg],
                d2rdesign.mat = d2rdesign.mat[valid.c.r.agg,,drop=FALSE], 
                rdesign.mat = rdesign.mat[valid.c.r.agg,valid.rf,drop=FALSE], 
                hessian = FALSE, 
                method = "L-BFGS-B", control = list(maxit=1000), 
                lower=-23, upper=23)
        )
        rfit$par <- pmax(pmin(rfit$par, 23), -23)
        r.par <- rfit$par
        names(r.par) <- NULL
        
        ## update iteration convergence reporters
        llold <- llnew
        llnew <- -sum(rfit$value + vapply(dfits, function(x) x$value, 0.0))
        iter <- iter + 1
        if(iter == MAXITER & llnew > llold-llold*RELTOL) {
            converged <- FALSE
        }
        if(print.progress) {
            message("iter:", iter, "\tlog-likelihood:", llnew)
        }
    }
    
    d.coef <- d.par
    d.df <- vapply(dfits, function(x) length(x$fit), 0)
    
    r.coef <- rep(NA, NCOL(rdesign.mat))
    r.coef[valid.rf] <- r.par[-1]
    r.df <- length(r.par) - 1
    return(lapply(seq_len(NROW(d.coef)), function(i){
        list(
            d.coef = d.coef[i,], d.se = NULL,      d.df = d.df[i],
            r.coef = r.coef,     r.se = NULL,      r.df = r.df,
            r.ctrl.coef = NULL,  r.ctrl.se = NULL, r.ctrl.df = 0,
            converged = converged,
            ll = NA)
    }))
}
