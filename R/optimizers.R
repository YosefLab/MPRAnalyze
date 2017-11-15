#' fit a log normal GLM to the DNA count data
#'
#' @param dcounts the DNA counts (integer, N)
#' @param depth library size correction factors (numeric, N)
#' @param design.mat the design matrix (logical, N x K)
#' @param compute.hessian if TRUE (default), compute the Hessian matrix of the
#' coefficients to facilitate coefficient-based hypothesis testing
#'
#' @return a list with components:
#' \itemize {
#'     \item fitval: the fitted values
#'     \item est: estimate of DNA levels. Differs from the fitval since this
#'     does not account for library depth. When using DNA estimates for fitting
#'     RNA models, this value should be used.
#'     \item coef: the fitted coefficients for the parameters
#'     \item se: the standard errors of the coefficients
#'     \item ll: the log likelihood of the model
#' }
fit.lnDNA <- function(dcounts, depth, design.mat, compute.hessian=TRUE) {
    ##TODO: lm does a faster and just as accurate job if depth factors are not
    ##a consideration. If we can incorporate them or remove the need to account
    ##for them, using lm will be preferable.
    
    ## filter invalid counts (NAs) from data and design
    valid.c <- (!is.na(dcounts) & dcounts != 0)
    dcounts.valid <- dcounts[valid.c]
    logdepth.valid <- log(depth[valid.c])
    
    ## clean design matrix from unused factors: note that these should be
    valid.f <- apply(design.mat[valid.c,], 2, function(x) !all(x==0))
    dmat.valid <- design.mat[valid.c,valid.f]
    
    ## Initialize parameter vector with a guess
    guess <- rep(1, NCOL(dmat.valid) + 1)
    
    ## fit the model
    fit <- optim(par = guess,
                 fn = ll.lnDNA,
                 dcounts = dcounts.valid,
                 logdepth = logdepth.valid,
                 design.mat = dmat.valid,
                 method = "BFGS",
                 hessian = compute.hessian)
    
    ## extract parameters
    est <- rep(NA, length(dcounts))
    est[valid.c] <- exp(dmat.valid %*% fit$par[-1])
    
    fv <- est * depth
    
    coef <- rep(NA, 1 + NCOL(design.mat))
    coef[c(1, which(valid.f) + 1)] <- fit$par
    names(coef) <- c("dispersion", colnames(design.mat))
    
    if(compute.hessian) {
        se <- rep(NA, 1 + NCOL(design.mat))
        se[c(1, 1 + which(valid.f))] <- sqrt(diag(solve(fit$hessian)))
        names(se) <- c("dispersion", colnames(design.mat))
    } else {
        se <- NULL
    }
    
    return(list(fitval = fv, est = est, coef = coef, se = se, ll = -fit$value))
}

#' fit a Gamma GLM to the DNA count data
#'
#' @param dcounts the DNA counts (integer, N)
#' @param depth library size correction factors (numeric, N)
#' @param design.mat the design matrix (logical, N x K)
#' @param compute.hessian if TRUE (default), compute the Hessian matrix of the
#' coefficients to facilitate coefficient-based hypothesis testing
#'
#' @return a list with components:
#' \itemize {
#'     \item fitval: the fitted values
#'     \item est: estimates of actual DNA levels, does not account for library
#'     depth factors. Should be used in downstream estimation instead of fitval
#'     \item coef: the fitted coefficients for the parameters
#'     \item se: the standard errors of the coefficients
#'     \item ll: the log likelihood of the model
#' }
fit.gammaDNA <- function(dcounts, depth, design.mat, compute.hessian=TRUE) {
    ## filter invalid counts (NAs) from data and design
    valid.c <- (!is.na(dcounts) & dcounts != 0)
    dcounts.valid <- dcounts[valid.c]
    logdepth.valid <- log(depth[valid.c])
    
    ## clean design matrix from unused factors: note that these should be
    valid.f <- apply(design.mat[valid.c,], 2, function(x) !all(x==0))
    dmat.valid <- design.mat[valid.c,valid.f]
    
    ## Initialize parameter vector with a guess
    guess <- rep(1, NCOL(dmat.valid) + 1)
    
    ## fit the model
    fit <- optim(par = guess,
                 fn = ll.gammaDNA,
                 dcounts = dcounts.valid,
                 logdepth = logdepth.valid,
                 design.mat = dmat.valid,
                 method = "BFGS",
                 hessian = compute.hessian)
    
    ## extract parameters
    est <- rep(NA, length(dcounts))
    est[valid.c] <- exp((dmat.valid %*% fit$par[-1]) + fit$par[1])
    
    fv <- est * depth
    
    coef <- rep(NA, 1 + NCOL(design.mat))
    coef[c(1, which(valid.f) + 1)] <- fit$par
    ##TODO: name the parameters!
    
    if(compute.hessian) {
        se <- rep(NA, 1 + NCOL(design.mat))
        se[c(1, which(valid.f) + 1)] <- sqrt(diag(solve(fit$hessian)))
        ##TODO: name the parameters!
    } else {
        se <- NULL
    }
    
    return(list(fitval = fv, est = est, coef = coef, se = se, ll = -fit$value))
}

#' fit a negative binomial model to the RNA counts using a point estimator for
#' the DNA counts
#'
#' @param rcounts RNA counts
#' @param depth library size correction factors
#' @param d.est point estimates of the DNA base counts
#' @param design.mat the RNA design matrix
#' @param compute.hessian if TRUE (default), compute the Hessian matrix of the
#' coefficients to facilitate coefficient-based hypothesis testing
#'
#' @return a list with components:
#' \itemize {
#'     \item fitval: the fitten values
#'     \item est: the estimated RNA levels (controled for library size)
#'     \item coed: the fitted coefficients for the parameters
#'     \item se: the standard errors of the coefficients
#'     \item ll: the log likelihood of the model
#' }
fit.nbRNA <- function(rcounts, depth, d.est, design.mat, compute.hessian=TRUE) {
    ## filter invalid counts (NAs) from data and design
    valid.c <- !is.na(rcounts) & !is.na(d.est)
    rcounts.valid <- rcounts[valid.c]
    logdepth.valid <- log(depth[valid.c])
    logd.est.valid <- log(d.est[valid.c])
    
    ## clean design matrix from unused factors: note that these should be
    valid.f <- apply(design.mat[valid.c,], 2, function(x) !all(x==0))
    dmat.valid <- design.mat[valid.c,valid.f]
    
    ## Initialize parameter vector with a guess
    guess <- rep(0, NCOL(dmat.valid) + 1)
    
    ## fit the model
    fit <- optim(par = guess,
                 fn = ll.nbRNA,
                 rcounts = rcounts.valid,
                 logdepth = logdepth.valid,
                 log.dest = logd.est.valid,
                 design.mat = dmat.valid,
                 method = "BFGS",
                 hessian = compute.hessian)
    
    ## extract parameters
    est <- rep(NA, length(dcounts))
    est[valid.c] <- exp((dmat.valid %*% fit$par[-1]) + logd.est.valid)
    
    fv <- est * depth
    
    coef <- rep(NA, 1 + NCOL(design.mat))
    coef[c(1, which(valid.f) + 1)] <- fit$par
    
    if(compute.hessian) {
        se <- rep(NA, 1 + NCOL(design.mat))
        se[c(1, which(valid.f) + 1)] <- sqrt(diag(solve(fit$hessian)))
    } else {
        se <- NULL
    }
    
    return(list(fitval = fv, est = est, coef = coef, se = se, ll = -fit$value))
    
}

#' Fit a gamma-poisson mixture model to the DNA-RNA counts. A gamma distribution
#' is fitted to the DNA counts and a Poisson noise is added for the RNA counts
#' which results in a Negative Binomial distribution for the RNA. The model is
#' fit simultaneously.
#'
#' @param dcounts the DNA counts
#' @param rcounts the RNA counts
#' @param ddepth DNA library size correction factors
#' @param rdepth RNA library size correction factors
#' @param ddesign.mat the design matrix for the DNA
#' @param rdesign.mat the design matrix for the RNA
#' @param compute.hessian if TRUE (default), compute the Hessian matrix of the
#' coefficients to facilitate coefficient-based hypothesis testing
#'
#' @return a list with components:
#' \itemize {
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
fit.dnarna.noctrlobs <- function(dcounts, rcounts,
                              ddepth, rdepth,
                              ddesign.mat, rdesign.mat,
                              compute.hessian) {
    
    ## set cost function
    if(model == "gamma.pois"){
        costfnDNA <- ll.dna.gamma.pois
        costfnRNA <- ll.rna.gamma.pois
    } else if(model == "ln.nb"){
        costfnDNA <- ll.dna.ln.nb
        costfnRNA <- ll.rna.ln.nb
    } else {
        stop("model ", model, " not available in fit.dnarna.wctrl.iter()")
    }
    
    ## filter invalid counts (NAs) from data and design
    valid.c <- (dcounts > 0) & !is.na(dcounts) & !is.na(rcounts)
    dcounts.valid <- dcounts[valid.c]
    rcounts.valid <- rcounts[valid.c]
    log.ddepth.valid <- log(ddepth[valid.c])
    log.rdepth.valid <- log(rdepth[valid.c])
    
    ## clean design matrix from unused factors: note that these should be
    valid.df <- apply(ddesign.mat[valid.c,,drop=FALSE], 2, function(x) !all(x==0))
    valid.rf <- apply(rdesign.mat[valid.c,,drop=FALSE], 2, function(x) !all(x==0))
    
    ddmat.valid <- ddesign.mat[valid.c,valid.df,drop=FALSE]
    rdmat.valid <- rdesign.mat[valid.c,valid.rf,drop=FALSE]
    
    ## Initialize parameter vector with a guess
    guess <- rep(0, 1 + NCOL(ddmat.valid) + NCOL(rdmat.valid))
    
    fit <- optim(par = guess,
                 fn = ll.dnarna.noctrl,
                 costfnDNA = costfnDNA,
                 costfnRNA = costfnRNA,
                 dcounts = dcounts.valid,
                 rcounts = rcounts.valid,
                 log.ddepth = log.ddepth.valid,
                 log.rdepth = log.rdepth.valid,
                 rctrlscale = rctrlscale.valid,
                 ddesign.mat = ddmat.valid,
                 rdesign.mat = rdmat.valid,
                 rctrldesign.mat = rdmat.ctrl.valid,
                 method = "BFGS",
                 hessian = compute.hessian)
    
    ## split parameters to the two parts of the model
    d.par <- fit$par[2:(NCOL(ddmat.valid) + 1)]
    r.par <- fit$par[-(1:(NCOL(ddmat.valid) + 1))]
    
    d.est <- rep(NA, length(dcounts))
    d.est[valid.c] <- exp(fit$par[1] + (ddmat.valid %*% d.par))
    
    d.fitval <- d.est * ddepth
    
    r.est <- rep(NA, length(rcounts))
    r.est[valid.c] <- exp(rdmat.valid %*% r.par)
    r.est <- r.est * d.est
    
    r.fitval <- r.est * rdepth
    
    d.coef <- c(fit$par[1], rep(NA, NCOL(ddesign.mat)))
    d.coef[1 + which(valid.df)] <- d.par
    
    r.coef <- c(fit$par[1], rep(NA, NCOL(rdesign.mat)))
    r.coef[1 + which(valid.rf)] <- r.par
    
    ## standard error of the estimates
    if (compute.hessian) {
        se <- sqrt(diag(solve(fit$hessian)))
        
        d.se <- rep(NA, 1 + NCOL(ddesign.mat))
        d.se[c(1, 1 + which(valid.df))] <- se[1:(NCOL(ddmat.valid) + 1)]
        
        r.se <- rep(NA, 1 + NCOL(rdesign.mat))
        r.se[c(1, 1 + which(valid.rf))] <- se[-(2:(NCOL(rdmat.valid) + 1))]
    } else {
        d.se <- NULL
        r.se <- NULL
    }
    
    return(list(
        d.fitval = d.fitval, d.est = d.est, d.coef = d.coef, d.se = d.se,
        r.fitval = r.fitval, r.est = r.est, r.coef = r.coef, r.se = r.se,
        ll = -fit$value))
}

#' Fit a gamma-poisson mixture model to the DNA-RNA counts 
#' with control enhancers. 
#' 
#' A gamma distribution
#' is fitted to the DNA counts and a Poisson noise is added for the RNA counts
#' which results in a Negative Binomial distribution for the RNA. The model is
#' fit simultaneously.
#'
#' @param dcounts the DNA counts of case and control enhancers
#' (matrix enhancers x samples)
#' @param rcounts the RNA counts of case and control enhancers
#' (matrix enhancers x samples)
#' @param ddepth DNA library size correction factors
#' @param rdepth RNA library size correction factors
#' @param ddesign.mat the design matrix for the DNA
#' @param rdesign.mat the design matrix for the RNA
#' @param rdesign.ctrl.mat the design matrix with additional factors for
#' the control RNA enhancers, beyond rdesign.mat
#' @param compute.hessian if TRUE (default), compute the Hessian matrix of the
#' coefficients to facilitate coefficient-based hypothesis testing
#'
#' @return a list with components:
#' \itemize {
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
fit.dnarna.wctrlobs.iter <- function(model, dcounts, rcounts,
                                  ddepth, rdepth, rctrlscale=NULL,
                                  ddesign.mat, rdesign.mat, rdesign.ctrl.mat,
                                  theta.d.ctrl.prefit,
                                  compute.hessian) {
    ## set cost function
    if(model == "gamma.pois"){
        costfnDNA <- ll.dna.gamma.pois
        costfnRNA <- ll.rna.gamma.pois
    } else if(model == "ln.nb"){
        costfnDNA <- ll.dna.ln.nb
        costfnRNA <- ll.rna.ln.nb
    } else {
        stop("model ", model, " not available in fit.dnarna.wctrl.iter()")
    }
    
    ## filter invalid counts (NAs) from data and design
    # this operation is performed on the case enhancer, the first row
    valid.c <- (dcounts[1,] > 0) & !is.na(dcounts[1,]) & !is.na(rcounts[1,])
    dcounts.valid <- dcounts[,valid.c]
    rcounts.valid <- rcounts[,valid.c]
    log.ddepth.valid <- log(ddepth[valid.c])
    log.rdepth.valid <- log(rdepth[valid.c])
    
    ## clean design matrix from unused factors: note that these should be
    valid.df <- apply(ddesign.mat[valid.c,,drop=FALSE], 2, 
                      function(x) !all(x==0))
    valid.rf <- apply(rdesign.mat[valid.c,,drop=FALSE], 2, 
                      function(x) !all(x==0))
    valid.rf.ctrl <- apply(rdesign.ctrl.mat[valid.c,,drop=FALSE], 2, 
                           function(x) !all(x==0))
    
    ddmat.valid <- ddesign.mat[valid.c,valid.df,drop=FALSE]
    rdmat.valid <- rdesign.mat[valid.c,valid.rf,drop=FALSE]
    rdmat.ctrl.valid <- rdesign.ctrl.mat[valid.c,valid.rf,drop=FALSE]
    
    ## Iterative parameter estimation: coordinate ascent
    # Iterate DNA and RNA model estimation
    # Initialize DNA model parameter vector with a guess
    d.par <- rep(0, 1 + NCOL(ddmat.valid))
    r.par <- rep(0, 1 + NCOL(rdmat.valid))
    r.ctrl.par <- rep(0, 1 + NCOL(rdmat.ctrl.valid))
    
    llold <- -Inf
    llnew <- 0
    iter <- 1
    converged <- TRUE
    RELTOL <- 10^(-4)
    MAXITER <- 1000
    while(llnew > llold+llold*RELTOL & iter < MAXITER) {
        ## estimate dna model condition on rna model
        fit <- optim(par = d.par,
                     fn = ll.dna.wctrl,
                     costfnDNA = costfnDNA,
                     costfnRNA = costfnRNA,
                     theta.d.ctrl.prefit = theta.d.ctrl.prefit,
                     theta.r = r.par,
                     theta.r.ctrl = r.ctrl.par,
                     dcounts = dcounts.valid,
                     rcounts = rcounts.valid,
                     log.ddepth = log.ddepth.valid,
                     log.rdepth = log.rdepth.valid,
                     rctrlscale = rctrlscale.valid,
                     ddesign.mat = ddmat.valid,
                     rdesign.mat = rdmat.valid,
                     rdesign.ctrl.mat = rdmat.ctrl.valid,
                     method = "BFGS",
                     hessian = compute.hessian)
        
        d.par <- fit$par[seq(1, 1+length(d.par))]
        
        ## estimate rna model conditioned on dna model
        fit <- optim(par = c(r.par, r.ctrl.par),
                     fn = ll.rna.wctrl,
                     costfnDNA = costfnDNA,
                     costfnRNA = costfnRNA,
                     theta.d = d.par,
                     theta.d.ctrl.prefit = theta.d.ctrl.prefit,
                     rcounts = rcounts.valid,
                     log.rdepth = log.rdepth.valid,
                     rctrlscale = rctrlscale.valid,
                     ddesign.mat = ddmat.valid,
                     rdesign.mat = rdmat.valid,
                     rdesign.ctrl.mat = rdmat.ctrl.valid,
                     method = "BFGS",
                     hessian = compute.hessian)
        
        r.par <- fit$par[seq(1, 1+length(r.par))]
        r.ctrl.par <- fit$par[seq(1+length(r.par)+1, 
                                  1+length(r.par)+length(r.ctrl.par))]
        
        ## update iteration convergence reporters
        llold <- llnew
        # need to compute because ll.rna.gammaPoisson.wctrl does not evaluate 
        # the full likelihood
        llnew <- fncost3(
            theta = c(d.par, r.par, r.par.ctrl), 
            theta.d.ctrl.prefit = theta.d.ctrl.prefit,
            dcounts = dcounts.valid, rcounts = rcounts.valid,
            log.ddepth = log.ddepth.valid,
            log.rdepth = log.rdepth.valid,
            rctrlscale = rctrlscale.valid,
            ddesign.mat = ddmat.valid,
            rdesign.mat = rdmat.valid,
            rdesign.ctrl.mat = rdmat.ctrl.valid)
        iter <- iter + 1
        if(iter == MAXITER & llnew > llold+llold*RELTOL) {
            converged <- FALSE
        }
    }
    
    d.est <- rep(NA, length(dcounts))
    d.est[valid.c] <- exp(fit$par[1] + (ddmat.valid %*% d.par))
    
    d.fitval <- d.est * ddepth
    
    r.est <- rep(NA, length(rcounts))
    r.est[valid.c] <- exp(rdmat.valid %*% r.par)
    r.est <- r.est * d.est
    
    r.fitval <- r.est * rdepth
    
    d.coef <- c(fit$par[1], rep(NA, NCOL(ddesign.mat)))
    d.coef[1 + which(valid.df)] <- d.par
    
    r.coef <- c(fit$par[1], rep(NA, NCOL(rdesign.mat)))
    r.coef[1 + which(valid.rf)] <- r.par
    
    ## standard error of the estimates
    if (compute.hessian) {
        se <- sqrt(diag(solve(fit$hessian)))
        
        d.se <- rep(NA, 1 + NCOL(ddesign.mat))
        d.se[c(1, 1 + which(valid.df))] <- se[1:(NCOL(ddmat.valid) + 1)]
        
        r.se <- rep(NA, 1 + NCOL(rdesign.mat))
        r.se[c(1, 1 + which(valid.rf))] <- se[-(2:(NCOL(rdmat.valid) + 1))]
    } else {
        d.se <- NULL
        r.se <- NULL
    }
    
    return(list(
        d.fitval = d.fitval, d.est = d.est, d.coef = d.coef, d.se = d.se,
        r.fitval = r.fitval, r.est = r.est, r.coef = r.coef, r.se = r.se,
        converged = converged,
        ll = -fit$value))
}

#' Fit DNA models of control enhancers
#' 
#' For later usage in case enhancer models.
#'
#' @param dcounts the DNA counts of case and control enhancers
#' (matrix ctrl enhancers x samples)
#' @param rcounts the RNA counts of case and control enhancers
#' (matrix ctrl enhancers x samples)
#' @param ddepth DNA library size correction factors
#' @param rdepth RNA library size correction factors
#' @param ddesign.mat the design matrix for the DNA
#' @param rdesign.mat the design matrix for the RNA
#' @param rdesign.ctrl.mat the design matrix with additional factors for
#' the control RNA enhancers, beyond rdesign.mat
#' @param compute.hessian if TRUE (default), compute the Hessian matrix of the
#' coefficients to facilitate coefficient-based hypothesis testing
#'
#' @return a list with components:
#' \itemize {
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
fit.dnactrl.wctrl.iter <- function(model, dcounts, rcounts,
                                   ddepth, rdepth, 
                                   ddesign.mat, rdesign.mat) {
    
    ## set cost function
    if(model == "gamma.pois"){
        fncost1 <- ll.dna.gammaPoisson.wctrl
        fncost2 <- ll.rna.gammaPoisson.wctrl
        fncost3 <- ll.dnarna.gammaPoisson.wctrl
    } else {
        stop("model ", model, " not available in fit.dnactlr.wctrl.iter()")
    }
    
    ## filter invalid counts (NAs) from data and design
    # this operation is performed on the case enhancer, the first row
    valid.c <- (dcounts[1,] > 0) & !is.na(dcounts[1,]) & !is.na(rcounts[1,])
    dcounts.valid <- dcounts[,valid.c]
    rcounts.valid <- rcounts[,valid.c]
    log.ddepth.valid <- log(ddepth[valid.c])
    log.rdepth.valid <- log(rdepth[valid.c])
    
    ## clean design matrix from unused factors: note that these should be
    valid.df <- apply(ddesign.mat[valid.c,,drop=FALSE], 2, 
                      function(x) !all(x==0))
    valid.rf <- apply(rdesign.mat[valid.c,,drop=FALSE], 2, 
                      function(x) !all(x==0))
    
    ddmat.valid <- ddesign.mat[valid.c,valid.df,drop=FALSE]
    rdmat.valid <- rdesign.mat[valid.c,valid.rf,drop=FALSE]
    
    ## Iterative parameter estimation: coordinate ascent
    # Iterate DNA and RNA model estimation
    # Initialize DNA model parameter vector with a guess
    d.par <- rep(0, 1 + NCOL(ddmat.valid))
    r.par <- rep(0, 1 + NCOL(rdmat.valid))
    
    llold <- -Inf
    llnew <- 0
    iter <- 1
    converged <- TRUE
    RELTOL <- 10^(-4)
    MAXITER <- 1000
    while(llnew > llold+llold*RELTOL & iter < MAXITER) {
        ## estimate dna model for each control enhancer
        fits <- bplappyl(seq_len(NROW(dcounts)), function(i) {
            optim(par = d.par[i,],
                  fn = fncost1,
                  theta.r = r.par,
                  dcounts = dcounts.valid[i,],
                  rcounts = rcounts.valid[i,],
                  log.ddepth = log.ddepth.valid,
                  log.rdepth = log.rdepth.valid,
                  ddesign.mat = ddmat.valid,
                  rdesign.mat = rdmat.valid,
                  method = "BFGS",
                  hessian = compute.hessian)
        })
        
        d.par <- do.call(cbind, lapply(fits, function(x) {
            $par[seq(1, NCOL(d.par))]
        }))
        
        ## estimate rna model conditioned on dna model
        fit <- optim(par = r.par,
                     fn = fncost2,
                     theta.d = d.par,
                     theta.d.ctrl.prefit = theta.d.ctrl.prefit,
                     rcounts = rcounts.valid,
                     log.rdepth = log.rdepth.valid,
                     rctrlscale = rctrlscale.valid,
                     ddesign.mat = ddmat.valid,
                     rdesign.mat = rdmat.valid,
                     method = "BFGS",
                     hessian = compute.hessian)
        
        r.par <- fit$par[seq(1, 1+length(r.par))]
        
        ## update iteration convergence reporters
        llold <- llnew
        # need to compute because ll.rna.gammaPoisson.wctrl does not evaluate 
        # the full likelihood
        llnew <- fncost3(
            theta = c(d.par, r.par, r.par.ctrl), 
            dcounts = dcounts.valid, rcounts = rcounts.valid,
            log.ddepth = log.ddepth.valid,
            log.rdepth = log.rdepth.valid,
            ddesign.mat = ddmat.valid,
            rdesign.mat = rdmat.valid)
        iter <- iter + 1
        if(iter == MAXITER & llnew > llold+llold*RELTOL) {
            converged <- FALSE
        }
    }
    
    d.est <- rep(NA, length(dcounts))
    d.est[valid.c] <- exp(fit$par[1] + (ddmat.valid %*% d.par))
    
    d.fitval <- d.est * ddepth
    
    r.est <- rep(NA, length(rcounts))
    r.est[valid.c] <- exp(rdmat.valid %*% r.par)
    r.est <- r.est * d.est
    
    r.fitval <- r.est * rdepth
    
    d.coef <- c(fit$par[1], rep(NA, NCOL(ddesign.mat)))
    d.coef[1 + which(valid.df)] <- d.par
    
    r.coef <- c(fit$par[1], rep(NA, NCOL(rdesign.mat)))
    r.coef[1 + which(valid.rf)] <- r.par
    
    return(list(
        d.fitval = d.fitval, d.est = d.est, d.coef = d.coef, d.se = NULL,
        r.fitval = r.fitval, r.est = r.est, r.coef = r.coef, r.se = NULL,
        converged = converged,
        ll = -fit$value))
}

#' A wrapper function for fitting a point-estimation based model, where a
#' point estimate for DNA levels is fitted from DNA counts, then used for RNA
#' estimation, without joint fitting.
#'
#' @param dnaFn the fitting function for the DNA counts. Should have the same API
#' as fit.gammaDNA (default) and fit.lnDNA
#' @param rnaFn the fitting function for the RNA counts, should have the same API
#' as fit.nbRNA (default)
#' @param dcounts the DNA count data
#' @param rcounts the RNA count data
#' @param ddepth the DNA library size correction factors
#' @param rdepth the RNA library size correction factors
#' @param ddesign.mat the DNA design matrix
#' @param rdesign.mat the RNA design matrix
#' @param compute.hessian if TRUE (default), compute the Hessian matrix of the
#' coefficients to facilitate coefficient-based hypothesis testing
#'
#' @return a list:
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
fit.separate <- function(dnaFn=fit.gammaDNA, rnaFn=fit.nbRNA,
                         dcounts, rcounts,
                         ddepth, rdepth,
                         ddesign.mat, rdesign.mat,
                         compute.hessian) {
    
    dna.fit <- dnaFn(dcounts = dcounts,
                     depth = ddepth,
                     design.mat = ddesign.mat,
                     compute.hessian = compute.hessian)
    
    rna.fit <- rnaFn(rcounts = rcounts,
                     depth = rdepth,
                     d.est = dna.fit$est,
                     design.mat = rdesign.mat,
                     compute.hessian = compute.hessian)
    
    return(list(
        d.fitval = dna.fit$fitval, d.est = dna.fit$est,
        d.coef = dna.fit$coef, d.se = dna.fit$se,
        r.fitval = rna.fit$fitval, r.est = rna.fit$est,
        r.coef = rna.fit$coef, r.se = rna.fit$se,
        ll = -(rna.fit$ll + dna.fit$ll)))
}
