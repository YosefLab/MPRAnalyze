#' Fit dna and rna model to a given enhancer
#' 
#' Optim wrapper that performs numerical optimisation and extracts results.
#' Depending on the function chosen, a single optimisation or iterative estimation
#' are performed.
#' 
#' @name fit.dnarna
#' @rdname fit.dnarna
#' 
#' @aliases 
#' fit.dnarna.noctrlobs
#' fit.dnarna.wctrlobs.iter
#' fit.dnarna.onlyctrl.iter
#' 
#' @param model noise model
#' @param dcounts the DNA counts
#' @param rcounts the RNA counts
#' @param ddepth dna library size correction vector (numeric, samples)
#' @param rdepth rna library size correction vector (numeric, samples)
#' @param ddesign.mat the dna model design matrix (logical, samples x dna parameters)
#' @param rdesign.mat the rna model design matrix (logical, samples x rna parameters)
#' @param rdesign.ctrl.mat the control rna model design matrix 
#' (logical, samples x rna parameters)
#' @param theta.d.ctrl.prefit ctrl dna model parameters to condition likelihood on
#' (numeric, control enhancers x dna parameters)
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
NULL

#' @rdname fit.dnarna
fit.dnarna.noctrlobs <- function(model,
                                 dcounts, rcounts,
                                 ddepth, rdepth, rctrlscale,
                                 ddesign.mat, rdesign.mat, rdesign.ctrl.mat,
                                 theta.d.ctrl.prefit,
                                 compute.hessian) {
    
    ## set cost function
    if(model == "gamma.pois"){
        llfnDNA <- ll.dna.gamma.pois
        llfnRNA <- ll.rna.gamma.pois
    } else if(model == "ln.nb"){
        llfnDNA <- ll.dna.ln.nb
        llfnRNA <- ll.rna.ln.nb
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
    guess <- matrix(0, nrow=1, ncol=1 + NCOL(ddmat.valid) + NCOL(rdmat.valid))
    
    fit <- optim(par = guess,
                 fn = cost.dnarna, llfnDNA = llfnDNA, llfnRNA = llfnRNA,
                 dcounts = dcounts.valid, rcounts = rcounts.valid,
                 log.ddepth = log.ddepth.valid, log.rdepth = log.rdepth.valid,
                 rctrlscale = rctrlscale.valid,
                 ddesign.mat = ddmat.valid, rdesign.mat = rdmat.valid,
                 rctrldesign.mat = rdmat.ctrl.valid,
                 method = "BFGS", hessian = compute.hessian)
    
    ## split parameters to the two parts of the model
    d.par <- fit$par[seq(1, 1+NCOL(ddmat.valid))]
    r.par <- fit$par[seq(1+NCOL(ddmat.valid)+1,
                         1+NCOL(ddmat.valid)+NCOL(rdmat.valid)+1)]
    
    #d.est <- rep(NA, length(dcounts))
    #d.est[valid.c] <- exp(d.par[1] + (ddmat.valid %*% d.par[-1]))
    
    #d.fitval <- d.est * ddepth
    
    #r.est <- rep(NA, length(rcounts))
    #r.est[valid.c] <- exp(rdmat.valid %*% r.par)
    #r.est <- r.est * d.est
    
    #r.fitval <- r.est * rdepth
    
    d.coef <- c(fit$par[1], rep(NA, NCOL(ddesign.mat)))
    d.coef[1 + which(valid.df)] <- d.par
    d.df <- length(d.par)
    
    r.coef <- rep(NA, NCOL(rdesign.mat))
    r.coef[which(valid.rf)] <- r.par
    r.df <- length(r.par)
    
    ## standard error of the estimates
    if (compute.hessian) {
        se <- sqrt(diag(solve(fit$hessian)))
        
        d.se <- rep(NA, 1 + NCOL(ddesign.mat))
        d.se[c(1, 1 + which(valid.df))] <- se[seq(1, 1+NCOL(ddmat.valid))]
        
        r.se <- rep(NA, NCOL(rdesign.mat))
        r.se[which(valid.rf)] <- se[seq(1+NCOL(ddmat.valid)+1,
                                        1+NCOL(ddmat.valid)+NCOL(rdmat.valid)+1)]
    } else {
        d.se <- NULL
        r.se <- NULL
    }
    
    return(list(
        #d.fitval = d.fitval, d.est = d.est, d.coef = d.coef, d.se = d.se,
        #r.fitval = r.fitval, r.est = r.est, r.coef = r.coef, r.se = r.se,
        d.coef = d.coef, d.se = d.se, d.df = d.df,
        r.coef = r.coef, r.se = r.se, r.df = r.df,
        r.ctrl.coef = NULL, r.ctrl.se = NULL, r.ctrl.df = NULL,
        converged = fit$convergence,
        ll = -fit$value))
}

#' @rdname fit.dnarna
fit.dnarna.wctrlobs.iter <- function(model,
                                     dcounts, rcounts,
                                     ddepth, rdepth, rctrlscale,
                                     ddesign.mat, rdesign.mat, rdesign.ctrl.mat,
                                     theta.d.ctrl.prefit,
                                     compute.hessian) {
    ## set cost function
    if(model == "gamma.pois"){
        llfnDNA <- ll.dna.gamma.pois
        llfnRNA <- ll.rna.gamma.pois
    } else if(model == "ln.nb"){
        llfnDNA <- ll.dna.ln.nb
        llfnRNA <- ll.rna.ln.nb
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
    rdmat.ctrl.valid <- rdesign.ctrl.mat[valid.c,valid.rf.ctrl,drop=FALSE]
    
    ## Iterative parameter estimation: coordinate ascent
    # Iterate DNA and RNA model estimation
    # Initialize DNA model parameter vector with a guess
    d.par <- matrix(0, nrow=1, ncol=1 + NCOL(ddmat.valid))
    r.par <- matrix(0, nrow=1, ncol=NCOL(rdmat.valid))
    r.ctrl.par <- matrix(0, nrow=1, ncol=NCOL(rdmat.ctrl.valid))
    
    llold <- -Inf
    llnew <- 0
    iter <- 1
    converged <- TRUE
    RELTOL <- 10^(-8)
    MAXITER <- 1000
    while(llnew > llold+llold*RELTOL & iter < MAXITER) {
        ## estimate dna model condition on rna model
        dfit <- optim(par = d.par, fn = cost.dna.wctrl,
                     llfnDNA = llfnDNA, llfnRNA = llfnRNA,
                     theta.d.ctrl.prefit = theta.d.ctrl.prefit, theta.r = r.par,
                     theta.r.ctrl = r.ctrl.par,
                     dcounts = dcounts.valid, rcounts = rcounts.valid,
                     log.ddepth = log.ddepth.valid, log.rdepth = log.rdepth.valid,
                     rctrlscale = rctrlscale.valid,
                     ddesign.mat = ddmat.valid, rdesign.mat = rdmat.valid,
                     rdesign.ctrl.mat = rdmat.ctrl.valid,
                     method = "BFGS", hessian = compute.hessian)
        
        d.par <- dfit$par[seq(1, length(d.par))]
        
        ## estimate rna model conditioned on dna model
        rfit <- optim(par = c(r.par, r.ctrl.par), 
                      fn = cost.rna.wctrl, llfnRNA = llfnRNA,
                      theta.d = d.par, theta.d.ctrl.prefit = theta.d.ctrl.prefit,
                      rcounts = rcounts.valid,
                      log.rdepth = log.rdepth.valid, rctrlscale = rctrlscale.valid,
                      ddesign.mat = ddmat.valid, rdesign.mat = rdmat.valid,
                      rdesign.ctrl.mat = rdmat.ctrl.valid,
                      method = "BFGS", hessian = compute.hessian)
        
        r.par <- rfit$par[seq(1, 1+length(r.par))]
        r.ctrl.par <- rfit$par[seq(1+length(r.par)+1, 
                                   1+length(r.par)+length(r.ctrl.par))]
        
        ## update iteration convergence reporters
        llold <- llnew
        llnew <- -cost.dnarna.wctrl(
            theta = c(d.par, r.par, r.ctrl.par),
            theta.d.ctrl.prefit = theta.d.ctrl.prefit,
            llfnDNA = llfnDNA, llfnRNA = llfnRNA,
            dcounts = dcounts, rcounts = rcounts,
            log.ddepth = log.ddepth, log.rdepth = log.rdepth, 
            ddesign.mat = ddesign.mat, rdesign.mat = rdesign.mat, 
            rdesign.ctrl.mat = rdesign.ctrl.mat)
        iter <- iter + 1
        if(iter == MAXITER & llnew > llold+llold*RELTOL) {
            converged <- FALSE
        }
    }
    
    #d.est <- rep(NA, NCOL(dcounts))
    #d.est[valid.c] <- exp(fit$par[1] + (ddmat.valid %*% d.par))
    
    #d.fitval <- d.est * ddepth
    
    #r.est <- rep(NA, length(rcounts))
    #r.est[valid.c] <- exp(rdmat.valid %*% r.par)
    #r.est <- r.est * d.est
    
    #r.fitval <- r.est * rdepth
    
    d.coef <- c(d.par[1], rep(NA, NCOL(ddesign.mat)))
    d.coef[1 + which(valid.df)] <- d.par[seq(2, length(d.par))]
    d.df <- length(d.par)
    
    r.coef <- rep(NA, NCOL(rdesign.mat))
    r.coef[which(valid.rf)] <- r.par
    r.df <- length(r.par)
    
    r.ctrl.coef <- rep(NA, NCOL(rdesign.ctrl.mat))
    r.ctrl.coef[which(valid.rf.ctrl)] <- r.ctrl.par
    r.ctrl.df <- length(r.ctrl.par)
    
    ## standard error of the estimates
    if (compute.hessian) {
        se <- sqrt(diag(solve(fit$hessian)))
        
        d.se <- rep(NA, 1 + NCOL(ddesign.mat))
        d.se[c(1, 1 + which(valid.df))] <- se[seq(1, 1+NCOL(ddmat.valid))]
        
        r.se <- rep(NA, NCOL(rdesign.mat))
        r.se[which(valid.rf)] <- 
            se[seq(1+NCOL(ddmat.valid)+1,
                   1+NCOL(ddmat.valid)+NCOL(rdmat.valid))]
        
        r.ctrl.se <- rep(NA, NCOL(rdesign.ctrl.mat))
        r.ctrl.se[which(valid.rf.ctrl)] <- 
            se[seq(1+NCOL(ddmat.valid)+NCOL(rdmat.valid)+1,
                   1+NCOL(ddmat.valid)+NCOL(rdmat.valid)+NCOL(rdmat.ctrl.valid))]
    } else {
        d.se <- NULL
        r.se <- NULL
    }
    
    return(list(
        #d.fitval = d.fitval, d.est = d.est, d.coef = d.coef, d.se = d.se,
        #r.fitval = r.fitval, r.est = r.est, r.coef = r.coef, r.se = r.se,
        d.coef = d.coef, d.se = d.se, d.df = d.df,
        r.coef = r.coef, r.se = r.se, r.df = r.df,
        r.ctrl.coef = r.ctrl.coef, r.ctrl.se = r.ctrl.se, r.ctrl.df = r.ctrl.df,
        converged = converged,
        ll = llnew))
}

#' @rdname fit.dnarna
fit.dnarna.onlyctrl.iter <- function(model, dcounts, rcounts,
                                     ddepth, rdepth, rctrlscale=NULL,
                                     ddesign.mat, rdesign.mat, rdesign.ctrl.mat=NULL,
                                     theta.d.ctrl.prefit=NULL, compute.hessian=NULL) {
    
    ## set cost function
    if(model == "gamma.pois"){
        llfnDNA <- ll.dna.gamma.pois
        llfnRNA <- ll.rna.gamma.pois
    } else if(model == "ln.nb"){
        llfnDNA <- ll.dna.ln.nb
        llfnRNA <- ll.rna.ln.nb
    } else {
        stop("model ", model, " not available in fit.dnarna.wctrl.iter()")
    }
    
    ## filter invalid counts (NAs) from data and design
    valid.c.d <- (dcounts > 0) & !is.na(dcounts) & !is.na(rcounts) 
    valid.c.r <- apply((dcounts > 0) & !is.na(dcounts) & !is.na(rcounts), 2, any) 
    #dcounts.valid <- dcounts[,valid.c]
    #rcounts.valid <- rcounts[,valid.c]
    #log.ddepth.valid <- log(ddepth[valid.c])
    #log.rdepth.valid <- log(rdepth[valid.c])
    log.ddepth <- log(ddepth)
    log.rdepth <- log(rdepth)
    
    ## clean design matrix from unused factors: note that these should be
    valid.rf <- apply(rdesign.mat[valid.c.r,,drop=FALSE], 2, 
                      function(x) !all(x==0))
    
    #ddmat.valid <- ddesign.mat[valid.c,valid.df,drop=FALSE]
    #rdmat.valid <- rdesign.mat[valid.c,valid.rf,drop=FALSE]
    
    ## Iterative parameter estimation: coordinate ascent
    # Iterate DNA and RNA model estimation
    # Initialize DNA model parameter vector with a guess
    d.par <- matrix(0, nrow=NROW(dcounts), ncol=1+NCOL(ddesign.mat))
    r.par <- rep(0, NCOL(rdesign.mat))
    
    llold <- -Inf
    llnew <- 0
    iter <- 1
    converged <- TRUE
    RELTOL <- 10^(-4)
    MAXITER <- 1000
    while(llnew > llold+llold*RELTOL & iter < MAXITER) {
        ## estimate dna model for each control enhancer
        # TODO: parallelize
        dfits <- lapply(seq_len(NROW(dcounts)), function(i) {
            valid.df <- apply(ddesign.mat[valid.c.d[i,],,drop=FALSE], 2, 
                              function(x) !all(x==0))
            fit <- optim(par = d.par[i,], 
                         fn = cost.dna, theta.r = r.par,
                         llfnDNA = llfnDNA, llfnRNA = llfnRNA,
                         dcounts = dcounts[i,valid.c.d[i,]], 
                         rcounts = rcounts[i,valid.c.d[i,]], 
                         log.ddepth = log.ddepth[valid.c.d[i,]], 
                         log.rdepth = log.rdepth[valid.c.d[i,]],
                         ddesign.mat = ddesign.mat[valid.c.d[i,],valid.df,drop=FALSE], 
                         rdesign.mat = rdesign.mat[valid.c.d[i,],valid.rf,drop=FALSE], 
                         method = "BFGS", hessian = FALSE)
            return(fit)
        })
        
        
        d.par <- matrix(0, nrow=length(dfits),
                        ncol=NCOL(ddesign.mat)+1)
        d.par[,1] <- sapply(dfits, function(x) x$par[1] )
        for(i in seq_len(NROW(dcounts))){
            valid.df <- apply(ddesign.mat[valid.c.d[i,],,drop=FALSE], 2, 
                              function(x) !all(x==0))
            d.par[i, 1 + which(valid.df)] <- dfits[[i]]$par[-1]
        }
        
        ## estimate rna model conditioned on dna model
        rfit <- optim(par = r.par, fn = cost.rna, theta.d = d.par,
                     llfnRNA = llfnRNA, 
                     rcounts = rcounts[,valid.c.r], 
                     log.rdepth = log.rdepth[valid.c.r],
                     ddesign.mat = ddesign.mat[valid.rf,,drop=FALSE], 
                     rdesign.mat = rdesign.mat[valid.rf,,drop=FALSE], 
                     method = "BFGS", hessian = FALSE)
        
        r.par <- rfit$par
        names(r.par) <- NULL
        
        ## update iteration convergence reporters
        llold <- llnew
        llnew <- -sum(sapply(seq_len(NROW(dcounts)), function(i) {cost.dnarna(
            theta = c(d.par[i,], theta.r = r.par),
            llfnDNA = llfnDNA, llfnRNA = llfnRNA,
            dcounts = dcounts[i,], 
            rcounts = rcounts[i,], 
            log.ddepth = log.ddepth.valid, log.rdepth = log.rdepth.valid, 
            rctrlscale = NULL,
            ddesign.mat = ddmat.valid, rdesign.mat = rdmat.valid,
            rctrldesign.mat = NULL)
        }))
        iter <- iter + 1
        if(iter == MAXITER & llnew > llold+llold*RELTOL) {
            converged <- FALSE
        }
    }
    
    #d.est <- matrix(NA, nrow = NROW(dcounts), ncol = NCOL(dcounts))
    #d.est[valid.c] <- exp(fit$par[1] + (ddmat.valid %*% d.par))
    
    #d.fitval <- d.est * ddepth
    
    #r.est <- matrix(NA, nrow = NROW(rcounts), ncol = NCOL(rcounts))
    #r.est[valid.c] <- exp(rdmat.valid %*% r.par)
    #r.est <- r.est * d.est
    
    #r.fitval <- r.est * rdepth
    
    d.coef <- d.par
    d.df <- sapply(dfits, function(x) length(x$fit) )
    
    r.coef <- rep(NA, NCOL(rdesign.mat))
    r.coef[which(valid.rf)] <- r.par
    r.df <- length(r.par)
    
    return(lapply(seq_len(NROW(d.coef)), function(i){
        list(
            #d.fitval = d.fitval, d.est = d.est, d.coef = d.coef, d.se = NULL,
            #r.fitval = r.fitval, r.est = r.est, r.coef = r.coef, r.se = NULL,
            d.coef = d.coef[i,], d.se = NULL, d.df = d.df[i],
            r.coef = r.coef, r.se = NULL, r.df = r.df,
            r.ctrl.coef = NULL, r.ctrl.se = NULL, r.ctrl.df = NULL,
            converged = converged,
            ll = NA)
    }))
}
