
#' fit models for a differential activity analysis
#' @param obj the MpraObject
#' @param model the model to use (TODO: if nul?)
#' @param dnaDesign the design to use for the DNA model (see details)
#' @param rnaDesign the design to use for the RNA model (see details)
#' @param condition the condition to do the differential analysis over. Only
#' used if design is provided as formulas.
#' @param dnaDesignNull used to explicitly specify the design of the null DNA model
#' @param rnaDesignNull used to explicitly specify the design of the null RNA model
#'
fit.differential <- function(obj, model=NULL, dnaDesign=NULL, rnaDesign=NULL,
                            condition=NULL,
                            dnaDesignNull=NULL, rnaDesignNull=NULL) {
    ## TODO: if condition is specified and design input is formulas, create design matrices
    ## for full (as is) and reduced (remove condition term). Also make sure the
    ## condition is the first term in the formula. This will be for cases where
    ## the comparison is all-vs-one

    ## TODO: if designs are all specified and condition is not, either create model matrices
    ##is the designs are formulas, or use the design matrices as is.
}

#' fit model for quantitative activity analysis
#' @param obj the MpraObject
#' @param model the model to fit
#' @param dnaDesign the design for the DNA model
#' @param rnaDesign the design for the RNA model
fit.quantitative <- function(obj, model=NULL, dnaDesign=NULL, rnaDesign=NULL) {
    ##TODO: fit a single model per enhancer
}

#' fit a log normal GLM to the DNA count data
#'
#' @param dcounts the DNA counts (integer, N)
#' @param depth library size correction factors (numeric, N)
#' @param design.mat the design matrix (logical, N x K)
#' @param control list of parameter to be passed to `optim` as
#' the control argument. NOTE: since this is a maximization problem,
#' control$fnscale should be set to a negative number (idealy -1)
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
fit.lnDNA <- function(dcounts, depth, design.mat) {
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
                 hessian = TRUE,
                 control = list(trace=TRUE))

    ## extract parameters
    est <- rep(NA, length(dcounts))
    est[valid.c] <- exp(dmat.valid %*% fit$par[-1])

    fv <- est * depth

    coef <- rep(NA, 1 + NCOL(design.mat))
    coef[c(1, which(valid.f) + 1)] <- fit$par
    names(coef) <- c("dispersion", colnames(design.mat))

    se <- rep(NA, 1 + NCOL(design.mat))
    se[c(1, 1 + which(valid.f))] <- sqrt(diag(solve(fit$hessian)))
    names(se) <- c("dispersion", colnames(design.mat))

    return(list(fitval = fv, est = est, coef = coef, se = se, ll = -fit$value))
}

#' fit a Gamma GLM to the DNA count data
#'
#' @param dcounts the DNA counts (integer, N)
#' @param depth library size correction factors (numeric, N)
#' @param design.mat the design matrix (logical, N x K)
#' @param control list of parameter to be passed to `optim` as
#' the control argument. NOTE: since this is a maximization problem,
#' control$fnscale should be set to a negative number (idealy -1)
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
fit.gammaDNA <- function(dcounts, depth, design.mat) {
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
                 hessian = TRUE,
                 control = list(trace=TRUE))

    ## extract parameters
    est <- rep(NA, length(dcounts))
    est[valid.c] <- exp((dmat.valid %*% fit$par[-1]) + fit$par[1])

    fv <- est * depth

    coef <- rep(NA, 1 + NCOL(design.mat))
    coef[c(1, which(valid.f) + 1)] <- fit$par
    ##TODO: name the parameters!

    se <- rep(NA, 1 + NCOL(design.mat))
    se[c(1, which(valid.f) + 1)] <- sqrt(diag(solve(fit$hessian)))
    ##TODO: name the parameters!

    return(list(fitval = fv, est = est, coef = coef, se = se, ll = -fit$value))
}

#' fit a negative binomial model to the RNA counts using a point estimator for
#' the DNA counts
#'
#' @param rcounts RNA counts
#' @param depth library size correction factors
#' @param d.est point estimates of the DNA base counts
#' @param design.mat the RNA design matrix
#'
#' @return a list with components:
#' \itemize {
#'     \item fitval: the fitten values
#'     \item est: the estimated RNA levels (controled for library size)
#'     \item coed: the fitted coefficients for the parameters
#'     \item se: the standard errors of the coefficients
#'     \item ll: the log likelihood of the model
#' }
fit.nbRNA <- function(rcounts, depth, d.est, design.mat) {
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
                 hessian = TRUE,
                 control = list(trace=TRUE))

    ## extract parameters
    est <- rep(NA, length(dcounts))
    est[valid.c] <- exp((dmat.valid %*% fit$par[-1]) + logd.est.valid)

    fv <- est * depth

    coef <- rep(NA, 1 + NCOL(design.mat))
    coef[c(1, which(valid.f) + 1)] <- fit$par

    se <- rep(NA, 1 + NCOL(design.mat))
    se[c(1, which(valid.f) + 1)] <- sqrt(diag(solve(fit$hessian)))

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
#' @param control control arguments to pass to `optim`. NOTE: since this is a
#' maximization problem, control$fnscale should be set to a negative number.
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
fit.mixture.gammaPoisson <- function(dcounts, rcounts,
                                     ddepth, rdepth,
                                     ddesign.mat, rdesign.mat) {
    ## filter invalid counts (NAs) from data and design
    valid.c <- (dcounts > 0) & !is.na(dcounts) & !is.na(rcounts)
    dcounts.valid <- dcounts[valid.c]
    rcounts.valid <- rcounts[valid.c]
    log.ddepth.valid <- log(ddepth[valid.c])
    log.rdepth.valid <- log(rdepth[valid.c])

    ## clean design matrix from unused factors: note that these should be
    valid.df <- apply(ddesign.mat[valid.c,], 2, function(x) !all(x==0))
    valid.rf <- apply(rdesign.mat[valid.c,], 2, function(x) !all(x==0))

    ddmat.valid <- ddesign.mat[valid.c,valid.df]
    rdmat.valid <- rdesign.mat[valid.c,valid.rf]

    ## Initialize parameter vector with a guess
    guess <- rep(0, 1 + NCOL(ddmat.valid) + NCOL(rdmat.valid))

    fit <- optim(par = guess,
                 fn = ll.mixture.gammaPoisson,
                 dcounts = dcounts.valid,
                 rcounts = rcounts.valid,
                 log.ddepth = log.ddepth.valid,
                 log.rdepth = log.rdepth.valid,
                 ddesign.mat = ddmat.valid,
                 rdesign.mat = rdmat.valid,
                 method = "BFGS",
                 hessian = TRUE,
                 control = list(trace=TRUE))

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
    se <- sqrt(diag(solve(fit$hessian)))

    d.se <- rep(NA, 1 + NCOL(ddesign.mat))
    d.se[c(1, 1 + which(valid.df))] <- se[1:(NCOL(ddmat.valid) + 1)]

    r.se <- rep(NA, 1 + NCOL(rdesign.mat))
    r.se[c(1, 1 + which(valid.rf))] <- se[-(2:(NCOL(rdmat.valid) + 1))]

    return(list(
        d.fitval = d.fitval, d.est = d.est, d.coef = d.coef, d.se = d.se,
        r.fitval = r.fitval, r.est = r.est, r.coef = r.coef, r.se = r.se,
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
                         ddesign.mat, rdesign.mat) {

    dna.fit <- dnaFn(dcounts = dcounts,
                     depth = ddepth,
                     design.mat = ddesign.mat)

    rna.fit <- rnaFn(rcounts = rcounts,
                     depth = rdepth,
                     d.est = dna.fit$est,
                     design.mat = rdesign.mat)

    return(list(
        d.fitval = dna.fit$fitval, d.est = dna.fit$est,
        d.coef = dna.fit$coef, d.se = dna.fit$se,
        r.fitval = rna.fit$fitval, r.est = rna.fit$est,
        r.coef = rna.fit$coef, r.se = rna.fit$se,
        ll = -(rna.fit$ll + dna.fit$ll)))
}
