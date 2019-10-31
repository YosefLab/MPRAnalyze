ll.rna.scale.nb <- function(theta, dcounts, rcounts, 
                            log.ddepth, log.rdepth, rdesign.mat) {
    log.d.est <- log(dcounts) - log.ddepth
    log.r.est <- log.d.est + rep(((rdesign.mat %*% theta[-1]) + log.rdepth), 
                                    NCOL(dcounts))
    
    ll <- sum(dnbinom(x = rcounts,
                      size = exp(exp(theta[1])),
                      mu = exp(log.r.est),
                      log = TRUE))
    return(-ll)
}

cost.rna.scale <- function(theta, dcounts, rcounts,
                        llfnRNA,
                        log.ddepth, log.rdepth, rctrlscale=NULL,
                        rdesign.mat,
                        rdesign.ctrl.mat=NULL) {

    ## compute likelihood
    # likelihood of case rna observations
    r.ll <- ll.rna.scale.nb(theta = c(theta, rctrlscale),
                    dcounts = dcounts,
                    rcounts = rcounts,
                    log.ddepth = log.ddepth,
                    log.rdepth = log.rdepth,
                    rdesign.mat = cbind(rdesign.mat, rdesign.ctrl.mat)) 
    
    return(r.ll)
}

fit.dnarna.noctrlobs.scale <- function(model, dcounts, rcounts,
                                 ddepth, rdepth, rctrlscale,
                                 rdesign.mat, rdesign.ctrl.mat, 
                                 compute.hessian) {
    # ## get likelihood function
    # ll.funs <- get.ll.functions(model)
    
    ## filter invalid counts (NAs) from data and design
    valid.c <- (dcounts > 0) & !is.na(dcounts) & (rcounts > 0) & !is.na(rcounts)
    
    dcounts.valid <- dcounts[valid.c]
    rcounts.valid <- rcounts[1,valid.c,drop=FALSE]
    log.ddepth.valid <- log(ddepth[valid.c])
    log.rdepth.valid <- log(rdepth[valid.c])
    
    ## clean design matrix from unused factors
    valid.rf <- apply(rdesign.mat[valid.c,,drop=FALSE], 2,
                      function(x) sum(x!=0) > 1)
    if(!any(valid.rf)) {
        valid.rf <- apply(rdesign.mat[valid.c,,drop=FALSE], 2,
                          function(x) !all(x==0))
    }
    rdmat.valid <- rdesign.mat[valid.c,valid.rf,drop=FALSE]
    rdmat.ctrl.valid <- rdesign.ctrl.mat[valid.c,,drop=FALSE]
    
    ## Initialize parameter vector with a guess
    guess <- rep(0, 1 + NCOL(rdmat.valid))
    
    ## optimize
    suppressWarnings(
        fit <- optim(
            par = guess,
            fn = cost.rna.scale, 
            llfnRNA = ll.rna.scale.nb,
            dcounts = dcounts.valid, 
            rcounts = rcounts.valid,
            log.ddepth = log.ddepth.valid, 
            log.rdepth = log.rdepth.valid,
            rctrlscale = rctrlscale,
            rdesign.mat = rdmat.valid,
            rdesign.ctrl.mat = rdmat.ctrl.valid,
            hessian = compute.hessian,
            method = "L-BFGS-B", control = list(maxit=1000), 
            lower=-23, upper=23)
    )

    r.coef <- c(fit$par[1], rep(NA, NCOL(rdesign.mat)))
    r.coef[1 + which(valid.rf)] <- fit$par[-1]
    r.df <- length(fit$par)
    
    ## standard error of the estimates
    se <- NULL
    r.se <- NULL
    if (compute.hessian) {
        se.comp <- tryCatch({
            se <- sqrt(diag(solve(fit$hessian)))
        }, error = function(e) return(e))
        
        if(!inherits(se.comp, "error")) {
            r.se <- c(se[1], rep(NA, NCOL(rdesign.mat)))
            r.se[1 + which(valid.rf)] <- se[seq(2, 1+NCOL(rdmat.valid))]
        } else {
            message(se.comp$message)
        }
    }
    
    return(list(
        r.coef = r.coef, r.se = r.se, r.df = r.df,
        r.ctrl.coef = NULL, r.ctrl.se = NULL, r.ctrl.df = 0,
        converged = fit$convergence,
        ll = -fit$value))
}

analyzeComparative.scale <- function(obj, rnaDesign, fit.se=FALSE, 
                               reducedDesign=NULL, correctControls=TRUE, 
                               verbose=TRUE) {
    if(!fit.se & is.null(reducedDesign)) {
        stop("Comparative analysis requires either a reduced design or fitting \
             the SE")
    }
    if(!is.null(reducedDesign)) {
        if (!isNestedDesign(full=rnaDesign, reduced=reducedDesign)) {
            stop("reduced design must be nested within the full RNA design")
        }
    }
    if(length(dnaDepth(obj)) != NCOL(dnaCounts(obj))) {
        obj <- estimateDepthFactors(obj, which.lib = "dna")
    }
    if(length(rnaDepth(obj)) != NCOL(rnaCounts(obj))) {
        obj <- estimateDepthFactors(obj, which.lib = "rna")
    }
    if(length(model(obj)) == 0) {
        obj <- autoChooseModel(obj)
    }
    if(!all(dim(dnaCounts(obj)) == dim(rnaCounts(obj)))) {
        stop("For scalable analysis, counts matrices must be of equal 
             dimensions and observations must be matched.")
    }
    
    obj@designs@rnaFull <- getDesignMat(design=rnaDesign, 
                                        annotations=rnaAnnot(obj))
    
    ## if controls are to be used: fit the control model
    # if(correctControls & any(controls(obj))) {
    #     message("Fitting controls-based background model...")
    #     obj@modelPreFits.dna.ctrl <- reformatModels(fit.dnarna.onlyctrl.iter(
    #         model=model(obj),
    #         dcounts = dnaCounts(obj)[controls(obj),],
    #         rcounts = rnaCounts(obj)[controls(obj),],
    #         ddepth=dnaDepth(obj),
    #         rdepth=rnaDepth(obj),
    #         ddesign.mat=obj@designs@dna,
    #         rdesign.mat=obj@designs@rnaFull, 
    #         d2rdesign.mat=obj@designs@dna2rna,
    #         BPPARAM = obj@BPPARAM,
    #         print.progress = verbose))
    #     
    #     obj@designs@rnaCtrlFull <- obj@designs@rnaFull
    #     obj@designs@rnaCtrlRed <- obj@designs@rnaFull
    #     obj@rnaCtrlScale <- obj@modelPreFits.dna.ctrl$r.coef[1,]
    #     theta.d.ctrl.prefit <- t(obj@modelPreFits.dna.ctrl$d.coef)
    # } else {
    theta.d.ctrl.prefit <- NULL
    # }
    
    ## Fit the full model (with SE extraction if fit.SE is on)
    message("Fitting model...")
    pb <- progress_bar$new(format = "[:bar] :percent (:current/:total)", 
                           total = NROW(dnaCounts(obj)), clear = !verbose)
    models <- bplapply(rownames(dnaCounts(obj)), function(rn) {
        pb$tick()
        tryCatch({
            return(fit.dnarna.noctrlobs.scale(
                model=model(obj),
                dcounts=dnaCounts(obj)[rn,,drop=FALSE],
                rcounts=rnaCounts(obj)[rn,,drop=FALSE],
                ddepth=dnaDepth(obj),
                rdepth=rnaDepth(obj),
                rctrlscale=obj@rnaCtrlScale,
                rdesign.mat=obj@designs@rnaFull,
                rdesign.ctrl.mat=obj@designs@rnaCtrlFull,
                compute.hessian=fit.se))
        }, error = function(err) {message("error fitting ", rn, ": ", err)})
    }, BPPARAM = obj@BPPARAM)
    names(models) <- rownames(dnaCounts(obj))
    obj@modelFits <-reformatModels.scale(models)
    
    if(!is.null(reducedDesign)) {
        obj@designs@rnaRed <- getDesignMat(design=reducedDesign, 
                                           annotations=rnaAnnot(obj))
        
        message("Fitting reduced model...")
        pb <- progress_bar$new(format = "[:bar] :percent (:current/:total)", 
                               total = NROW(dnaCounts(obj)), clear = verbose)
        models <- bplapply(rownames(dnaCounts(obj)), function(rn) {
            pb$tick()
            tryCatch({
                return(fit.dnarna.noctrlobs.scale(
                    model=model(obj),
                    dcounts=dnaCounts(obj)[rn,,drop=FALSE],
                    rcounts=rnaCounts(obj)[rn,,drop=FALSE], 
                    ddepth=dnaDepth(obj),
                    rdepth=rnaDepth(obj),
                    rctrlscale=obj@rnaCtrlScale,
                    rdesign.mat=obj@designs@rnaRed,
                    rdesign.ctrl.mat=obj@designs@rnaCtrlRed,
                    compute.hessian=FALSE))
            }, error = function(err) {message("error fitting: ", 
                                              rn, ": ", err)})
        }, BPPARAM = obj@BPPARAM)
        names(models) <- rownames(dnaCounts(obj))
        obj@modelFits.red <- reformatModels.scale(models)
    }
    
    message("Analysis Done!")
    return(obj)
}

reformatModels.scale <- function(models) {
    valid <- !vapply(models, is.null, TRUE)
    res <- list(
        ll = extractProp(models, "ll", valid),
        converged = extractProp(models, "converged", valid),
        
        r.coef = extractProp(models, "r.coef", valid),
        r.df = extractProp(models, "r.df", valid),
        r.se = extractProp(models, "r.se", valid),
        
        r.ctrl.coef = extractProp(models, "r.ctrl.coef", valid),
        r.ctrl.df = extractProp(models, "r.ctrl.df", valid),
        r.ctrl.se = extractProp(models, "r.ctrl.se", valid)
    )
    
    return(res)
}

testLrt.scale <- function(obj) {
    if(length(obj@modelFits) == 0 | length(obj@modelFits.red) == 0) {
        stop("An LRT analysis must be performed before computing the test")
    }
    
    message("Performing Likelihood Ratio Test...")
    
    ll.full <- obj@modelFits$ll
    ll.red <- obj@modelFits.red$ll
    df.rna.full <- obj@modelFits$r.df + obj@modelFits$r.ctrl.df + 
        length(obj@rnaCtrlScale)
    df.rna.red <- obj@modelFits.red$r.df + obj@modelFits.red$r.ctrl.df + 
        length(obj@rnaCtrlScale)
    df.full <- df.rna.full 
    df.red <- df.rna.red 
    
    lrt <- 2*(ll.full - ll.red)
    df <- df.full-df.red
    pval <- pchisq(lrt, df=df, lower.tail=FALSE)
    fdr <- p.adjust(pval, 'BH')
    
    res <- data.frame(statistic=lrt, pval=pval, fdr=fdr, df.test=df,
                      df.rna.full=df.rna.full, 
                      df.rna.red=df.rna.red)
    ## if condition is single term, extract the coefficient as logFC
    condition.name <- (colnames(obj@designs@rnaFull)
                       [!(colnames(obj@designs@rnaFull) %in% 
                              colnames(obj@designs@rnaRed))])
    if(length(condition.name) == 1) {
        ## single coefficient is the log Fold Change
        res$logFC <- getModelParameters_RNA(obj)[,condition.name]
    }
    
    return(res)
}