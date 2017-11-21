#' differential activity analysis by condition
#' 
#' test condition effect
#' 
#' @name analyse.condition
#' @rdname analyse.condition
#' 
#' @aliases 
#' analyse.condition.lrt
#' analyse.condition.ttest
#'
#' @param obj the MpraObject
#' @param model the model to use (TODO: if nul?)
#' @param dnaDesign the design to use for the DNA model (see details)
#' @param rnaDesign the design to use for the RNA model (see details)
#' @param condition_totest the part of rnaDesign to test for (see details)
#' @param mode options are "scaled" (default") or "full", see details
#'
#' @details TODO - explain input designs (either formula or model matrices)
#' 
#' @details TODO rna control usage
#'
#' @details TODO - explain mode
#'
#' @return the MpraObject with fitted models and condition test
NULL

#' @rdname analyse.condition
#' @export
analyse.condition.lrt <- function(obj, model="gamma.pois", mode=NULL,
                                  dnaDesign=NULL, rnaDesign=NULL, condition_totest=NULL) {
    ## check depth is set
    if(length(obj@dnaDepth) == 0 | length(obj@rnaDepth) == 0) {
        stop("Library depth factors must be estimated or manually set")
    }
    ## check mode
    if(!is.null(obj@controls)) {
        # set default
        if(is.null(mode)) {
            mode <- "scaled"
        } else if(!mode %in% c("scaled", "full") ) {
            stop("Only mode 'scaled' and 'full' are supported if control enhancers are supplied")
        }
    } else {
        # set default
        if(is.null(mode)) {
            mode <- "full"
        }
        else if(!mode %in% c("full") ) {
            stop("Only mode 'full' is supported and sensible if no control enhancers are supplied")
        }
    }
    obj@mode <- mode
    obj@model <- model
    
    ## get design matrices
    obj@designs@dna <- getDesignMat(obj, dnaDesign)
    obj@designs@rnaFull <- getDesignMat(obj, rnaDesign, rna=TRUE)
    if(is.null(obj@controls)) {
        obj@designs@rnaRed <- getDesignMat(obj, rnaDesign, condition_totest, rna=TRUE)
        obj@designs@rnaCtrlFull <- NULL
        obj@designs@rnaCtrlRed <- NULL
        obj@modelPreFits.dna.ctrl <- NULL
        obj@controls.forfit <- NULL
        fitfun <- fit.dnarna.noctrlobs
    } else {
        message("Fit control enhancer background models")
        obj@modelPreFits.dna.ctrl <- fit.dnarna.onlyctrl.iter(
            model=model,
            dcounts = obj@dnaCounts[obj@controls,], 
            rcounts = obj@rnaCounts[obj@controls,],
            ddepth=obj@dnaDepth,
            rdepth=obj@rnaDepth,
            ddesign.mat=obj@designs@dna,
            rdesign.mat=obj@designs@rnaFull, 
            BPPARAM = obj@BPPARAM)
        if(obj@mode == "scaled") {
            # H1: rna - prefit ~ 1 + condition + batch
            # H0: rna - prefit ~ 1 + batch
            obj@designs@rnaRed <- getDesignMat(obj, rnaDesign, condition_totest, rna=TRUE)
            obj@designs@rnaCtrlFull <- obj@designs@rnaFull # used to correct prefit 
            obj@designs@rnaCtrlRed <- obj@designs@rnaFull # used to correct prefit 
            obj@rnaCtrlScale <- obj@modelPreFits.dna.ctrl[[1]]$r.coef
            obj@controls.forfit <- NULL
            fitfun <- fit.dnarna.noctrlobs
        } else if(obj@mode == "full") {
            # cs: case-control identity of enhancer
            # H1: rna ~ 1 + cs + condition + batch
            # H0: rna ~ 1 + condition + batch
            obj@designs@rnaRed <- obj@designs@rnaFull
            obj@designs@rnaCtrlFull <- getDesignMat(obj, rnaDesign, rna=TRUE)
            obj@designs@rnaCtrlRed <- NULL
            obj@rnaCtrlScale <- NULL
            obj@controls.forfit <- obj@controls
            fitfun <- fit.dnarna.wctrlobs.iter
        }
    }
    
    ## fit full models
    message("Fit full models")
    obj@modelFits <- bplapply(rownames(obj@dnaCounts), function(rn) {
        return(fitfun(model=model,
                      dcounts=obj@dnaCounts[rn,,drop=FALSE],
                      rcounts=rbind(obj@rnaCounts[rn,], obj@rnaCounts[obj@controls.forfit,]),
                      ddepth=obj@dnaDepth,
                      rdepth=obj@rnaDepth,
                      rctrlscale=obj@rnaCtrlScale,
                      ddesign.mat=obj@designs@dna,
                      rdesign.mat=obj@designs@rnaFull,
                      rdesign.ctrl.mat=obj@designs@rnaCtrlFull,
                      theta.d.ctrl.prefit=
                          do.call(rbind, lapply(obj@modelPreFits.dna.ctrl,
                                                function(x) x$d.coef)),
                      compute.hessian=FALSE))
    }, BPPARAM = obj@BPPARAM)
    names(obj@modelFits) <- rownames(obj@dnaCounts)
    
    ## fit reduced models
    message("Fit reduced models")
    obj@modelFits.red <- bplapply(rownames(obj@dnaCounts), function(rn) {
        return(fitfun(model=model,
                      dcounts=obj@dnaCounts[rn,,drop=FALSE],
                      rcounts=rbind(obj@rnaCounts[rn,], obj@rnaCounts[obj@controls.forfit,]),
                      ddepth=obj@dnaDepth,
                      rdepth=obj@rnaDepth,
                      rctrlscale=obj@rnaCtrlScale,
                      ddesign.mat=obj@designs@dna,
                      rdesign.mat=obj@designs@rnaRed,
                      rdesign.ctrl.mat=obj@designs@rnaCtrlRed,
                      theta.d.ctrl.prefit=
                          do.call(rbind, lapply(obj@modelPreFits.dna.ctrl,
                                                function(x) x$d.coef)),
                      compute.hessian=FALSE))
    }, BPPARAM = obj@BPPARAM)
    names(obj@modelFits.red) <- rownames(obj@dnaCounts)
    
    ## run lrt
    obj@results <- test.lrt(obj)
    
    return(obj)
}

#' @rdname analyse.condition
#' @export
analyse.condition.ttest <- function(obj, model="gamma.pois", mode="ttest",
                                    dnaDesign=NULL, rnaDesign=NULL, condition_totest=NULL) {
    ## check depth is set
    if(length(obj@dnaDepth) == 0 | length(obj@rnaDepth) == 0) {
        stop("Library depth factors must be estimated or manually set")
    }
    ## check mode
    if(!is.null(obj@controls)) {
        # set default
        if(is.null(mode)) {
            mode <- "scaled"
        } else if(!mode %in% c("scaled", "full") ) {
            stop("Only mode 'scaled' and 'full' are supported if control enhancers are supplied")
        }
    } else {
        # set default
        if(is.null(mode)) {
            mode <- "full"
        }
        else if(!mode %in% c("full") ) {
            stop("Only mode 'full' is supported and sensible if no control enhancers are supplied")
        }
    }
    obj@mode <- mode
    obj@model <- model
    
    ## get design matrices
    obj@designs@dna <- getDesignMat(obj, dnaDesign)
    obj@designs@rnaFull <- getDesignMat(obj, rnaDesign)
    if(is.null(obj@controls)) {
        obj@designs@rnaRed <- getDesignMat(obj, rnaDesign, condition_totest)
        obj@designs@rnaCtrlFull <- NULL
        obj@designs@rnaCtrlRed <- NULL
        obj@modelPreFits.dna.ctrl <- NULL
        obj@controls.forfit <- NULL
        fitfun <- fit.dnarna.noctrlobs
    } else {
        message("Fit control enhancer background models")
        obj@modelPreFits.dna.ctrl <- fit.dnarna.onlyctrl.iter(
            model=model,
            dcounts = obj@dnaCounts[obj@controls,], 
            rcounts = obj@rnaCounts[obj@controls,],
            ddepth=obj@dnaDepth,
            rdepth=obj@rnaDepth,
            ddesign.mat=obj@designs@dna,
            rdesign.mat=obj@designs@rnaFull, 
            BPPARAM = obj@BPPARAM)
        if(obj@mode == "scaled") {
            # H1: rna - prefit ~ 1 + condition + batch
            # H0: rna - prefit ~ 1 + batch
            obj@designs@rnaRed <- getDesignMat(obj, rnaDesign, condition_totest)
            obj@designs@rnaCtrlFull <- obj@designs@rnaFull # used to correct prefit 
            obj@designs@rnaCtrlRed <- obj@designs@rnaFull # used to correct prefit 
            obj@rnaCtrlScale <- obj@modelPreFits.dna.ctrl[[1]]$r.coef
            obj@controls.forfit <- NULL
            fitfun <- fit.dnarna.noctrlobs
        } else if(obj@mode == "full") {
            # cs: case-control identity of enhancer
            # H1: rna ~ 1 + cs + condition + batch
            # H0: rna ~ 1 + condition + batch
            obj@designs@rnaRed <- obj@designs@rnaFull
            obj@designs@rnaCtrlFull <- getDesignMat(obj, rnaDesign)
            obj@designs@rnaCtrlRed <- NULL
            obj@rnaCtrlScale <- NULL
            obj@controls.forfit <- obj@controls
            fitfun <- fit.dnarna.wctrlobs.iter
        }
    }
    
    ## fit full models
    message("Fit full models")
    obj@modelFits <- bplapply(rownames(obj@dnaCounts), function(rn) {
        return(fitfun(model=model,
                      dcounts=obj@dnaCounts[rn,,drop=FALSE],
                      rcounts=rbind(obj@rnaCounts[rn,], obj@rnaCounts[obj@controls.forfit,]),
                      ddepth=obj@dnaDepth,
                      rdepth=obj@rnaDepth,
                      rctrlscale=obj@rnaCtrlScale,
                      ddesign.mat=obj@designs@dna,
                      rdesign.mat=obj@designs@rnaFull,
                      rdesign.ctrl.mat=obj@designs@rnaCtrlFull,
                      theta.d.ctrl.prefit=
                          do.call(rbind, lapply(obj@modelPreFits.dna.ctrl,
                                                function(x) x$d.coef)),
                      compute.hessian=FALSE))
    }, BPPARAM = obj@BPPARAM)
    names(obj@modelFits) <- rownames(obj@dnaCounts)
    
    ## TODO run test on coefficients 
    
    return(obj)
}

#' quantitative/case-control analysis
#' 
#' test for difference between case and control enhancers or 
#' compare estimated magnitude against distribution of effect
#' on lower quantile.
#' 
#' @name analyse.casectrl
#' @rdname analyse.casectrl
#' 
#' @aliases 
#' analyse.casectrl.lrt
#' analyse.quant
#' 
#' @param obj the MpraObject
#' @param model the model to fit
#' @param dnaDesign the design for the DNA model
#' @param rnaDesign the design for the RNA model
#' 
#' @return the MpraObject with fitted models and condition test
NULL

#' @rdname analyse.casectrl
#' @export
analyse.casectrl.lrt <- function(obj, mode="scaled", model=NULL, dnaDesign=NULL, rnaDesign=~1) {
    
    ## check depth is set
    if(length(obj@dnaDepth) == 0 | length(obj@rnaDepth) == 0) {
        stop("Library depth factors must be estimated or manually set")
    }
    ## check mode
    if(!is.null(obj@controls)) {
        # set default
        if(is.null(mode)) {
            mode <- "scaled"
        } else if(!mode %in% c("scaled", "full") ) {
            stop("Only mode 'scaled' and 'full' are supported if control enhancers are supplied")
        }
    } else {
        # set default
        if(is.null(mode)) {
            mode <- "quant"
        }
        else if(!mode %in% c("quant") ) {
            stop("Only mode 'quant' is sensible if no control enhancers are supplied")
        }
    }
    obj@mode <- mode
    obj@model <- model
    
    ## get design matrices
    obj@designs@dna <- getDesignMat(obj, dnaDesign)
    obj@designs@rnaFull <- getDesignMat(obj, rnaDesign)
    if(is.null(obj@controls)) {
        stop("no control enhancers supplied")
    } else {
        message("Fit control enhancer background models")
        obj@modelPreFits.dna.ctrl <- fit.dnarna.onlyctrl.iter(
            model=model,
            dcounts = obj@dnaCounts[obj@controls,], 
            rcounts = obj@rnaCounts[obj@controls,],
            ddepth=obj@dnaDepth,
            rdepth=obj@rnaDepth,
            ddesign.mat=obj@designs@dna,
            rdesign.mat=obj@designs@rnaFull, 
            BPPARAM = obj@BPPARAM)
        if(obj@mode == "scaled") {
            # H1: rna - prefit ~ 1 + batch
            # H0: rna - prefit ~ 0
            obj@designs@rnaRed <- NULL
            obj@designs@rnaCtrlFull <- obj@designs@rnaFull # used to correct prefit 
            obj@designs@rnaCtrlRed <- obj@designs@rnaFull # used to correct prefit 
            obj@rnaCtrlScale <- obj@modelPreFits.dna.ctrl[[1]]$r.coef
            obj@controls.forfit <- NULL
            fitfun <- fit.dnarna.noctrlobs
        } else if(obj@mode == "full") {
            # cs: case-control identity of enhancer
            # H1: rna ~ 1 + cs + batch
            # H0: rna ~ 1 + batch
            obj@designs@rnaRed <- obj@designs@rnaFull
            obj@designs@rnaCtrlFull <- getDesignMat(obj, rnaDesign)
            obj@designs@rnaCtrlRed <- NULL
            obj@rnaCtrlScale <- NULL
            obj@controls.forfit <- obj@controls
            fitfun <- fit.dnarna.wctrlobs.iter
        }
    }
    
    ## fit full models
    message("Fit full models")
    obj@modelFits <- bplapply(rownames(obj@dnaCounts), function(rn) {
        return(fitfun(model=model,
                      dcounts=obj@dnaCounts[rn,,drop=FALSE],
                      rcounts=rbind(obj@rnaCounts[rn,], obj@rnaCounts[obj@controls.forfit,]),
                      ddepth=obj@dnaDepth,
                      rdepth=obj@rnaDepth,
                      rctrlscale=obj@rnaCtrlScale,
                      ddesign.mat=obj@designs@dna,
                      rdesign.mat=obj@designs@rnaFull,
                      rdesign.ctrl.mat=obj@designs@rnaCtrlFull,
                      theta.d.ctrl.prefit=
                          do.call(rbind, lapply(obj@modelPreFits.dna.ctrl,
                                                function(x) x$d.coef)),
                      compute.hessian=FALSE))
    }, BPPARAM = obj@BPPARAM)
    names(obj@modelFits) <- rownames(obj@dnaCounts)
    
    ## fit reduced models
    message("Fit reduced models")
    obj@modelFits.red <- bplapply(rownames(obj@dnaCounts), function(rn) {
        return(fitfun(model=model,
                      dcounts=obj@dnaCounts[rn,,drop=FALSE],
                      rcounts=rbind(obj@rnaCounts[rn,], obj@rnaCounts[obj@controls.forfit,]),
                      ddepth=obj@dnaDepth,
                      rdepth=obj@rnaDepth,
                      rctrlscale=obj@rnaCtrlScale,
                      ddesign.mat=obj@designs@dna,
                      rdesign.mat=obj@designs@rnaRed,
                      rdesign.ctrl.mat=obj@designs@rnaCtrlRed,
                      theta.d.ctrl.prefit=
                          do.call(rbind, lapply(obj@modelPreFits.dna.ctrl,
                                                function(x) x$d.coef)),
                      compute.hessian=FALSE))
    }, BPPARAM = obj@BPPARAM)
    names(obj@modelFits.red) <- rownames(obj@dnaCounts)
    
    ## run lrt
    obj@results <- test.lrt(obj)
    
    # TODO rank based on strength of effect
    return(obj)
}

#' @rdname analyse.casectrl
#' @export
analyse.quant <- function(obj, mode="quant", model=NULL, dnaDesign=NULL, rnaDesign=~1) {
    
    ## check depth is set
    if(length(obj@dnaDepth) == 0 | length(obj@rnaDepth) == 0) {
        stop("Library depth factors must be estimated or manually set")
    }
    ## check mode
    if(!mode %in% c("quant")) {
        stop("Only mode 'quant' is sensible if no control enhancers are supplied")
    }
    obj@mode <- mode
    obj@model <- model
    
    ## get design matrices
    obj@designs@dna <- getDesignMat(obj, dnaDesign)
    obj@designs@rnaFull <- getDesignMat(obj, rnaDesign)
    
    ## fit full models
    fitfun <- fit.dnarna.noctrlobs
    obj@modelFits <- bplapply(rownames(obj@dnaCounts), function(rn) {
        return(fitfun(model=model,
                      dcounts=obj@dnaCounts[rn,],
                      rcounts=obj@rnaCounts[rn,],
                      ddepth=obj@dnaDepth,
                      rdepth=obj@rnaDepth,
                      rctrlscale=NULL,
                      ddesign.mat=obj@designs@dna,
                      rdesign.mat=obj@designs@rnaFull,
                      rdesign.ctrl.mat=NULL,
                      theta.d.ctrl.prefit=NULL,
                      compute.hessian=FALSE))
    }, BPPARAM = obj@BPPARAM)
    
    # TODO add test against lower quantile here (zscore)
    
    # TODO rank based on strength of effect
}

#' get a design matrix from the input design
#' 
#' adds an extra first column if this is the rna design matrix and
#' the rna noise model requires a separate variance parameter which is not
#' used for the mean model, ie. if the first parameter to estimate is not 
#' part of the GLM that defines the mean parameter.
#'
#' @param obj the MpraObject
#' @param design the input design. A matrix is returned as is, a formula is
#' expanded to a model matrix using the object colAnnot. If design is NULL, an
#' intercept-only design matrix is created and returned
#' @param testcondition condition to substract from formula
#' @param rna whether this is a rna model. Together with the noise model set
#' in obj@model, this defines whether an extra column for the variance link
#' parameter is allowed in the design matrix.
#'
#' @return a design matrix
getDesignMat <- function(obj, design, testcondition=NULL, rna=FALSE) {
    # depreceate matrix entry?
    if (is.matrix(design)) {
        dmat <- design
    } else if (is.null(design)) {
        dmat <- matrix(rep(1,NCOL(obj@dnaCounts)), ncol = 1,
                       dimnames = list(colnames(obj@dnaCounts), "(intercept)"))
    } else if (is(design, "formula")) {
        if(!is.null(testcondition)){
            # substract this condition from formula
            terms <- attr(terms.formula(design), "term.labels")
            termsnew <- terms[terms != testcondition]
            if(length(termsnew) >= 1) {
                design <- as.formula(paste0("~", paste(termsnew, collapse="+")))
            } else {
                design <- ~1
            }
        }
        dmat <- model.matrix(design, obj@colAnnot)
    } else {
        stop("invalid design")
    }
    if(obj@model %in% c("ln.nb") & rna==TRUE) {
        # the first column with only zero entries allows
        # treating the parameter vector in the GLM desing matrix
        # matrix multiplication vector even though the first parameter
        # just defines the variance:
        # ln.nb: the first parameter is the dispersion parameter
        dmat <- as.matrix(cbind(variance_padding=0, dmat))
    }
    return(dmat)
}
