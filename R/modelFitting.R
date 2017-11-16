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
    obj@designs@rnaFull <- getDesignMat(obj, rnaDesign)
    if(is.null(obj@controls)) {
        obj@designs@rnaRed <- getDesignMat(obj, rnaDesign, condition_totest)
        obj@designs@rnaCtrlFull <- NULL
        obj@designs@rnaCtrlRed <- NULL
        obj@modelPreFits.dna.ctrl <- NULL
        fitfun <- fit.dnarna.noctrlobs
    } else {
        message("Fit control enhancer background models")
        obj@designs@rnaRed <- obj@designs@rnaFull
        obj@designs@rnaCtrlFull <- getDesignMat(obj, rnaDesign, condition_totest) #?
        obj@designs@rnaCtrlRed <- NULL
        obj@modelPreFits.dna.ctrl <- fit.dnarna.onlyctrl.iter(
            model=model,
            dcounts = obj@dnaCounts[obj@controls,], 
            rcounts = obj@rnaCounts[obj@controls,],
            ddepth=obj@dnaDepth,
            rdepth=obj@rnaDepth,
            ddesign.mat=obj@designs@dna,
            rdesign.mat=obj@designs@rnaFull)
        if(obj@mode == "scaled") {
            obj@rnaCtrlScale <- obj@modelPreFits.dna.ctrl[[1]]$r.par
            fitfun <- fit.dnarna.noctrlobs
        } else if(obj@mode == "full") {
            obj@rnaCtrlScale <- NULL
            fitfun <- fit.dnarna.wctrlobs.iter
        }
    }
    
    ## fit full models
    message("Fit full models")
    obj@modelFits <- lapply(rownames(obj@dnaCounts), function(rn) {
        return(fitfun(model=model,
                      dcounts=obj@dnaCounts[rn,,drop=FALSE],
                      rcounts=obj@rnaCounts[rn,,drop=FALSE],
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
    })#, BPPARAM = obj@BPPARAM)
    
    ## fit reduced models
    message("Fit reduced models")
    obj@modelFits.red <- bplapply(rownames(obj@dnaCounts), function(rn) {
        return(fitfun(model=model,
                      dcounts=obj@dnaCounts[rn,,drop=FALSE],
                      rcounts=obj@rnaCounts[rn,,drop=FALSE],
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
        obj@designs@rnaCtrlFull <- NULL
        obj@modelPreFits.dna.ctrl <- NULL
        fitfun <- fit.dnarna.noctrlobs
    } else {
        obj@designs@rnaCtrlFull <- getDesignMat(obj, rnaDesign, condition_totest) #?
        obj@modelPreFits.dna.ctrl <- fit.dnarna.onlyctrl.iter(
            model=model,
            dcounts = obj@dnaCounts[obj@controls,], 
            rcounts = obj@rnaCounts[obj@controls,],
            ddepth=obj@dnaDepth,
            rdepth=obj@rnaDepth,
            ddesign.mat=obj@designs@dna,
            rdesign.mat=obj@designs@rnaFull)
        if(obj@mode == "scaled") {
            obj@rnaCtrlScale <- obj@modelPreFits.dna.ctrl[[1]]$r.par
            fitfun <- fit.dnarna.noctrlobs
        } else if(obj@mode == "full") {
            obj@rnaCtrlScale <- NULL
            fitfun <- fit.dnarna.wctrlobs.iter
        }
    }
    
    ## fit full models
    obj@modelFits <- bplapply(rownames(obj@dnaCounts), function(rn) {
        return(fitfun(model=model,
                      dcounts=obj@dnaCounts[rn,],
                      rcounts=obj@rnaCounts[rn,],
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
        stop("nothing to test against")
    } else {
        obj@designs@rnaRed <- NULL
        obj@designs@rnaCtrlFull <- obj@designs@rnaFull
        obj@designs@rnaCtrlRed <- NULL
        obj@modelPreFits.dna.ctrl <- fit.dnarna.onlyctrl.iter(
            model=model,
            dcounts = obj@dnaCounts[obj@controls,], 
            rcounts = obj@rnaCounts[obj@controls,],
            ddepth=obj@dnaDepth,
            rdepth=obj@rnaDepth,
            ddesign.mat=obj@designs@dna,
            rdesign.mat=obj@designs@rnaFull)
        if(obj@mode == "scaled") {
            obj@rnaCtrlScale <- obj@modelPreFits.dna.ctrl[[1]]$r.par
            fitfun <- fit.dnarna.noctrlobs
        } else if(obj@mode == "full") {
            obj@rnaCtrlScale <- NULL
            fitfun <- fit.dnarna.wctrlobs.iter
        }
    }
    
    ## fit full models
    obj@modelFits <- bplapply(rownames(obj@dnaCounts), function(rn) {
        return(fitfun(model=model,
                      dcounts=obj@dnaCounts[rn,],
                      rcounts=obj@rnaCounts[rn,],
                      ddepth=obj@dnaDepth,
                      rdepth=obj@rnaDepth,
                      rctrlscale=obj@rnaCtrlScale,
                      ddesign.mat=obj@designs@dna,
                      rdesign.mat=obj@designs@rnaFull,
                      rdesign.ctrl.mat=obj@designs@rnaFull,
                      theta.d.ctrl.prefit=
                          do.call(rbind, lapply(obj@modelPreFits.dna.ctrl,
                                                function(x) x$d.coef)),
                      compute.hessian=FALSE))
    }, BPPARAM = obj@BPPARAM)
    
    ## fit reduced models
    obj@modelFits.red <- bplapply(rownames(obj@dnaCounts), function(rn) {
        return(fitfun(model=model,
                      dcounts=obj@dnaCounts[rn,],
                      rcounts=obj@rnaCounts[rn,],
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
    
    ## run lrt
    obj@results <- test.lrt(obj)
    
    # TODO rank based on strength of effect
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
#' @param obj the MpraObject
#' @param design the input design. A matrix is returned as is, a formula is
#' expanded to a model matrix using the object colAnnot. If design is NULL, an
#' intercept-only design matrix is created and returned
#' @param testcondition condition to substract from formulae
#'
#' @return a design matrix
getDesignMat <- function(obj, design, testcondition=NULL) {
    if (is.matrix(design)) {
        return(design)
    } else if (is.null(design)) {
        return(matrix(rep(1,NCOL(obj@dnaCounts)), ncol = 1,
                      dimnames = list(colnames(obj@dnaCounts), "(intercept)")))
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
        return(model.matrix(design, obj@colAnnot))
    } else {
        stop("invalid design")
    }
}
