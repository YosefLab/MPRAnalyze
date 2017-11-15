#' easy-access list relating model id character to the fitting wrappers
#' @include optimizers.R
FITTING_FUNCTIONS <- list(
    "gamma.pois" = fit.dnarna.gammaPoisson,
    "lnDNA.nbRNA" = function(...) {
        fit.separate(dnaFn=fit.lnDNA, rnaFn=fit.nbRNA, ...)
    },
    "gammaDNA.nbRNA" = function(...) {
        fit.separate(dnaFn=fit.gammaDNA, rnaFn=fit.nbRNA, ...)
    })

#' fit models for a differential activity analysis. 
#' 
#' This function fit both full and reduced model and 
#' performs the hypothesis test.
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
#' @detils TODO rna control usage
#'
#' @details TODO - explain mode
#'
#' @return the MpraObject with fitted models
fit.differential.lrt <- function(obj, model="gamma.pois", mode=NULL
                                 dnaDesign=NULL, rnaDesign=NULL, condition_totest=NULL, ) {
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
    
    obj <- estimateDepthFactors(obj)
    
    ## get design matrices
    obj@designs@dna <- getDesignMat(obj, dnaDesign)
    obj@designs@rnaFull <- getDesignMat(obj, rnaDesign)
    if(is.null(obj@controls)) {
        obj@designs@rnaRed <- getDesignMat(obj, rnaDesign, condition_totest)
        obj@designs@rnaCtrlFull <- NULL
        obj@modelPreFits.dna.ctrl <- NULL
        fitfun <- fit.dnarna.noctrlobs
    } else {
        obj@designs@rnaRed <- obj@designs@rnaFull
        obj@designs@rnaCtrlFull <- getDesignMat(obj, condition_totest) #?
        if(obj@mode == "scaled") {
            obj@rnaCtrlScale <- prefit.rctrlscale(obj)
            obj@modelPreFits.dna.ctrl <- NULL
            fitfun <- fit.dnarna.noctrlobs
        } else if(obj@mode == "full") {
            obj@rnaCtrlScale <- NULL
            obj@modelPreFits.dna.ctrl <- fit.dnarna.onlyctrl.iter(
                model=model,
                dcounts = obj@dnaCounts[obj@controls,], 
                rcounts = obj@rnaCounts[obj@controls,],
                ddepth=obj@dnaDepth,
                rdepth=obj@rnaDepth,
                ddesign.mat=obj@designs@dna,
                rdesign.mat=obj@designs@rnaFull)
            fitfun <- fit.dnarna.wctrlobs.iter
        }
    }
    
    obj@model <- model
    
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
                          do.call(cbind, lappyl(obj@modelPreFits.dna.ctrl,
                                                function(x) x$d.coef)),
                      compute.hessian=FALSE))
    }, BPPARAM = obj@BPPARAM)
    
    ## fit reduced models
    obj@modelFits <- bplapply(rownames(obj@dnaCounts), function(rn) {
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
                          do.call(cbind, lappyl(obj@modelPreFits.dna.ctrl,
                                                function(x) x$d.coef)),
                      compute.hessian=FALSE))
    }, BPPARAM = obj@BPPARAM)
    
    ## run lrt
    obj@modelFits <- test.lrt(obj)
    
    return(obj)
}

#' fit model for quantitative activity analysis
#' @param obj the MpraObject
#' @param model the model to fit
#' @param dnaDesign the design for the DNA model
#' @param rnaDesign the design for the RNA model
fit.quantitative <- function(obj, model=NULL, dnaDesign=NULL, rnaDesign=NULL) {
    obj@mode = "Quant"
    ##TODO: fit a single model per enhancer
}

#' get a design matrix from the input design
#'
#' @param obj the MpraObject
#' @param design the input design. A matrix is returned as is, a formula is
#' expanded to a model matrix using the object colAnnot. If design is NULL, an
#' intercept-only design matrix is created and returned
#' @param testcondition condition to substract from formulae TODO
#'
#' @return a design matrix
getDesignMat <- function(obj, design, testcondition=NULL) {
    if (is.matrix(design)) {
        return(design)
    } else if (is.null(design)) {
        return(matrix(rep(1,NCOL(obj@dnaCounts)), ncol = 1,
                      dimnames = list(colnames(obj@dnaCounts), "(intercept)")))
    } else if (is(design, "formula")) {
        return(model.matrix(design, obj@colAnnot))
    } else {
        stop("invalid design")
    }
}
