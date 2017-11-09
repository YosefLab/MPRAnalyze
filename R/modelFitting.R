
#' easy-access list relating model id character to the fitting wrappers
#' @include optimizers.R
FITTING_FUNCTIONS <- list(
    gammaPoisson = fit.mixture.gammaPoisson,
    lnDNA.nbRNA = function(...) {
        fit.separate(dnaFn=fit.lnDNA, rnaFn=fit.nbRNA, ...)
        },
    gammaDNA.nbRNA = function(...) {
        fit.separate(dnaFn=fit.gammaDNA, rnaFn=fit.nbRNA, ...)
    })


#' fit models for a differential activity analysis. This function does not compute
#' the test tatistics, only performs the initial model fitting necessary.
#'
#' @param obj the MpraObject
#' @param model the model to use (TODO: if nul?)
#' @param dnaDesign the design to use for the DNA model (see details)
#' @param rnaDesign the design to use for the RNA model (see details)
#' @param mode options are "LRT" (default") or "T-Test", see details
#'
#' @details TODO - explain input designs (either formula or model matrices)
#'
#' @details TODO - explain difference between fitting for T-test and fitting for
#' LRT, refer to vignette for additional details
#'
#' @return the MpraObject with fitted models
fit.differential <- function(obj, model=NULL, dnaDesign=NULL, rnaDesign=NULL,
                            mode="LRT") {
    ## check mode
    if(mode=="LRT") {
        obj@mode = "LRT"
        compHess = FALSE
    } else if (mode=="T-Test") {
        obj@mode = "T-Test"
        compHess = TRUE
    } else {
        stop("Only 'LRT' and 'T-Test' are currently supported")
    }

    obj <- estimateDepthFactors(obj)

    ## get design matrices
    obj@designs@dnaFull <- getDesignMat(obj, dnaDesign)
    obj@designs@rnaFull <- getDesignMat(obj, rnaDesign)

    ## get fitting function
    if(is.null(model)) {
        model <- autoChooseModel(obj)
    }
    fitfun <- FITTING_FUNCTIONS[[model]]
    if(is.null(fitfun)) {
        stop(paste0(model, " is not a supported model"))
    }
    obj@model <- model

    ## fit models
    obj@modelFits <- bplapply(rownames(obj@dnaCounts), function(rn) {
        return(fitfun(dcounts=obj@dnaCounts[rn,],
                      rcounts=obj@rnaCounts[rn,],
                      ddepth=obj@dnaDepth,
                      rdepth=obj@rnaDepth,
                      ddesign.mat=obj@designs@dnaFull,
                      rdesign.mat=obj@designs@rnaFull,
                      compute.hessian=compHess))
    }, BPPARAM = obj@BPPARAM)

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
#'
#' @return a design matrix
getDesignMat <- function(obj, design) {
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
