#' Fit the model for a likelihood-ratio based testing of the significance of the
#' effect of the given condition.
#' 
#' @export
#' 
#' @details comparative LRT-based analysis can be run in three possible modes.
#' \itemize{
#'   \item "nocontrol": this is the default mode for scenarios where control 
#'   sequences (scrambles) are not available. In this mode, The full model folows
#'   the given design, and the reduced model is similar but with the condition
#'   factor removed.
#'   \item "scaled": this is the default mode if controls sequences are available.
#'   In this mode, a separate model is fitted to the control sequences, and the
#'   fitten parameters are incorporated in the data model to account for unknown
#'   sources of unqanted variation. The control model is used to scale the data 
#'   model parameters, similar to library depth factors.
#'   \item "full": this mode answers a slightly different statistical question.
#'   instead of measuring the effect of the condition explicitly, this mode test
#'   if the sequence at hand behaves differently from the control enhancers. In
#'   this mode a separate model is fit for the control enhancers, and incorporated
#'   fully into the data model for each sequence in th full model, or separately
#'   for the reduced model. This mode would detect any sequence with statistically
#'   significant deviance from the control distribution.
#' }
#' 
#' @param obj the MPRAnalyze object
#' @param condition the condition to test. Must be a valid column name in the 
#' annotation of the object. Ignored in "full" mode.
#' @param mode the LRT mode to run in. See details
#' @param dnaDesign the design of the DNA model. See details
#' @param rnaDesign the design of the RNA model. See details
#' 
#' @return the MPRAnalyze object, populated with fitted models
analyze.comparative.lrt <- function(obj, condition=NULL, mode=NULL, 
                                  dnaDesign=NULL, rnaDesign= ~condition){
    if(length(obj@dnaDepth) == 0){
        stop("library depth factors must be estimated or provided before analysis")
    }
    if(length(obj@model) == 0) {
        obj <- autoChooseModel(obj)
    }
    
    ## get LRT mode
    if(is.null(mode)) {
        if(!is.null(obj@controls)) {
            mode <- "scaled"
        } else {
            mode <- "nocontrols"
        }
    } else if (!(mode %in% c("nocontrols", "scaled", "full"))){
        stop("mode ", mode, " is not supported")
    } else if (!(mode %in% c("full")) & is.null(condition)) {
        stop("condition must be provided for mode ", mode)
    } else if (!(mode %in% c("nocontrols")) & is.null(obj@controls)) {
        stop("mode ", mode, " is not supported if no controls are provided")
    }
    obj@mode <- paste0("comparative.lrt.", mode)
    
    obj@designs@dna <- getDesignMat(obj=obj, design=dnaDesign)
    obj@designs@rnaFull <- getDesignMat(obj=obj, design=rnaDesign)
    
    return(analyze.lrt(obj=obj, condition=condition, 
                       dnaDesign=dnaDesign, rnaDesign=rnaDesign))
}

#' TODO:
#' 
#' @export
analyze.comparative.ttest <- function(obj, dnaDesign, rnaDesign, condition){
    if(length(obj@model) == 0) {
        obj <- autoChooseModel(obj)
    }
    ## check if there is an intercept or not
    ## make sure condition is first term in formula
    ## mode: comparative.ttest.1ref (has intercept, coeff = offset)
    ## mode: comparative.ttest.mrefs (no intercept, all pairs possible)
    ## fit model
    ## extract slopes, reformat for results data.frame
    ## TODO: test.ttest:
    ##    if .1ref: get one argmuent, compute ttest
    ##    if .mrefs: get two arguments, compute ttest
    ##    return results and append them to the results data.frame
}

#' Perform quantitative analysis on the MPRA data. This analysis aims to determine
#' which sequences have a regulatory function, when no condition is being tested.
#' 
#' @details quantitative analysis can be performed in several ways, depending on
#' the experimental design and user preference:
#' \itemize{
#'   \item empirical: default mode if negative controls are provided. After the
#'   fit, the 'slope' factor is extracted from the model (this is the factor 
#'   that captures transcriptional rate). An empricial pvalue is computed for 
#'   each candidate based on the distribution of the control enhancers.
#'   \item lrt: only available if negative controls are provided. A likelihood
#'   ratio test is used, with the null hypothesis a joint model of the controls
#'   and a given candidate sequence, and the alternative model being a joint model
#'   with an additional factor allowing for a separation of controls and the candidate.
#'   \item idr: default mode if no controls are available, but multiple libraries
#'   exist in the experiment. In this mode, the library factor used for depth
#'   estimation is used. We fit the model, then the library-specific slope is
#'   extracted from the model and a ranking vector of the slopes is computed.
#'   The Irreproducible-Discovery-Rate (IDR) method is then used to estimate 
#'   significance.
#'   \item globalScale: this is the default method for scenarios in which the 
#'   experiment has a single library and no controls (this is also the only
#'   method that supports that scenario). TODO
#' }
#' 
#' @param obj the MpraObject
#' @param mode the analysis mode, see details
#' @param dnaDesign the design of the DNA counts
#' @param rnaDesign the design of the RNA counts
#' 
#' @return the MpraObject, with populated models and results
#' 
#' @export
analyze.quantitative <- function(obj, mode=NULL, dnaDesign=~1, rnaDesign=~1){
    if(length(obj@dnaDepth) == 0){
        stop("library depth factors must be estimated or provided before analysis")
    }
    if(length(obj@model) == 0) {
        obj <- autoChooseModel(obj)
    }
    if(is.null(mode)) {
        if(is.null(obj@controls)) {
            if(length(levels(obj@lib.factor)) > 1) {
                mode <- "IDR"
            } else {
                mode <- "globalScale"
            }
        } else {
            mode <- "empirical"
        }
    } else if (!(mode %in% c("IDR", "lrt", "globalScale", "empirical"))) {
        stop("mode ", mode, " is not supported")
    }
    
    obj@mode <- paste0("quantitative.", mode)
    
    obj@designs@dna <- getDesignMat(obj=obj, design=dnaDesign)
    obj@designs@rnaFull <- getDesignMat(obj=obj, design=rnaDesign)
    
    return(QUANT_ANALYSIS[[obj@mode]](obj, dnaDesign, rnaDesign))
}

analyze.quantitative.emppval <- function(obj, dnaDesign=NULL, rnaDesign=NULL){
    ## validate
    if(!checkForIntercept(obj@designs@rnaFull)) {
        stop("Design matrix must have an intercept in empirical mode")
    }

    ## fit model
    obj@designs@rnaRed <- NULL
    obj@designs@rnaCtrlFull <- NULL
    obj@designs@rnaCtrlRed <- NULL
    obj@modelPreFits.dna.ctrl <- NULL
    obj@controls.forfit <- NULL
    
    message("Fitting Model...")    
    obj@modelFits <- bplapply(rownames(obj@dnaCounts), function(rn) {
        return(fit.dnarna.noctrlobs(model=obj@model,
                      dcounts=obj@dnaCounts[rn,,drop=FALSE],
                      rcounts=obj@rnaCounts[rn,,drop=FALSE],
                      ddepth=obj@dnaDepth,
                      rdepth=obj@rnaDepth,
                      rctrlscale=NULL,
                      ddesign.mat=obj@designs@dna,
                      rdesign.mat=obj@designs@rnaFull,
                      rdesign.ctrl.mat=NULL,
                      theta.d.ctrl.prefit=NULL,
                      compute.hessian=FALSE))
    }, BPPARAM = obj@BPPARAM)
    names(obj@modelFits) <- rownames(obj@dnaCounts)
    
    ## extract slope - second coefficient of the rna model
    slopes <- vapply(obj@modelFits, function(x) x$r.coef[2], 0.0)
    names(slopes) <- names(obj@modelFits)
    
    ## create ecdf and compute epvalues
    ctrls <- slopes[obj@controls]
    ctrl.cdf <- ecdf(ctrls)
    epval <- 1 - ctrl.cdf(slopes)
    obj@results <- data.frame(row.names=names(slopes),
                              coef=slopes,
                              pval=epval,
                              fdr=p.adjust(epval, 'BH'))
    
    return(obj)
}

analyze.quantitative.globalscale <- function(obj, dnaDesign=NULL, rnaDesign=NULL){
    ## fit model
    ## take median of slopes as global scaler
    ## look for deviance from 1
    ## TODO: look at DESeq2 to see what exactly it is they do
}

analyze.quantitative.idr <- function(obj, dnaDesign=NULL, rnaDesign=NULL){
    ## fit model
    ## extract library-specific slopes
    ## feed rankings to idr and get significance estimates
}

#' run LRT-based analysis
#' @param obj the MPRAnalyze object
#' @param condition the condition to test. Must be a valid column name in the 
#' annotation of the object. Ignored in "full" mode.
#' @param mode the LRT mode to run in. See details
#' @param dnaDesign the design of the DNA model. See details
#' @param rnaDesign the design of the RNA model. See details
#' 
#' @return the MPRAnalyze object, populated with fitted models
analyze.lrt <- function(obj, condition=NULL, 
                        dnaDesign=NULL, rnaDesign= ~condition) {
    ##TODO: verify mode quantitative.lrt!!! works
    
    ## get distributional model
    if(is.null(obj@model)) {
        obj@model <- autoChooseModel(obj)
    }
    
    if(obj@mode == "comparative.lrt.nocontrols") {
        obj@designs@rnaRed <- getDesignMat(obj, rnaDesign, condition)
        obj@designs@rnaCtrlFull <- NULL
        obj@designs@rnaCtrlRed <- NULL
        obj@modelPreFits.dna.ctrl <- NULL
        obj@controls.forfit <- NULL
        fitfun <- fit.dnarna.noctrlobs
    } else {
        message("Fitting controls-based background model...")
        obj@modelPreFits.dna.ctrl <- fit.dnarna.onlyctrl.iter(
            model=obj@model,
            dcounts = obj@dnaCounts[obj@controls,],
            rcounts = obj@rnaCounts[obj@controls,],
            ddepth=obj@dnaDepth,
            rdepth=obj@rnaDepth,
            ddesign.mat=obj@designs@dna,
            rdesign.mat=obj@designs@rnaFull,
            BPPARAM = obj@BPPARAM)
        if(obj@mode == "comparative.lrt.scaled") {
            obj@designs@rnaRed <- getDesignMat(obj, rnaDesign, condition)
            obj@designs@rnaCtrlFull <- obj@designs@rnaFull
            obj@designs@rnaCtrlRed <- obj@designs@rnaFull
            obj@rnaCtrlScale <- obj@modelPreFits.dna.ctrl[[1]]$r.coef
            obj@controls.forfit <- NULL
            fitfun <- fit.dnarna.noctrlobs
        } else { ## comparative.lrt.full / quantitative.lrt ?
            obj@designs@rnaRed <- obj@designs@rnaFull
            obj@designs@rnaCtrlFull <- obj@designs@rnaFull
            obj@designs@rnaCtrlRed <- NULL
            obj@rnaCtrlScale <- NULL
            obj@controls.forfit <- obj@controls
            fitfun <- fit.dnarna.wctrlobs.iter
        }
    }
    
    theta.d.ctrl.prefit <- do.call(cbind, lapply(obj@modelPreFits.dna.ctrl,
                                                 function(x) x$d.coef))
    
    message("Fitting full model...")
    obj@modelFits <- bplapply(rownames(obj@dnaCounts), function(rn) {
        return(fitfun(model=obj@model,
                      dcounts=obj@dnaCounts[rn,,drop=FALSE],
                      rcounts=rbind(obj@rnaCounts[rn,], 
                                    obj@rnaCounts[obj@controls.forfit,]),
                      ddepth=obj@dnaDepth,
                      rdepth=obj@rnaDepth,
                      rctrlscale=obj@rnaCtrlScale,
                      ddesign.mat=obj@designs@dna,
                      rdesign.mat=obj@designs@rnaFull,
                      rdesign.ctrl.mat=obj@designs@rnaCtrlFull,
                      theta.d.ctrl.prefit=theta.d.ctrl.prefit,
                      compute.hessian=FALSE))
    }, BPPARAM = obj@BPPARAM)
    names(obj@modelFits) <- rownames(obj@dnaCounts)
    
    message("Fitting reduced model...")
    obj@modelFits.red <- bplapply(rownames(obj@dnaCounts), function(rn) {
        return(fitfun(model=obj@model,
                      dcounts=obj@dnaCounts[rn,,drop=FALSE],
                      rcounts=rbind(obj@rnaCounts[rn,], 
                                    obj@rnaCounts[obj@controls.forfit,]),
                      ddepth=obj@dnaDepth,
                      rdepth=obj@rnaDepth,
                      rctrlscale=obj@rnaCtrlScale,
                      ddesign.mat=obj@designs@dna,
                      rdesign.mat=obj@designs@rnaRed,
                      rdesign.ctrl.mat=obj@designs@rnaCtrlRed,
                      theta.d.ctrl.prefit=theta.d.ctrl.prefit,
                      compute.hessian=FALSE))
    }, BPPARAM = obj@BPPARAM)
    names(obj@modelFits.red) <- rownames(obj@dnaCounts)
    
    message("Computing statistical test...")
    obj@results <- test.lrt(obj)
    
    return(obj)
}

QUANT_ANALYSIS <- list(quantitative.IDR = analyze.quantitative.idr,
                       quantitative.lrt = analyze.lrt,
                       quantitative.globalScale = analyze.quantitative.globalscale,
                       quantitative.empirical = analyze.quantitative.emppval)
