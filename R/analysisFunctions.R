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

#' Fit the model to enable coefficient-based testing (using the test.coef).
#' This method should be used when multiple conditions in the data are to be
#' tested using the same model.
#' For example, various stimulations compared with unstimulated.
#' 
#' @note Currently, this mode only supports a single reference (intercept term).
#' 
#' @param obj the MpraObject
#' @param dnaDesign the design of the DNA model
#' @param rnaDesign the design of the RNA model
#' @param use.controls if the experiment has negative controls, thes can be 
#' included in the model and used to correct for unwanted variation. Default is 
#' TRUE, and this is ignored if no controls exist in the experiment.
#' 
#' @return the MpraObject with the fitted model, that can be passed to test.coef 
#' to test for significance of any coefficient in the model.
#' 
#' @export
analyze.comparative.coef <- function(obj, dnaDesign, rnaDesign, 
                                     use.controls=TRUE){
    if(length(obj@dnaDepth) == 0) {
        stop("library depth factors must be estimated or provided before analysis")
    }
    if(!checkForIntercept(rnaDesign)) {
        stop("only designs with an intercept term are currently supported in this mode")
    }
    
    if(length(obj@model) == 0) {
        obj <- autoChooseModel(obj)
    }
    
    obj@designs@dna <- getDesignMat(obj=obj, design=dnaDesign)
    obj@designs@rnaFull <- getDesignMat(obj=obj, design=rnaDesign)
    obj@mode <- "comparative.coef"
    
    if(use.controls & !is.null(obj@controls)) {
        message("Fitting controls-based background model...")
        obj@modelPreFits.dna.ctrl <- reformatModels(fit.dnarna.onlyctrl.iter(
            model=obj@model,
            dcounts = obj@dnaCounts[obj@controls,],
            rcounts = obj@rnaCounts[obj@controls,],
            ddepth=obj@dnaDepth,
            rdepth=obj@rnaDepth,
            ddesign.mat=obj@designs@dna,
            rdesign.mat=obj@designs@rnaFull,
            BPPARAM = obj@BPPARAM))
        
        obj@designs@rnaCtrlFull <- obj@designs@rnaFull
        obj@designs@rnaCtrlRed <- obj@designs@rnaFull
        obj@rnaCtrlScale <- obj@modelPreFits.dna.ctrl$r.coef[1,]
        obj@controls.forfit <- NULL
        fitfun <- fit.dnarna.noctrlobs
        obj@rnaCtrlScale <- obj@modelPreFits.dna.ctrl$r.coef[1,]
        if(!is.null(obj@modelPreFits.dna.ctrl)) {
            theta.d.ctrl.prefit <- t(obj@modelPreFits.dna.ctrl$d.coef)
        } else {
            theta.d.ctrl.prefit <- NULL
        }
    } else {
        obj@designs@rnaCtrlFull <- NULL
        obj@designs@rnaCtrlRed <- NULL
        obj@modelPreFits.dna.ctrl <- NULL
        obj@controls.forfit <- NULL
        theta.d.ctrl.prefit <- NULL
    }
    ## fit model
    message("Fitting model...")
    models <- bplapply(rownames(obj@dnaCounts), function(rn) {
        return(fit.dnarna.noctrlobs(model=obj@model,
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
                      compute.hessian=TRUE))
    }, BPPARAM = obj@BPPARAM)
    names(models) <- rownames(obj@dnaCounts)
    obj@modelFits <-reformatModels(models)
    
    message("Analysis done!")
    return(obj)
}

#' Perform quantitative analysis on the MPRA data. This analysis aims to determine
#' which sequences have a regulatory function, when no condition is being tested.
#' 
#' @details quantitative analysis can be performed in several ways, depending on
#' the experimental design and user preference. This decision effects the model
#' fitting and may cause certain types of hypothesis testing to be enabled or 
#' disabled.
#' \itemize{
#'   \item epirical: the model is fitted as specified, enabling future empirical
#'   testing (either empirical p-value if negative controls are provided, or a
#'   global devience analysis, see details in `test.empirical`)
#'   \item lrt: only available if negative controls are provided. A likelihood
#'   ratio test is used, with the null hypothesis a joint model of the controls
#'   and a given candidate sequence, and the alternative model being a separate
#'   model for controls and candidates.
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
        mode <- "empirical"
    } else if (!(mode %in% c("lrt", "empirical"))) {
        stop("mode ", mode, " is not supported")
    }
    
    obj@mode <- paste0("quantitative.", mode)
    
    obj@designs@dna <- getDesignMat(obj=obj, design=dnaDesign)
    obj@designs@rnaFull <- getDesignMat(obj=obj, design=rnaDesign)

    return(QUANT_ANALYSIS[[obj@mode]](obj, dnaDesign, rnaDesign))
}

#' Fit the model for quantitative analysis
#' @param obj the MpraObject
#' @param dnaDesign the design of the DNA model
#' @param rnaDesign the design of the RNA model
#' @return the MpraObject with fitted models
analyze.quantitative.empirical <- function(obj, dnaDesign=NULL, rnaDesign=NULL){
    ## fit model
    obj@designs@rnaRed <- NULL
    obj@designs@rnaCtrlFull <- NULL
    obj@designs@rnaCtrlRed <- NULL
    obj@modelPreFits.dna.ctrl <- NULL
    obj@controls.forfit <- NULL
    
    message("Fitting model...")
    models <- bplapply(rownames(obj@dnaCounts), function(rn) {
        tryCatch({
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
        }, error = function(err) {message("error fitting: ", rn)})
    }, BPPARAM = obj@BPPARAM)
    names(models) <- rownames(obj@dnaCounts)
    obj@modelFits <- reformatModels(models)
    
    message("Analysis done!")
    return(obj)
}


#' run LRT-based analysis
#' @param obj the MPRAnalyze object
#' @param condition the condition to test. Must be a valid column name in the 
#' annotation of the object. Ignored in "full" mode.
#' @param dnaDesign the design of the DNA model. See details
#' @param rnaDesign the design of the RNA model. See details
#' 
#' @return the MPRAnalyze object, populated with fitted models
analyze.lrt <- function(obj, condition=NULL, 
                        dnaDesign=NULL, rnaDesign= ~condition) {
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
        obj@modelPreFits.dna.ctrl <- reformatModels(
                                     fit.dnarna.onlyctrl.iter(
            model=obj@model,
            dcounts = obj@dnaCounts[obj@controls,],
            rcounts = obj@rnaCounts[obj@controls,],
            ddepth=obj@dnaDepth,
            rdepth=obj@rnaDepth,
            ddesign.mat=obj@designs@dna,
            rdesign.mat=obj@designs@rnaFull,
            BPPARAM = obj@BPPARAM))
        if(obj@mode == "comparative.lrt.scaled") {
            obj@designs@rnaRed <- getDesignMat(obj, rnaDesign, condition)
            obj@designs@rnaCtrlFull <- obj@designs@rnaFull
            obj@designs@rnaCtrlRed <- obj@designs@rnaFull
            obj@rnaCtrlScale <- obj@modelPreFits.dna.ctrl$r.coef[1,]
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
    
    if(!is.null(obj@modelPreFits.dna.ctrl)) {
        theta.d.ctrl.prefit <- t(obj@modelPreFits.dna.ctrl$d.coef)
    } else {
        theta.d.ctrl.prefit <- NULL
    }
    
    
    message("Fitting model...")
    models <- bplapply(rownames(obj@dnaCounts), function(rn) {
        tryCatch({
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
        }, error = function(err) {message("error fitting: ", rn)})
    }, BPPARAM = obj@BPPARAM)
    names(models) <- rownames(obj@dnaCounts)
    obj@modelFits <- reformatModels(models)
    
    message("Fitting reduced model...")
    models <- bplapply(rownames(obj@dnaCounts), function(rn) {
        tryCatch({
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
        }, error = function(err) {message("error fitting: ", rn)})
    }, BPPARAM = obj@BPPARAM)
    names(models) <- rownames(obj@dnaCounts)
    obj@modelFits.red <- reformatModels(models)
    
    message("Analysis done!")
    return(obj)
}

analyze.quantitative.lrt <- function(obj, dnaDesign=NULL, rnaDesign=NULL) {
    ## Fit Full model: DNARNA joint fit per enhancer
    message("Fitting full model...")    
    models <- bplapply(rownames(obj@dnaCounts), function(rn) {
        tryCatch({
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
        }, error = function(err) {message("error fitting: ", rn)})
    })
    names(models) <- rownames(obj@dnaCounts)
    obj@modelFits <- reformatModels(models)
    
    ## fit control model
    message("Fitting controls-based model...")
    obj@modelPreFits.dna.ctrl <- reformatModels(fit.dnarna.onlyctrl.iter(
        model=obj@model,
        dcounts = obj@dnaCounts[obj@controls,],
        rcounts = obj@rnaCounts[obj@controls,],
        ddepth=obj@dnaDepth,
        rdepth=obj@rnaDepth,
        ddesign.mat=obj@designs@dna,
        rdesign.mat=obj@designs@rnaFull,
        BPPARAM = obj@BPPARAM))
    
    ctrl.rna <- obj@modelPreFits.dna.ctrl$r.coef[1,]
    ## TODO: ctrl models don't save the dispersion parameter, so padding is needed
    ctrl.rna <- c(0, ctrl.rna)
    
    ## fit Reduced model: DNA per-enhancer, conditioned on control RNA
    message("Fitting reduced model...")
    models <- bplapply(rownames(obj@dnaCounts), function(rn) {
        tryCatch({
        return(fit.dna.controlrna(model = obj@model,
                                  dcounts = obj@dnaCounts[rn,,drop=FALSE],
                                  rcounts = obj@rnaCounts[rn,,drop=FALSE],
                                  r.coef = ctrl.rna,
                                  ddepth=obj@dnaDepth,
                                  rdepth=obj@rnaDepth,
                                  ddesign.mat=obj@designs@dna,
                                  rdesign.mat=obj@designs@rnaFull))
        }, error = function(err) {message("error fitting: ", rn)})
    }, BPPARAM = obj@BPPARAM)
    names(models) <- rownames(obj@dnaCounts)
    obj@modelFits.red <- reformatModels(models)
    obj@modelFits.red$r.df <- 0
    
    message("Analysis done!")
    return(obj)
}

QUANT_ANALYSIS <- list(quantitative.lrt = analyze.quantitative.lrt,
                       quantitative.empirical = analyze.quantitative.empirical)
