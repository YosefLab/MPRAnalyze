#' Run a comparative analysis between conditions
#' @param obj the MpraObject
#' @param dnaDesign the design for the DNA model. Only terms that are matched
#' with the RNA design should be included.
#' @param rnaDesign the design for the RNA model.
#' @param fit.se logical, if TRUE the standard errors of the coefficients
#' are extracted from the model. These are necessary for computing coefficient-
#' based testing, but make the model fitting slower. Deafult: FALSE
#' @param reducedDesign the design for the reduced RNA model, for a likelihood-
#' ratio testing scheme. The Reduced design must be nested within the full 
#' design (i.e all terms in the reduced must be included in the full).
#' @param correctControls if TRUE (default), use the negative controls to 
#' establish the null hypothesis, correcting for systemic bias in the data
#' @param verbose print progress reports (default: TRUE)
#' @import progress
#' @export
#' @return the MpraObject with fitted models for the input enhancers
#' @examples
#' data <- simulateMPRA(tr = rep(2,10), da=c(rep(2,5), rep(2.5,5)), 
#'                      nbatch=2, nbc=20)
#' obj <- MpraObject(dnaCounts = data$obs.dna, 
#'                   rnaCounts = data$obs.rna, 
#'                   colAnnot = data$annot)
#' obj <- estimateDepthFactors(obj, lib.factor = "batch", which.lib = "both")
#' ## run an LRT-based analysis, as recommnded:
#' obj <- analyzeComparative(obj, dnaDesign = ~ batch + barcode + condition, 
#'                               rnaDesign = ~ condition, reducedDesign = ~ 1)
#'                               
#' ## alternatively, run a coefficient-based analysis:
#' obj <- analyzeComparative(obj, dnaDesign = ~ batch + barcode + condition, 
#'                               rnaDesign = ~ condition, fit.se = TRUE)
analyzeComparative <- function(obj, dnaDesign, rnaDesign, fit.se=FALSE, 
                                reducedDesign=NULL, correctControls=TRUE, 
                                verbose=TRUE) {
    ##TODO: if full becomes operational, add 'fullLRT' as a logical argument
    
    if(!fit.se & is.null(reducedDesign)) {
        stop("Comparative analysis requires either a reduced design or fitting \
             the SE")
    }
    if(!is.null(reducedDesign)) {
        if (!isNestedDesign(full=rnaDesign, reduced=reducedDesign)) {
            stop("reduced design must be nested within the full RNA design")
        }
    }
    if(length(obj@dnaDepth) != NCOL(obj@dnaCounts)) {
        obj <- estimateDepthFactors(obj, which.lib = "dna")
    }
    if(length(obj@rnaDepth) != NCOL(obj@rnaCounts)) {
        obj <- estimateDepthFactors(obj, which.lib = "rna")
    }
    if(length(obj@model) == 0) {
        obj <- autoChooseModel(obj)
    }
    
    obj@designs@dna <- getDesignMat(design=dnaDesign, annotations=obj@dnaAnnot)
    obj@designs@dna2rna <- getDesignMat(design=dnaDesign, 
                                        annotations=obj@rnaAnnot)
    obj@designs@rnaFull <- getDesignMat(design=rnaDesign, 
                                        annotations=obj@rnaAnnot)
    
    ##TODO: if full model, call a different function.
    
    ## if controls are to be used and fullModel not: fit the control model
    if(correctControls & all(!is.na(obj@controls))) {
        message("Fitting controls-based background model...")
        obj@modelPreFits.dna.ctrl <- reformatModels(fit.dnarna.onlyctrl.iter(
            model=obj@model,
            dcounts = obj@dnaCounts[obj@controls,],
            rcounts = obj@rnaCounts[obj@controls,],
            ddepth=obj@dnaDepth,
            rdepth=obj@rnaDepth,
            ddesign.mat=obj@designs@dna,
            rdesign.mat=obj@designs@rnaFull, 
            d2rdesign.mat=obj@designs@dna2rna,
            BPPARAM = obj@BPPARAM,
            print.progress = verbose))
        
        obj@designs@rnaCtrlFull <- obj@designs@rnaFull
        obj@designs@rnaCtrlRed <- obj@designs@rnaFull
        obj@rnaCtrlScale <- obj@modelPreFits.dna.ctrl$r.coef[1,]
        theta.d.ctrl.prefit <- t(obj@modelPreFits.dna.ctrl$d.coef)
    } else {
        theta.d.ctrl.prefit <- NULL
    }
    
    ## Fit the full model (with SE extraction if fit.SE is on)
    message("Fitting model...")
    pb <- progress_bar$new(format = "[:bar] :percent (:current/:total)", 
                        total = NROW(obj@dnaCounts), clear = !verbose)
    models <- bplapply(rownames(obj@dnaCounts), function(rn) {
        pb$tick()
        tryCatch({
            return(fit.dnarna.noctrlobs(
                model=obj@model,
                dcounts=obj@dnaCounts[rn,,drop=FALSE],
                rcounts=obj@rnaCounts[rn,,drop=FALSE],
                ddepth=obj@dnaDepth,
                rdepth=obj@rnaDepth,
                rctrlscale=obj@rnaCtrlScale,
                ddesign.mat=obj@designs@dna,
                rdesign.mat=obj@designs@rnaFull,
                d2rdesign.mat=obj@designs@dna2rna,
                rdesign.ctrl.mat=obj@designs@rnaCtrlFull,
                theta.d.ctrl.prefit=theta.d.ctrl.prefit,
                compute.hessian=fit.se))
        }, error = function(err) {message("error fitting ", rn, ": ", err)})
    }, BPPARAM = obj@BPPARAM)
    names(models) <- rownames(obj@dnaCounts)
    obj@modelFits <-reformatModels(models)
    
    if(!is.null(reducedDesign)) {
        obj@designs@rnaRed <- getDesignMat(design=reducedDesign, 
                                        annotations=obj@rnaAnnot)
        
        message("Fitting reduced model...")
        pb <- progress_bar$new(format = "[:bar] :percent (:current/:total)", 
                            total = NROW(obj@dnaCounts), clear = verbose)
        models <- bplapply(rownames(obj@dnaCounts), function(rn) {
            pb$tick()
            tryCatch({
                return(fit.dnarna.noctrlobs(
                    model=obj@model,
                    dcounts=obj@dnaCounts[rn,,drop=FALSE],
                    rcounts=obj@rnaCounts[rn,,drop=FALSE], 
                    ddepth=obj@dnaDepth,
                    rdepth=obj@rnaDepth,
                    rctrlscale=obj@rnaCtrlScale,
                    ddesign.mat=obj@designs@dna,
                    rdesign.mat=obj@designs@rnaRed,
                    d2rdesign.mat=obj@designs@dna2rna,
                    rdesign.ctrl.mat=obj@designs@rnaCtrlRed,
                    theta.d.ctrl.prefit=theta.d.ctrl.prefit,
                    compute.hessian=FALSE))
            }, error = function(err) {message("error fitting: ", 
                                              rn, ": ", err)})
        }, BPPARAM = obj@BPPARAM)
        names(models) <- rownames(obj@dnaCounts)
        obj@modelFits.red <- reformatModels(models)
    }
    
    message("Analysis Done!")
    return(obj)
}


#' Perform quantitative analysis on the MPRA data. This analysis aims to 
#' determine which sequences have a regulatory function, when no condition is 
#' being tested.
#' 
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
#' @param dnaDesign the design of the DNA counts
#' @param rnaDesign the design of the RNA counts
#' 
#' @return the MpraObject, with populated models
#' @import progress
#' @export
#' 
#' @examples
#' data <- simulateMPRA(tr = rep(2,10), da=NULL, nbatch=2, nbc=20)
#' obj <- MpraObject(dnaCounts = data$obs.dna, 
#'                   rnaCounts = data$obs.rna, 
#'                   colAnnot = data$annot)
#' obj <- estimateDepthFactors(obj, lib.factor = "batch", which.lib = "both")
#' obj <- analyzeQuantification(obj, dnaDesign = ~ batch + barcode, 
#'                               rnaDesign = ~1)
analyzeQuantification <- function(obj, dnaDesign=~1, rnaDesign=~1){
    if(length(obj@dnaDepth) == 0){
        warning("No DNA library depth factors set")
        obj <- estimateDepthFactors(obj, which.lib = "dna")
    }
    if(length(obj@rnaDepth) == 0){
        warning("No RNA library depth factors set")
        obj <- estimateDepthFactors(obj, which.lib = "rna")
    }
    if(length(obj@model) == 0) {
        obj <- autoChooseModel(obj)
    }
    
    obj@designs@dna <- getDesignMat(design=dnaDesign, annotations=obj@dnaAnnot)
    obj@designs@dna2rna <- getDesignMat(design=dnaDesign, 
                                        annotations=obj@rnaAnnot)
    obj@designs@rnaFull <- getDesignMat(design=rnaDesign, 
                                        annotations=obj@rnaAnnot)
    
    ## fit model
    obj@designs@rnaRed <- NULL
    obj@designs@rnaCtrlFull <- NULL
    obj@designs@rnaCtrlRed <- NULL
    obj@modelPreFits.dna.ctrl <- NULL
    obj@controls.forfit <- NULL
    
    message("Fitting model...")
    pb <- progress_bar$new(format = "[:bar] :percent (:current/:total)", 
                           total = NROW(obj@dnaCounts))
    models <- bplapply(rownames(obj@dnaCounts), function(rn) {
        tryCatch({
            pb$tick()
            return(fit.dnarna.noctrlobs(
                model=obj@model,
                dcounts=obj@dnaCounts[rn,,drop=FALSE],
                rcounts=obj@rnaCounts[rn,,drop=FALSE],
                ddepth=obj@dnaDepth,
                rdepth=obj@rnaDepth,
                rctrlscale=NULL,
                ddesign.mat=obj@designs@dna,
                rdesign.mat=obj@designs@rnaFull,
                d2rdesign.mat=obj@designs@dna2rna,
                rdesign.ctrl.mat=NULL,
                theta.d.ctrl.prefit=NULL,
                compute.hessian=FALSE))
        }, error = function(err) {message("error fitting: ", rn, ": ", err)})
    }, BPPARAM = obj@BPPARAM)
    names(models) <- rownames(obj@dnaCounts)
    obj@modelFits <- reformatModels(models)
    
    message("Analysis done!")
    return(obj)

}