#' @import Biobase
#' @import BiocParallel
#' @importFrom compiler cmpfun
#' @import ggplot2
#' @importFrom grDevices dev.off graphics.off pdf
#' @import knitr
#' @import Matrix
#' @import methods
#' @importFrom stats dnbinom median optim p.adjust pchisq rnbinom rnorm runif sd time
#' @import SummarizedExperiment
#' @importFrom utils packageDescription
NULL

#library(compiler)
#library(BiocParallel)
#library(Biobase)
#library(methods)

#wdCWD <- getwd()
#setwd("~/gitDevelopment/MPRAnalyze")
#setwd("/data/yosef2/users/fischerd/software/MPRAnalyze")
#source("R/srcMPRAnalyze_classMPRAnalyzeObject.R")
#source("R/srcMPRAnalyze_CostFunctions.R")
#source("R/srcMPRAnalyze_estimateDepthFactors.R")
#source("R/srcMPRAnalyze_estimateDispersions.R")
#source("R/srcMPRAnalyze_extractModel.R")
#source("R/srcMPRAnalyze_fitModels.R")
#source("R/srcMPRAnalyze_getModelFits.R")
#source("R/srcMPRAnalyze_plotBoxplots.R")
#source("R/srcMPRAnalyze_prefitCtrlModels.R")
#source("R/srcMPRAnalyze_doHypothesisTests.R")
#setwd(wdCWD)

#' Wrapper function MPRAnalyze
#' 
#' Performs model fitting and differential expression analysis
#' steps required within MPRAnalyze framework.
#' 
#' Note that if UMI barcodes are used to label plasmids 
#' with the same enhancer and these barcodes differ between enhancers, 
#' the UMI count of the first barcode of each plasmid 
#' is treated as the first sample etc.
#' 
#' @param matRNACounts (integer count matrix enhancers x samples)
#' RNA counts per enhancer and sample. 
#' @param matDNACounts (integer count matrix enhancers x samples)
#' DNA counts per enhancer and sample.
#' @param dfAnnotation (data frame samples x meta data)
#' Contains sample names and associated batch labels. All model
#' factors listed below which are used to group samples into batches
#' must be column names in dfAnnotation.
#' @param vecCtrlIDs (vector of strings number of control enhancers)
#' [default NULL] IDs of enhancers which are scrambled/controls and
#' used to define a background set.
#' @param strModelFull (string) Not yet supported. 
#' formula() based model input.
#' Alternative to vecModelFacRNAFull and vecModelFacRNAFullCtrl.
#' @param strModelRed (string) Not yet supported.
#' formula() based model input.
#' Alternative to vecModelFacRNARed and vecModelFacRNARedCtrl.
#' @param vecModelFacRNAFull (vector of strings lenght number of factors
#' in full expression model)
#' [Default NULL]
#' The model factors (columns in dfAnnotation) which are to be used 
#' in the full expression model. The expression model will have batch
#' adjustment for these factors.
#' @param vecModelFacRNAFullCtrl (vector of strings lenght number of 
#' additional factors in full expression model for control enhancers)
#' [Default NULL]
#' The model factors (columns in dfAnnotation) which are to be used 
#' in the full control expression model. The expression model will have batch
#' adjustment for these factors.
#' @param vecModelFacRNARed (vector of strings lenght number of factors
#' in reduced expression model)
#' [Default NULL]
#' The model factors (columns in dfAnnotation) which are to be used 
#' in the reduced expression model. The expression model will have batch
#' adjustment for these factors.
#' @param vecModelFacRNARedCtrl (vector of strings lenght number of 
#' additional factors in reduced expression model for control enhancers)
#' [Default NULL]
#' The model factors (columns in dfAnnotation) which are to be used 
#' in the reduced control expression model. The expression model will have batch
#' adjustment for these factors.
#' @param vecModelFacDNA (vector of strings lenght number of factors
#' in DNA model)
#' [Default NULL]
#' The model factors (columns in dfAnnotation) which are to be used 
#' in the DNA model. The DNA model will have batch adjustment for these factors.
#' @param scaNProc (scalar) [Default 1]
#' Number of processes to be used during parallelisation.
#' @param vecDispersionsExternal Not yet supported. 
#' External dispersion factors.
#' @param vecRNADepthExternal (numeric vector length number of samples)
#' Library depth correction factors for RNA data set, one per sample.
#' @param vecDNADepthExternal (numeric vector length number of samples)
#' Library depth correction factors for DNA data set, one per sample.
#' @param boolPreFitCtrlDNA (bool) [Default False]
#' Whether to pre-fit the control DNA models.
#' Can only be used if strModel=="gammaDNApoisRNA".
#' The size of the estimation problem is of gammaDNApoisRNA
#' is very large if many control sequences are given. Setting
#' this pre-fitting to TRUE reduces the problem size by conditioning
#' full and reduced model on pre-estimated DNA models of the control sequences.
#' This likely comes at little cost in accuracy.
#' @param strModel (string) [Default "lnDNAnbRNA"]
#' {"lnDNAnbRNA", "gammaDNApoisRNA", "gammaDNApoisRNA_coordascent"}
#' Model and estimation scheme to be used:
#' lnDNAnbRNA: DNA point estimators are pre-estimated based on a lognormal model.
#' The expression model is then fit based on a negative binomial likelihood of the 
#' RNA counts.
#' gammaDNApoisRNA: DNA and expression model are co-estimated, exloiting a closed
#' form marginalisation. The DNA counts are modeled with a beta distribution
#' and the RNA counts with a poisson distribution.
#' gammaDNApoisRNA_coordascent: As gammaDNApoisRNA but can be used to estimate
#' models with many control enhancers. The problem is broken up into iterative
#' DNA and expression model estimation until it converges to a local optimum of
#' the likelihood function.
#' @param boolVerbose (bool) [Default TRUE]
#' Whether to print progress statements.
#' 
#' @return (MPRAnalyzeObject)
#' Output container class object with data frame 
#' with results available via list-like properties.
#' 
#' @author David Sebastian Fischer
#' 
#' @export
runMPRAnalyze <- function(
    matRNACounts, 
    matDNACounts,
    dfAnnotation,
    vecCtrlIDs=NULL,
    strModelFull=NULL,
    strModelRed=NULL,
    vecModelFacRNAFull=NULL,
    vecModelFacRNAFullCtrl=NULL,
    vecModelFacRNARed=NULL,
    vecModelFacRNARedCtrl=NULL,
    vecModelFacDNA=NULL,
    scaNProc=1, 
    vecDispersionsExternal=NULL,
    vecRNADepthExternal=NULL,
    vecDNADepthExternal=NULL,
    boolPreFitCtrlDNA=TRUE,
    strModel="gammaDNApoisRNA",
    boolVerbose=TRUE ){
    
    strMessage <- "MPRAnalyze v0.90 for MPRA data (DNA and RNA counts)"
    if(boolVerbose) message(strMessage)
    strReport <- strMessage
    
    # 1. Process input
    vecAllIDs <- setdiff(rownames(matRNACounts), vecCtrlIDs)
    lsModelFac <- extractModel(
        strModelFull=strModelFull,
        strModelRed=strModelRed,
        vecModelFacRNAFull=vecModelFacRNAFull,
        vecModelFacRNAFullCtrl=vecModelFacRNAFullCtrl,
        vecModelFacRNARed=vecModelFacRNARed,
        vecModelFacRNARedCtrl=vecModelFacRNARedCtrl)
    obj <- new(
        'MPRAnalyzeObject',
        dfMPRAnalyzeResults  = NULL,
        lsModelFitsFull      = NULL,
        lsModelFitsRed       = NULL,
        lsDNAModelFits       = NULL,
        lsDNAModelFitsCtrl   = NULL,
        vecDispersions       = vecDispersionsExternal,
        dfDispersions        = NULL,
        matRNACountsProc     = matRNACounts,
        matDNACountsProc     = matDNACounts,
        vecAllIDs            = vecAllIDs,
        vecCtrlIDs           = vecCtrlIDs,
        dfAnnotationProc     = dfAnnotation,
        vecRNADepth          = vecRNADepthExternal,
        vecDNADepth          = vecDNADepthExternal,
        vecModelFacRNAFull   = lsModelFac$vecModelFacRNAFull,
        vecModelFacRNAFullCtrl = lsModelFac$vecModelFacRNAFullCtrl,
        vecModelFacRNARed    = lsModelFac$vecModelFacRNARed,
        vecModelFacRNARedCtrl= lsModelFac$vecModelFacRNARedCtrl,
        vecModelFacDNA       = vecModelFacDNA,
        strModel             = strModel,
        scaNProc             = scaNProc,
        strReport            = strReport)
    
    # Set parallelisation
    if(scaNProc > 1){ register(MulticoreParam(workers=scaNProc)) 
    } else { register(SerialParam()) }
    
    strMessage <- "# Estimate library depth factors"
    if(boolVerbose) message(strMessage)
    obj@strReport <- paste0(obj@strReport, "\n", strMessage)
    obj <- estimateDepthFactors(obj=obj)
    
    # 2. Estimate dispersion parameters if not supplied
    # NOT USED
    if(is.null(obj@vecDispersions) &
       obj@strModel == "pointlnDNAnbRNA" & FALSE){
        strMessage <- "# Estimate over-dispersion parameters"
        if(boolVerbose) message(strMessage)
        obj@strReport <- paste0(obj@strReport, "\n", strMessage)
        lsDispEst <- estimateDispersion(
            objMPRAnalyze=obj,
            vecModelFac=obj@vecModelFacFull,
            boolVerbose=boolVerbose)
        obj@vecDispersions <- lsDispEst$vecDispersions
        obj@dfDispersions <- lsDispEst$dfDispersions
    }
    
    # 3. Fit full and reduced model
    # Pre-fit DNA models in marginalisation scenario
    # to reduce size of optimisation problems: gammaDNApoisRNA
    if(boolPreFitCtrlDNA & obj@strModel=="gammaDNApoisRNA"){
        strMessage <- "# Pre-Fit control DNA model"
        if(boolVerbose) message(strMessage)
        obj@strReport <- paste0(obj@strReport, "\n", strMessage)
        obj@lsDNAModelFitsCtrl <- prefitCtrlDNAModels(
            obj=obj, vecModelFacRNA=obj@vecModelFacRNAFull,
            vecModelFacDNA=obj@vecModelFacDNA,
            MAXIT=1000, boolVerbose=boolVerbose)
    }
    # DNA models are shared between null and full if expression
    # model is conditioned on DNA model: pointlnDNAnbRNA
    if(obj@strModel=="pointlnDNAnbRNA"){
        strMessage <- "# Fit DNA models"
        if(boolVerbose) message(strMessage)
        obj@strReport <- paste0(obj@strReport, "\n", strMessage)
        obj@lsDNAModelFits <- fitModels(
            obj=obj, vecModelFacRNA=NULL,
            vecModelFacRNACtrl=NULL,
            vecModelFacDNA=obj@vecModelFacDNA,
            boolFitDNA = TRUE,
            MAXIT=1000, boolVerbose=boolVerbose)
    }
    
    strMessage <- "# Fit full model"
    if(boolVerbose) message(strMessage)
    obj@strReport <- paste0(obj@strReport, "\n", strMessage)
    obj@lsModelFitsFull <- fitModels(
        obj=obj, vecModelFacRNA=obj@vecModelFacRNAFull,
        vecModelFacRNACtrl=obj@vecModelFacRNAFullCtrl,
        vecModelFacDNA=obj@vecModelFacDNA,
        MAXIT=1000, boolVerbose=boolVerbose)
    
    if(!is.null(vecModelFacRNARed)) {
        strMessage <- "# Fit null model"
        if(boolVerbose) message(strMessage)
        obj@strReport <- paste0(obj@strReport, "\n", strMessage)
        obj@lsModelFitsRed <- fitModels(
            obj=obj, vecModelFacRNA=obj@vecModelFacRNARed,
            vecModelFacRNACtrl=obj@vecModelFacRNARedCtrl,
            vecModelFacDNA=obj@vecModelFacDNA,
            MAXIT=1000, boolVerbose=boolVerbose)
        
        # 4. Perform differential expression analysis
        strMessage <- "# Perform differential expression analysis"
        if(boolVerbose) message(strMessage)
        obj@strReport <- paste0(obj@strReport, "\n", strMessage)
        obj <- doHypothesisTests(obj=obj)
    }
    
    return(obj)
}