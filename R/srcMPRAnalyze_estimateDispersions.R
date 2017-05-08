#' Estimate dispersion parameters on MPRA data set
#' 
#' Uses DESeq2 model and package components to estimate
#' a MAP (trend-smoothed) dispersion estimator for 
#' the RNA distribution per gene. The intial 
#' MLE is fitted by this package according to MPRA 
#' characteristics.
#' 
#' @return (data frame genes x characteristics)
#' Table containing details of dispersion fit per gene.
#' 
#' @author David Sebastian Fischer
#' 
#' @export
estimateDispersion <- function(
  objMPRAnalyze,
  vecModelFactors,
  boolVerbose=TRUE){
  
  # Create DESeq2 object from 
  objectDESeq2 <- DESeqDataSetFromMatrix(
    countData=objMPRAnalyze@matRNACountsProc, 
    colData=objMPRAnalyze@dfAnnotationProc, 
    design=as.formula(paste0(
      "~", paste0(vecModelFactors, collapse="+"))) )
  # Estimate size factors on DESeq2 object
  objectDESeq2 <- estimateSizeFactors(objectDESeq2)
  # Step 1 DESeq2 dispersion estimation: MLE
  # This step is replaced by custom MLE
  lsFits_Disp <- fitModels(
    objectMPRAnalyze=objMPRAnalyze,
    vecModelFactors=vecModelFactors,
    MAXIT=1000)
  # Substitute custom dispersion MLE into DESeq2 object
  vecDispMLE <- 1/sapply(lsFits_Disp, function(f) f$scaDisp )
  objectDESeq2 <- estimateDispersionsGeneEst(objectDESeq2) # to get correct structure, check this
  mcols(objectDESeq2)$dispGeneEst <- vecDispMLE # replace
  # Step 2 DESeq2 dispersion estimation: Estimate trend
  objectDESeq2 <- estimateDispersionsFit(
    object=objectDESeq2,
    fitType=c("parametric"),
    quiet=!boolVerbose)
  # Step 3 DESeq2 dispersion estimation: MAP
  objectDESeq2 <- estimateDispersionsMAP(
    object=objectDESeq2, 
    quiet=!boolVerbose)
  
  dfDispEst <- as.data.frame(mcols(objectDESeq2))
  vecDispersions <- 1/mcols(objectDESeq2)$dispersion
  return(list(dfDispersions=dfDispEst,
              vecDispersions=vecDispersions))
}