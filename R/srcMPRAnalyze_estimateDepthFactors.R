#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++     Compute library depth factors    ++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Compute a size factor for each sample
#' 
#' @seealso Called by \link{computeNormConst}.
#' 
#' @param matCountData (matrix genes x samples)
#'    Read count data.
#' 
#' @return vecSizeFactors (numeric vector number of samples) 
#'    Model scaling factors for each sample which take
#'    sequencing depth into account (size factors).
#' 
#' @author David Sebastian Fischer
computeSizeFactors <- function(matCountData){
  
  # Compute geometric count mean over replicates
  # for genes without zero observations: Samples
  # with more than half zero observations receive 
  # size factor =1 otherwise.
  vecboolZeroObs <- apply(matCountData,1,function(gene){
      !any(!is.na(gene) & gene==0)
  })
  # Take geometric mean
  vecGeomMean <- apply(matCountData[vecboolZeroObs,], 1, function(gene){
    ( prod(gene[!is.na(gene)]) )^( 1/sum(!is.na(gene)) )
  })
  
  # Chose median of ratios over genes as size factor
  vecSizeFactors <- apply(matCountData[vecboolZeroObs,], 2, function(sample){
    median(sample/vecGeomMean, na.rm=TRUE)
  })
  
  if(any(vecSizeFactors==0)){
    print("WARNING: Found size factors==0, setting these to 1.")
    vecSizeFactors[vecSizeFactors==0] <- 1
  }
  names(vecSizeFactors) <- colnames(matCountData)
  
  return(vecSizeFactors)
}

#' Set depth factors in objectMPRAnalyze object
#' 
#' @author David Sebastian Fischer
#' 
#' @export
estimateDepthFactors <- function(objectMPRAnalyze){
  
  # Compute size factors if not supplied by user.
  # a) RNA
  if(is.null(objectMPRAnalyze@vecRNADepth)){
    objectMPRAnalyze@vecRNADepth <- computeSizeFactors(matCountData=objectMPRAnalyze@matRNACountsProc)
  }
  # b) DNA
  if(is.null(objectMPRAnalyze@vecDNADepth)){
    objectMPRAnalyze@vecDNADepth <- computeSizeFactors(matCountData=objectMPRAnalyze@matDNACountsProc)
  }
  
  return(objectMPRAnalyze)
}