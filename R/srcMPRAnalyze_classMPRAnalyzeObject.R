################################################################################
#################     MPRAnalyze output container class     ####################
################################################################################

### 1. Define output container class

# Define class unions for slots
# DO NOT USE SPACE BETWEEN members=c(..)
# THAT CAUSES ERROR DURING R CMD BUILD

#' @importFrom methods setClassUnion
setClassUnion('numericORNULL', members=c('numeric', 'NULL'))
setClassUnion('characterORNULL', members=c('character', 'NULL'))
setClassUnion('listORNULL', members=c('list', 'NULL'))
setClassUnion('data.frameORNULL', members=c('data.frame', 'NULL'))

#' Container class for MPRAnalyze output
#' 
#' MPRAnalyze output and intermediate results including as model fits.
#' @name MPRAnalyzeObject-class
#'      
#' @author David Sebastian Fischer
setClass(
  'MPRAnalyzeObject',
  slots = c(
    dfMPRAnalyzeResults = "data.frameORNULL",
    lsModelFitsFull      = "listORNULL",
    lsModelFitsRed       = "listORNULL",
    lsDNAModelFitsCtrl   = "listORNULL",
    lsDNAModelFits       = "listORNULL",
    vecDispersions       = "numericORNULL",
    dfDispersions        = "data.frameORNULL",
    matDNACountsProc     = "matrix",
    matRNACountsProc     = "matrix",
    vecAllIDs            = "characterORNULL",
    vecCtrlIDs           = "characterORNULL",
    dfAnnotationProc     = "data.frame",
    vecDNADepth   = "numericORNULL",
    vecRNADepth   = "numericORNULL",
    vecModelFacRNAFull  = "characterORNULL",
    vecModelFacRNAFullCtrl   = "characterORNULL",
    vecModelFacRNARed  = "characterORNULL",
    vecModelFacRNARedCtrl  = "characterORNULL",
    vecModelFacDNA   = "characterORNULL",
    strModel             = "characterORNULL",
    scaNProc             = "numericORNULL",
    strReport            = "characterORNULL"
  )
)

### 1. Enable accession of main output via list-like
### properties of MPRAnalyzeObject

#' List-like accessor methods for MPRAnalyze
#' 
#' Allow usage of MPRAnalyze ouput object like a list with
#' respect to the main output: MPRAnalyze
#' 
#' @param x (MPRAnalyzeObject) MPRAnalyze output object.
#' @param i,name (idx or str) Name or index of core output element of MPRAnalyzeObject
#' @param j       Not used, only vectors.
#' @param ...     Not used.
#' 
#' @name list_accession
#' @aliases names.MPRAnalyzeObject
#' names,MPRAnalyzeObject-method          
#' $,MPRAnalyzeObject-method          
#' [[,MPRAnalyzeObject,character,missing-method
#' 
#' @author David Sebastian Fischer
NULL

#' @return Names of core output in MPRAnalyzeObject.
#' @name list_accession
#' @export
setMethod('names', 'MPRAnalyzeObject', function(x) {
  return( c("dfMPRAnalyzeResults") )
})

#' @return Target element from MPRAnalyzeObject
#' @name list_accession
#' @export
setMethod('[[', c('MPRAnalyzeObject', 'character', 'missing'), function(x, i, j, ...){
  if(identical(i, "dfMPRAnalyzeResults")){ 
      return(x@dfMPRAnalyzeResults)
  } else { 
      return(NULL) 
  }
})

#' @return Target element from MPRAnalyzeObject
#' @name list_accession
#' @export
setMethod('$', 'MPRAnalyzeObject', function(x, name) x[[name]] )
