setClassUnion('Design', members = c('matrix', 'formula'))

setClass("Designs", slots = c(
    dnaFull = "Design",
    rnaFull = "Design",

    ## only used in differential LRT mode
    dnaRed = "Design",
    rnaRed = "Design"
))

setClass("MpraObject", slots = c(
    ## provided by user
    dnaCounts = "matrix",
    rnaCounts = "matrix",
    colAnnot = "data.frame",
    controls = "integer", #idx of negative controls

    dnaDepth = "numeric",
    rnaDepth = "numeric",

    model = "character",
    designs = "Designs",
    modelFits = "list",
    reducedModelFits = "list", ##only used for LRT diff mode

    ##TODO: analysis results containers?

    BPPARAM = "BiocParallelParam"
))

##TODO: validity stuff, check rownames\colnames match in data, dimensions match, etc

#' Initialize a MpraObject object
#'
#' @import BiocParallel
#'
#' @param dnaCounts the DNA counts matrix
#' @param rnaCounts the RNA counts matrix
#' @param colAnnot column annotations
#' @param controls a vector specifying which enhancers are negative controls
#' (scrambles)
#' @param BPPARAM the biocParalell backend to use for parallelization throughout
#' the analysis
#'
#' @export
#'
#' @examples
#' ##TODO
MpraObject <- function(dnaCounts, rnaCounts, colAnnot=NULL, controls=NA_integer_,
                       BPPARAM=NULL) {
    if(is.null(BPPARAM)) {
        BPPARAM <- bpparam()
    }
    obj <- new("MpraObject", dnaCounts=dnaCounts, rnaCounts=rnaCounts,
               colAnnot=colAnnot, controls=controls, BPPARAM=BPPARAM)
    return(obj)
}
