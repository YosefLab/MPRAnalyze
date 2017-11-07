setClass("MpraObject", slots = c(
    ## provided by user
    dnaCounts = "matrix",
    rnaCounts = "matrix",
    colAnnot = "data.frame",
    controls = "integer", #idx of negative controls

    dneDepth = "numeric",
    rnaDepth = "numeric",

    mode = "character", ## quantitative, diff-one.vs.all, diff-all.vs.all
    model = "character",
    condition = "character",
    designFormulas = "list",
    designMats = "list",
    modelFits = "list"
    ##TODO: analysis results containers?
))

#' Initialize a MpraObject object

#' an initialize function for the MpraObject class
#'
#' @param dnaCounts the DNA counts matrix
#' @param rnaCounts the RNA counts matrix
#' @param colAnnot column annotations
#' @param controls a vector specifying which enhancers are negative controls
#' (scrambles)
#' @export
#'
#' @examples
#' ##TODO
MpraObject <- function(dnaCounts, rnaCounts, colAnnot=NULL, controls=NULL) {

    obj <- new("MpraObject", dnaCounts=dnaCounts, rnaCounts=rnaCounts,
               colAnnot=colAnnot, controls=controls)
    return(obj)
}
