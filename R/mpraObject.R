setClass("MpraObject", slots = c(
    ## provided by user
    dnaCounts = "matrix",
    rnaCounts = "matrix",
    colAnnot = "data.frame",
    dnaDesign = "Matrix", ##sparse
    rnaDesign = "Matrix", ##sparse
    controls = "integer", #idx of negative controls

    dneDepth = "numeric",
    rnaDepth = "numeric",
    model = "character",
    condition = "character",

    modelFits = "list"
))

#' an initialize function for the MpraObject class
#'
#' @param dnaCounts the DNA counts matrix
#' @param rnaCounts the RNA counts matrix
#' @param colAnnot column annotations
#' @param controls a vector specifying which enhancers are negative controls
#' (scrambles)
#' @param dnaDesign a formula or design matrix describing the design of the DNA
#' counts
#' @param rnaDesign a formula or design matrix describing the design of the DNA
#' counts
#' @model The model to use (optional), options are: TODO
#' @condition if differential activity analysis is desired, this parameter determines
#' which column in colAnnot determines the condition to test across
#'
#' @examples
#' TODO: generate data in two conditions, create dataframe of the conditions,
#' then create an object. Also have one of a single condition.
MpraObject <- function(dnaCounts, rnaCounts, colAnnot=NULL, controls=NULL,
                       dnaDesign=NULL, rnaDesign=NULL, model=NULL,
                       condition=NULL) {

    obj <- new("MpraObject", dnaCounts=dnaCounts, rnaCounts=rnaCounts,
               colAnnot=colAnnot, controls=controls)

    if(!is.null(condition)) {
        if(!(condition %in% colnames(colAnnot))) {
            stop("Condition must be a valid column name in colAnnot")
        }
        obj@condition <- condition
    }

    ## even if a condition is specified, it doesn't need to be the foremost
    ## term in the DNA formula, since the testing is done on the RNA model
    obj@dnaDesign <- getDesign(dnaDesign, colAnnot, condition=NULL)
    obj@rnaDesign <- getDesign(rnaDesign, colAnnot, condition)

    if(!is.null(model)) {
        ##TODO: make sure this is a valid model ID
        obj@model <- model
    }

    return(obj)
}

#' get a design Matrix from the user input
#'
#' TODO: decide if the matrix should be sparse or not. Few copies, mostly
#' multiplied, dense might be better despite sparsity.
#'
#' @param design the design input
#' @param colAnnot the column annotations
#' @param condition the condition of interest, if provided. Must be foremost
#' term in an input of type formula.
#'
#' @return the design matrix
getDesign <- function(design, colAnnot, condition) {
    if(is.null(design)) {
        design <- Matrix(rep(1, NROW(colAnnot)), ncol = 1)
    } else if(is(design, "Matrix")) {
        return(design)
    } else if (is.matrix(design)) {
        return(Matrix(design, sparse=TRUE))
    } else if (is(design, "formula")) {
        ##TODO: make sure condition is first if exists, and create model matrix
        ##TODO: if formula has response, remove it
        return(sparse.model.matrix(design, colAnnot))
    } else {
        stop("Invalid design input")
    }
}
