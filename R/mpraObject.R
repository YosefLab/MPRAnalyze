setClassUnion('Design', members = c('matrix', 'formula'))

setClass("Designs", slots = c(
    dnaFull = "Design",
    rnaFull = "Design",

    ## only used in differential LRT mode
    dnaRed = "Design",
    rnaRed = "Design"
))

#' validate an MpraObject
#'
#' @param object the MpraObject instance to validate
#'
#' @return TRUE if object is valid, otherwise vector of character strings
#' explaining how it's invalid
validateMpraObject <- function(object) {
    errors = character()
    if (dim(dnaCounts) != dim(rnaCounts)) {
        errors <- c(errors,
                    "DNA, RNA matrix must be of same dimensions")
    }
    if (is.null(rownnames(dnaCounts)) |
        any(rownames(dnaCounts) != rownames(rnaCounts))) {
        errors <- c(errors,
                    "RNA, DNA feature names either missing or don't match")
    }

    if(length(errors) > 0) {
        return(errors)
    } else {
        return(TRUE)
    }
}

##TODO: document this
setClass("MpraObject", validity = validateMpraObject,
         slots = c(
    ## provided by user
    dnaCounts = "matrix",
    rnaCounts = "matrix",
    colAnnot = "data.frame",
    controls = "integer", #idx of negative controls

    dnaDepth = "numeric",
    rnaDepth = "numeric",

    mode = "character",
    model = "character",
    designs = "Designs",
    modelFits = "list",
    modelFits.reduced = "list", ##only used for LRT diff mode

    ##TODO: analysis results containers?
    hyptestResults = "data.frame"

    BPPARAM = "BiocParallelParam"
))


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

    if(!is.na(controls)) {
        if(is(controls, "character")) {
            ctrl <- which(rownames(dnaCounts) %in% controls)
        } else if(is(controls, "logical")) {
            ctrl <- which(controls)
        } else if(is(controls, "integer")) {
            ctrl <- controls
        } else {
            stop("controls must be integer, character or logical vector")
        }
    }

    obj <- new("MpraObject", dnaCounts=dnaCounts, rnaCounts=rnaCounts,
               colAnnot=colAnnot, controls=ctrl, BPPARAM=BPPARAM)
    return(obj)
}
