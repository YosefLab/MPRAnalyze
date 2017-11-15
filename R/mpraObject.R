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
    if (any(dim(dnaCounts) != dim(rnaCounts))) {
        errors <- c(errors,
                    "DNA, RNA matrix must be of same dimensions")
    }
    if (is.null(rownames(dnaCounts)) |
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

    hyptestResults = "data.frame",

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

    if(is(controls, "character")) {
        controls <- which(rownames(dnaCounts) %in% controls)
    } else if(is(controls, "logical")) {
        controls <- which(controls)
    }

    obj <- new("MpraObject", dnaCounts=dnaCounts, rnaCounts=rnaCounts,
               colAnnot=colAnnot, controls=controls, BPPARAM=BPPARAM)
    return(obj)
}
