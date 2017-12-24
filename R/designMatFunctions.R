

#' get a design matrix from the input design
#' 
#' adds an extra first column if this is the rna design matrix and
#' the rna noise model requires a separate variance parameter which is not
#' used for the mean model, ie. if the first parameter to estimate is not 
#' part of the GLM that defines the mean parameter.
#'
#' @param obj the MpraObject
#' @param design the input design. A matrix is returned as is, a formula is
#' expanded to a model matrix using the object colAnnot. If design is NULL, an
#' intercept-only design matrix is created and returned
#' @param testcondition condition to substract from formula
#' @param rna whether this is a rna model. Together with the noise model set
#' in obj@model, this defines whether an extra column for the variance link
#' parameter is allowed in the design matrix.
#'
#' @return a design matrix
getDesignMat <- function(obj, design, testcondition=NULL) {
    # deprecate matrix entry?
    if (is.matrix(design)) {
        dmat <- design
    } else if (is.null(design)) {
        dmat <- matrix(rep(1,NCOL(obj@dnaCounts)), ncol = 1,
                       dimnames = list(colnames(obj@dnaCounts), "(intercept)"))
    } else if (is(design, "formula")) {
        if(!is.null(testcondition)){
            ## substract this condition from formula
            terms <- attr(terms.formula(design), "term.labels")
            termsnew <- terms[terms != testcondition]
            if(length(termsnew) >= 1) {
                design <- as.formula(paste0("~", paste(termsnew, collapse="+")))
            } else {
                design <- ~1
            }
        }
        dmat <- model.matrix(design, obj@colAnnot)
    } else {
        stop("invalid design")
    }
    return(dmat)
}

#' TODO
#' Return TRUE iff the given design has an intercept term
#' @param design either a formula or a design matrix
#' @return TRUE iff the design has an intercept term
checkForIntercept <- function(design) {
    if(is.matrix(design)) {
        return(all(design[,1] == 1))
    } else {
        return(as.logical(attr(terms(design), "intercept")))
    }
}