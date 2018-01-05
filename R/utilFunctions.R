

#' get a design matrix from the input design
#' 
#' adds an extra first column if this is the rna design matrix and
#' the rna noise model requires a separate variance parameter which is not
#' used for the mean model, ie. if the first parameter to estimate is not 
#' part of the GLM that defines the mean parameter.
#'
#' @import methods
#'
#' @param obj the MpraObject
#' @param design the input design. A matrix is returned as is, a formula is
#' expanded to a model matrix using the object colAnnot. If design is NULL, an
#' intercept-only design matrix is created and returned
#' @param testcondition condition to substract from formula
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

#' Reformat a list of models (return values of a fit.* function) to a list
#' of model parameters containing the corresponding parameters, for easier 
#' access and control
#' @param models the models to reformet, should be a list of results from a 
#' fit.* function
#' @return a list of formatted extacted properties 
reformatModels <- function(models) {
    res = list(ll = extractProp(models, "ll"),
               converged = extractProp(models, "converged"),
               
               d.coef = extractProp(models, "d.coef"),
               d.df = extractProp(models, "d.df"),
               d.se = extractProp(models, "d.se"),
               
               r.coef = extractProp(models, "r.coef"),
               r.df = extractProp(models, "r.df"),
               r.se = extractProp(models, "r.se"),
               
               r.ctrl.coef = extractProp(models, "r.ctrl.coef"),
               r.ctrl.df = extractProp(models, "r.ctrl.df"),
               r.ctrl.se = extractProp(models, "r.ctrl.se")
               )
    
    return(res)
}

#' extract the given property from the list of models
#' @param models the models to extract the propety from
#' @param prop the name of the property
#' @return the formatted extracted property (NULL, vector or matrix)
extractProp <- function(models, prop) {
    value <- models[[1]][[prop]]
    if(is.null(value)) {
        res <- NULL
    } else if (length(value) == 1) {
        res <- vapply(models, function(x) x[[prop]], value)
        names(res) <- names(models)
    } else {
        res <- do.call(rbind, lapply(models, function(x) {
            if(is.null(x[[prop]])) {
                return(rep(NA, length(value)))
            } else {
                return(x[[prop]])
            }}))
        rownames(res) <- names(models)
    }
    return(res)
}
