

#' get a design matrix from the input design
#' 
#' adds an extra first column if this is the rna design matrix and
#' the rna noise model requires a separate variance parameter which is not
#' used for the mean model, ie. if the first parameter to estimate is not 
#' part of the GLM that defines the mean parameter.
#'
#' @import methods
#' 
#' @param design the formula object describing the design
#' @param annotations the column annotations corresponding to the design
#' @param condition condition to substract from formula
#'
#' @return a design matrix
getDesignMat <- function(design, annotations, condition=NULL) {
    if (is.null(design)) {
        design <- ~ 1
    } else if(!is.null(condition)){
        ## substract this condition from formula
        design <- update.formula(design, formula(paste0("~.-", condition)))
    }
    dmat <- model.matrix(design, annotations)
    return(dmat)
}

#' Return TRUE iff the reduced design is nested within the full design
#' @param full the full design (formula)
#' @param reduced the reduced design (formula)
#' @return TRUE iff the reduced design is nested in the full design
isNestedDesign <- function(full, reduced) {
    return(all(attr(terms(reduced), "term.labels") %in% 
                   attr(terms(full), "term.labels")))
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
    valid <- !vapply(models, is.null, TRUE)
    res <- list(
        ll = extractProp(models, "ll", valid),
        converged = extractProp(models, "converged", valid),
        
        d.coef = extractProp(models, "d.coef", valid),
        d.df = extractProp(models, "d.df", valid),
        d.se = extractProp(models, "d.se", valid),
        
        r.coef = extractProp(models, "r.coef", valid),
        r.df = extractProp(models, "r.df", valid),
        r.se = extractProp(models, "r.se", valid),
        
        r.ctrl.coef = extractProp(models, "r.ctrl.coef", valid),
        r.ctrl.df = extractProp(models, "r.ctrl.df", valid),
        r.ctrl.se = extractProp(models, "r.ctrl.se", valid)
    )
    
    return(res)
}

#' extract the given property from the list of models
#' @param models the models to extract the propety from
#' @param prop the name of the property
#' @param valids indices of the valid models in the input list
#' @return the formatted extracted property (NULL, vector or matrix)
extractProp <- function(models, prop, valids) {
    value <- models[valids][[1]][[prop]]
    if(is.null(value)) {
        res <- NULL
    } else if (length(value) == 1) {
        res <- rep(NA, length(models))
        res[valids] <- vapply(models[valids], function(x) x[[prop]], value)
        names(res) <- names(models)
    } else {
        res <- matrix(NA, nrow = length(models), ncol = length(value))
        res[valids,] <- do.call(rbind, lapply(models[valids], function(x) {
            if(is.null(x[[prop]])) {
                return(rep(NA, length(value)))
            } else {
                return(x[[prop]])
            }}))
        rownames(res) <- names(models)
    }
    return(res)
}

#' return the fitted value for alpha (transcription rate).
#' @param obj the MpraObject
#' @param term the term to get the alpha of (see details)
#' @param value the value of the term to get the alpha of (See details)
#' @param full if true, return alpha of the full model (default), otherwise of
#' the reduced model (only applies if an LRT-based analysis was used)
#' 
#' @return a named vector 
#' 
#' @details return the estimate for transcription rate as fitted by the package.
#' If the design is intercepted, then by default the baseline (intercept) rate 
#' is returned. Otherwise, term and value must be provided, such that term is a
#' valid term in the design provided to the fit, and value is one of the 
#' levels in the term.
getSingleAlpha <- function(obj, term=NULL, value=NULL, full=TRUE) {
    coefs <- getModelParameters_RNA(obj, full = full)
    if(full) {
        des <- obj@designs@rnaFull
    } else {
        des <- obj@designs@rnaRed
        if(is.null(des)) {
            stop("Reduced model can only be used for LRT based models")
        }
    }
    if(checkForIntercept(des) & is.null(term) & is.null(value)) {
        ##return the intercept
        coef <- exp(coefs[,2])
        names(coef) <- rownames(coefs)
        return(coef)
    } 
    if(is.null(term) | is.null(value)) {
        stop("both term and value must be provided")
    }
    coef.id <- colnames(obj@designs@rnaFull) %in% paste0(term, value)
    if(!any(coef.id)) {
        stop("no matching coefficient for given arguments")
    }
    coef.id <- 1 + which(coef.id)
    
    if(checkForIntercept(des)) {
        coef <- exp(coefs[,2] + coefs[,coef.id])
    } else {
        coef <- exp(coefs[,coef.id])
    }
    
    ##add the intercept from the control model
    if(!is.null(obj@designs@rnaCtrlFull)) {
        coef <- coef * exp(obj@modelPreFits.dna.ctrl$r.coef[1,2])
    }
    
    names(coef) <- rownames(coefs)
    return(coef)
}