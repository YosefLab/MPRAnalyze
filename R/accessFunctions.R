#' Get DNA model fits from an MpraObject
#' 
#' @param obj MpraObject to extract from
#' @param enhancers enhancer to extract 
#' @param depth include depth correction
#' @param full whether to extract from full model
#' @param transition use the DNA->RNA transition matrix (deafult: FALSE)
#' 
#' @return DNA fits (numeric, enhancers x samples)
#' 
#' @export
getFits.DNA <- function(obj, enhancers=NULL, depth=TRUE, full=TRUE, 
                        transition=FALSE){
    if(is.null(enhancers)) {
        enhancers <- names(obj@modelFits$ll)
    } else if(is.numeric(enhancers)) {
        enhancers <- names(obj@modelFits$ll)[enhancers]
    }
    if(full == TRUE){
        fit <- obj@modelFits
    } else {
        fit <- obj@modelFits.red
    }
    
    coef.mat <- t(fit$d.coef[enhancers,,drop=FALSE])
    coef.mat[is.na(coef.mat)] <- 0
    
    dmat <- obj@designs@dna
    if(transition) {
        dmat <- obj@designs@dna2rna
    }
    
    if(obj@model == "gamma.pois") {
        # alpha / rate
        dfit <- exp(coef.mat[1,] + dmat %*% coef.mat[-1,,drop=FALSE])
    } else if(obj@model == "ln.nb" | obj@model == "ln.ln") {
        dfit <- exp(dmat %*% coef.mat[-1,,drop=FALSE])
    }
    
    if(depth == TRUE){
        dfit <- dfit * replicate(NCOL(dfit), obj@dnaDepth)
    }
    
    colnames(dfit) <- enhancers
    if(transition) {
        rownames(dfit) <- rownames(obj@rnaAnnot)
    } else {
        rownames(dfit) <- rownames(obj@dnaAnnot)
    }
    return(t(dfit))
}

#' Get RNA full model fits from an MpraObject
#' 
#' @param obj MpraObject to extract from
#' @param enhancers enhancer to extract 
#' @param depth include depth correction
#' @param full whether to extract from full model
#' @param rnascale whether to use prefit rna model factors
#' 
#' @return RNA fits (numeric, enhancers x samples)
#' 
#' @export
getFits.RNA <- function(obj, enhancers=NULL, depth=TRUE, full=TRUE, 
                        rnascale=TRUE){
    if(is.null(enhancers)) {
        enhancers <- names(obj@modelFits$ll)
    }
    if(full == TRUE){
        fit <- obj@modelFits
        rdesign <- obj@designs@rnaFull
        rctrldesign <- obj@designs@rnaCtrlFull
    } else {
        fit <- obj@modelFits.red
        rdesign <- obj@designs@rnaRed
        rctrldesign <- obj@designs@rnaCtrlRed
    }
    if(!is.null(obj@rnaCtrlScale) & rnascale) {
        rctrlscale <- obj@rnaCtrlScale
    } else {
        rctrldesign <- NULL # overwrite initialisation
        rctrlscale <- NULL
    }
    
    dfit <- getFits.DNA(obj=obj, enhancers=enhancers, depth=FALSE, 
                        full=full, transition=TRUE)
    
    joint.des.mat <- cbind(rdesign, rctrldesign)
    
    coef.mat <- t(fit$r.coef[enhancers,-1,drop=FALSE])
    if(!is.null(obj@rnaCtrlScale) & rnascale) {
        coef.mat <- rbind(coef.mat,
                        replicate(NCOL(coef.mat), obj@rnaCtrlScale))
    }
    rfit <- exp(joint.des.mat %*% coef.mat)
    
    if(depth == TRUE){
        rfit <- rfit * replicate(NCOL(rfit), obj@rnaDepth)
    }
    
    rfit <- t(rfit) * dfit
    
    rownames(rfit) <- enhancers
    colnames(rfit) <- rownames(obj@rnaAnnot)
    return(rfit)
}

#' extract the DNA model parameters
#' 
#' @rdname extractModelParameters
#' @aliases extractModelParameters.DNA
#' extractModelParameters.RNA
#' 
#' @param obj the MpraObject to extract the parameters from
#' @param features the features to extract the parameters from (by default, 
#' parameters will be returned for all features)
#' @param full if TRUE (default), return the parameters of the full model. 
#' Otherwise, return the parameters of the reduced model (only relevant for 
#' LRT-based analyses)
#' @return a data.frame of features (rows) by parameters (cols). By convension, the
#' first parameter is related to the second moment, and the interpretation of 
#' it depends on the distributional model used (`alpha` for `gamma.pois`, variance 
#' for `ln.nb` and `ln.ln`)
#' @export
#' @rdname extractModelParameters
getModelParameters.DNA <- function(obj, features=NULL, full=TRUE) {
    if(is.null(features)) {
        features <- 1:NROW(obj@dnaCounts)
    } else if (is.character(features)) {
        features <- which(rownames(obj@dnaCounts) %in% features)
    }
    
    if(is.null(obj@modelFits)){
        stop("can't extract model parameters before fitting a model. An analysis function must be called first.")
    }
    if(full) {
        coef.mat <- obj@modelFits$d.coef[features,,drop=FALSE] 
    } else if (!full & !is.null(obj@modelFits.red)) {
        coef.mat <- obj@modelFits.red$d.coef[features,,drop=FALSE] 
    } else {
        stop("Parameters can't be extracted from reduced model, since analysis did not include fitting a reduced model")
    }
    colnames(coef.mat) <- c("disp", colnames(obj@designs@dna))
    rownames(coef.mat) <- rownames(obj@dnaCounts)[features]
    return(as.data.frame(coef.mat))
}

#' @export
#' @rdname extractModelParameters
getModelParameters.RNA <- function(obj, features=NULL, full=TRUE) {
    if(is.null(obj@modelFits)){
        stop("can't extract model parameters before fitting a model")
    }
    if(is.null(features)) {
        features <- 1:NROW(obj@dnaCounts)
    } else if (is.character(features)) {
        features <- which(rownames(obj@dnaCounts) %in% features)
    }
    
    if(full) {
        coef.mat <- obj@modelFits$r.coef[features,,drop=FALSE] 
        colnames(coef.mat) <- c("disp", colnames(obj@designs@rnaFull))
    } else if (!full & !is.null(obj@modelFits.red)) {
        coef.mat <- obj@modelFits.red$r.coef[features,,drop=FALSE] 
        colnames(coef.mat) <- c("disp", colnames(obj@designs@rnaRed))
    } else {
        stop("Reduced model unavailable")
    }
    
    rownames(coef.mat) <- rownames(obj@dnaCounts)[features]
    return(as.data.frame(coef.mat))
}

#' Get model distribution parameters from an MpraObject
#' 
#' @rdname getDistrParam
#' @aliases getDistrParam.DNA
#' getDistrParam.RNA
#' 
#' @param obj MpraObject to extract from
#' @param enhancer enhancer to extract 
#' @param full whether to extract from full model
#' 
#' @return fit parameters (numeric, samples x parameters)
#' 
#' @export
#' @rdname getDistrParam
getDistrParam.DNA <- function(obj, enhancer=NULL, full=TRUE){
    
    if(full == TRUE){
        fit <- obj@modelFits
    } else {
        fit <- obj@modelFits.red
    }
    
    coef.mat <- t(fit$d.coef[enhancer,,drop=FALSE])
    coef.mat[is.na(coef.mat)] <- 0
    
    if(obj@model == "gamma.pois") {
        par.shape <- as.vector(exp(coef.mat[1,]))
        par.rate <- as.vector(exp(-obj@designs@dna %*% coef.mat[-1,,drop=FALSE]))
        par <- data.frame(shape = par.shape, rate = par.rate)
    } else if(obj@model == "ln.nb") {
        par.sdlog <- as.vector( exp(coef.mat[1,]))
        par.meanlog <- as.vector(obj@designs@dna %*% coef.mat[-1,,drop=FALSE])
        par <- data.frame(meanlog = par.meanlog, sdlog = par.sdlog)
    } else if(obj@model == "ln.ln") {
        par.sdlog <- as.vector( exp(coef.mat[1,]))
        par.meanlog <- as.vector(obj@designs@dna %*% coef.mat[-1,,drop=FALSE])
        par <- data.frame(meanlog = par.meanlog, sdlog = par.sdlog)
    }
    
    rownames(par) <- rownames(obj@dnaAnnot)
    return(par)
}

#' @export
#' @rdname getDistrParam
getDistrParam.RNA <- function(obj, enhancer=NULL, full=TRUE){
    
    if(full == TRUE){
        fit <- obj@modelFits
    } else {
        fit <- obj@modelFits.red
    }
    rfit <- as.vector(getFits.RNA(obj, enhancers=enhancer, depth=FALSE, 
                                full=full, rnascale=TRUE))
    
    if(obj@model == "gamma.pois") {
        par.size <- as.vector(exp(fit$r.coef[enhancer,1]))
        par.mu <- rfit
        par <- data.frame(size = par.size, mu = par.mu)
    } else if(obj@model == "ln.nb") {
        par.size <- as.vector(exp(fit$r.coef[enhancer,1]))
        par.mu <- rfit
        par <- data.frame(size = par.size, mu = par.mu)
    } else if(obj@model == "ln.ln") {
        par.sdlog <- as.vector(exp(fit$r.coef[enhancer,1]))
        par.meanlog <- log(rfit)
        par <- data.frame(meanlog = par.meanlog, sdlog = par.sdlog)
    }
    
    rownames(par) <- rownames(obj@rnaAnnot)
    return(par)
}

#' return the fitted value for the transcription rate.
#' @param obj the MpraObject to extract from, must be after model fitting
#' @param by.factor return a matrix of values, corresponding to the estimated
#' rates of transcription under different values of a factor included in the 
#' design. Value must be of these options:
#' NULL: (default) return only the intercept term, a single baseline rate for 
#' each enhancer
#' "all": will return the corresponding transcription rates for all values 
#' included in the model
#' factor name: must be a factor included in the RNA annotations and the rna 
#' design. Will return the corresponding rates for all values of the given factor
#' @param full if true, return rate of the full model (default), otherwise of
#' the reduced model (only applies if an LRT-based analysis was used)
#' @export
#' @return the estimate for transcription rate as fitted by the model
getAlpha <- function(obj, by.factor=NULL, full=TRUE) {
    des <- obj@designs@rnaFull
    if(!full) {
        des <- obj@designs@rnaRed
        if(is.null(des)) {
            stop("reduced model not available")
        }
    }
    
    if(is.null(by.factor)) {
        return(data.frame(TR=getSingleAlpha(obj, full=full)))
    } else if (by.factor=="all") {
        ## extract all parametres except for the dispersion
        alpha.mat <- getModelParameters.RNA(obj, full=full)[,-1]
        if(checkForIntercept(des)) {
            ## add the intercept to all other columns
            alpha.mat[,-1] <- alpha.mat[,-1] + alpha.mat[,1]
        }
        return(as.data.frame(exp(alpha.mat)))
    } else {
        l <- levels(obj@rnaAnnot[,by.factor])
        ## if intercepted: first level is different
        if(checkForIntercept(des)) {
            int.alpha <- getSingleAlpha(obj)
            l <- l[-1]
        } else {
            int.alpha <- NULL
        }
        
        alpha.mat <- do.call(cbind, lapply(l, function(v) {
            getSingleAlpha(obj, term = by.factor, value = v, full = full)
        }))
        alpha.mat <- cbind(int.alpha, alpha.mat)
        colnames(alpha.mat) <- levels(obj@rnaAnnot[,by.factor])
        
        return(as.data.frame(alpha.mat))
    }
}

#' extract library depth factors.
#' @param obj the object o extract the factors from
#' @return a list with two elements: 'dna' with the dna library depth factors, 
#' and 'rna' with the rna factors
#' @export
getLibraryDepthFactors <- function(obj) {
    return(list(dna=obj@dnaDepth, rna=obj@rnaDepth))
}
