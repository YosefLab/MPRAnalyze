
setClassUnion('listORNULL', members=c('list', 'NULL'))
setClassUnion('numericORNULL', members=c('numeric', 'NULL'))
setClassUnion('integerORNULL', members=c('integer', 'NULL'))
setClassUnion('Design', members = c('matrix', 'formula', 'NULL'))

setClass("Designs", slots = c(
    dna = "Design",
    
    rnaFull = "Design",
    rnaCtrlFull = "Design",
    
    ## only used in differential LRT mode
    rnaRed = "Design",
    rnaCtrlRed = "Design"
))

#' validate an MpraObject
#'
#' @param object the MpraObject instance to validate
#'
#' @return TRUE if object is valid, otherwise vector of character strings
#' explaining how it's invalid
validateMpraObject <- function(object) {
    errors = character()
    if (any(dim(object@dnaCounts) != dim(object@rnaCounts))) {
        errors <- c(errors,
                    "DNA, RNA matrix must be of same dimensions")
    }
    if (is.null(rownames(object@dnaCounts)) |
        any(rownames(object@dnaCounts) != rownames(object@rnaCounts))) {
        errors <- c(errors,
                    "RNA, DNA feature names either missing or don't match")
    }
    
    if(length(errors) > 0) {
        return(errors)
    } else {
        return(TRUE)
    }
}

#' The main container class for MPRAnalyze
#' @slot dnaCounts matrix of DNA counts
#' @slot rnaCounts matrix of RNA counts
#' @slot colAnnot column annotations (info on condition, batch, barcode, etc)
#' @slot controls indices of negative controls
#' @slot controls.forfit indices of negative controls used for fitting in `full` 
#' mode
#' @slot lib.factor a factor with a level for each library, used for depth 
#' estimation and in `idr` quantitative mode
#' @slot dnaDepth library depth correction factors for DNA libraries
#' @slot rnaDepth library depth correction factors for RNA libraries
#' @slot rnaCtrlScale control-based correction factors for `scaled` modes
#' @slot mode analysis mode (identifies the type of analysis performed)
#' @slot model id of the distributional model used
#' @slot designs a Designs object containing the various design matrices used
#' @slot modelFits fitted models, populated by an analysis function
#' @slot modelFits.red fitted reduced models, populated by an LRT analysis function
#' @slot modelPreFits.dna.ctrl fitted models for control DNA
#' @slot BPPARAM The BiocParallel parallelization backend to use throughout
setClass("MpraObject", validity = validateMpraObject,
         slots = c(
             ## provided by user
             dnaCounts = "matrix",
             rnaCounts = "matrix",
             colAnnot = "data.frame",
             controls = "integerORNULL", 
             controls.forfit = "integerORNULL", 
             
             lib.factor = "factor",
             dnaDepth = "numeric",
             rnaDepth = "numeric",
             rnaCtrlScale = "numericORNULL",
             
             mode = "character",
             model = "character",
             designs = "Designs",
             modelFits = "list",
             modelFits.red = "list", 
             modelPreFits.dna.ctrl = "listORNULL",
             
             results = "listORNULL",
             
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
MpraObject <- function(dnaCounts, rnaCounts, colAnnot=NULL, controls=NA_integer_,
                       BPPARAM=NULL) {
    if(is.null(BPPARAM)) {
        BPPARAM <- SerialParam()
    }
    
    if(is.logical(controls)) {
        controls <- which(controls)
    }
    
    obj <- new("MpraObject", dnaCounts=dnaCounts, rnaCounts=rnaCounts,
               colAnnot=colAnnot, controls=controls, BPPARAM=BPPARAM)
    return(obj)
}


#' Get DNA model fits from an MpraObject
#' 
#' @param obj MpraObject to extract from
#' @param enhancers enhancer to extract 
#' @param depth include depth correction
#' @param full whether to extract from full model
#' 
#' @return DNA fits (numeric, enhancers x samples)
#' 
#' @export
getDNAFits <- function(obj, enhancers=NULL, depth=FALSE, full=TRUE){
    if(is.null(enhancers)) {
        enhancers <- names(obj@modelFits)
    }
    if(full == TRUE){
        fit <- obj@modelFits
    } else {
        fit <- obj@modelFits.red
    }
    
    coef.mat <- t(fit$d.coef[enhancers,,drop=FALSE])
    coef.mat[is.na(coef.mat)] <- 0
    
    if(obj@model == "gamma.pois") {
        # alpha / rate
        dfit <- exp(coef.mat[1,] + obj@designs@dna %*% coef.mat[-1,,drop=FALSE])
    } else if(obj@model == "ln.nb") {
        dfit <- exp(obj@designs@dna %*% coef.mat[-1,,drop=FALSE])
    }
    
    if(depth == TRUE){
        dfit <- dfit * matrix(obj@dnaDepth, 
                              nrow=NROW(dfit), ncol=NCOL(dfit), byrow=TRUE)
    }
    
    colnames(dfit) <- enhancers
    rownames(dfit) <- rownames(obj@colAnnot)
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
getRNAFits <- function(obj, enhancers=NULL, depth=TRUE, full=TRUE, rnascale=TRUE){
    if(is.null(enhancers)) {
        enhancers <- names(obj@modelFits)
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
    
    dfit <- getDNAFits(obj=obj, enhancers=enhancers, depth=FALSE, full=full)
    
    joint.des.mat <- cbind(rdesign, rctrldesign)
    
    coef.mat <- t(fit$r.coef[enhancers,-1,drop=FALSE])
    if(!is.null(obj@rnaCtrlScale) & rnascale) {
        coef.mat <- rbind(coef.mat, 
                          replicate(NCOL(coef.mat), obj@rnaCtrlScale))
    }
    rfit <- exp(joint.des.mat %*% coef.mat)
        
    if(depth == TRUE){
        rfit <- rfit * matrix(obj@rnaDepth, nrow = NROW(rfit), 
                              ncol = NCOL(rfit), byrow = TRUE)
    }
    
    rfit <- t(rfit) * dfit
    
    rownames(rfit) <- enhancers
    colnames(rfit) <- rownames(obj@colAnnot)
    return(rfit)
}

#' extract the DNA model parameters
#' @param obj the MpraObject to extract the parameters from
#' @param features the features to extract the parameters from (be default, 
#' parameters will be returned for all features)
#' @param full if TRUE (default), return the parameters of the full model. 
#' Otherwise, return the parameters of the reduced model (only relevant for 
#' LRT-based analyses)
#' @return a data.frame of features (rows) by parameters (cols). By convension, the
#' first parameter is related to the second moment, and the interpretation of 
#' it depends on the distributional model used (`alpha` for `gamma.pois`, variance 
#' for `ln.nb`)
#' @export
extractModelParameters.DNA <- function(obj, features=NULL, full=TRUE) {
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

#' extract the RNA model parameters
#' @param obj the MpraObject to extract the parameters from
#' @param features the features to extract the parameters from (be default, 
#' parameters will be returned for all features)
#' @param full if TRUE (default), return the parameters of the full model. 
#' Otherwise, return the parameters of the reduced odel
#' @return a data.frame of features (rows) by parameters (cols). By convension, the
#' first parameter is related to the second moment, and the interpretation of 
#' it depends on the distributional model used (`alpha` for `gamma.pois`,  
#' `psi`for `ln.nb`)
#' @export
extractModelParameters.RNA <- function(obj, features=NULL, full=TRUE) {
    if(is.null(obj@modelFits)){
        stop("can't extract model parameters before fitting a model. An analysis function must be called first.")
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
        stop("Parameters can't be extracted from reduced model, since analysis did not include fitting a reduced model")
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
        par.sdlog <-as.vector( exp(coef.mat[1,]))
        par.meanlog <- as.vector(exp(obj@designs@dna %*% coef.mat[-1,,drop=FALSE]))
        par <- data.frame(meanlog = par.meanlog, sdlog = par.sdlog)
    }
    
    rownames(par) <- rownames(obj@colAnnot)
    return(par)
}

#' @rdname getDistrParam
getDistrParam.RNA <- function(obj, enhancer=NULL, full=TRUE){
    
    if(full == TRUE){
        fit <- obj@modelFits
    } else {
        fit <- obj@modelFits.red
    }
    rfit <- as.vector(getRNAFits(obj, enhancers=enhancer, depth=FALSE, 
                                 full=full, rnascale=TRUE))
    
    if(obj@model == "gamma.pois") {
        par.size <- as.vector(exp(fit$r.coef[enhancer,1]))
        par.mu <- rfit
        par <- data.frame(size = par.size, mu = par.mu)
    } else if(obj@model == "ln.nb") {
        par.size <- as.vector(exp(fit$r.coef[enhancer,1]))
        par.mu <- rfit
        par <- data.frame(size = par.size, mu = par.mu)
    }
    
    rownames(par) <- rownames(obj@colAnnot)
    return(par)
}

#' Resample observations of enhancer from fit distribution
#' 
#' @param obj MpraObject to extract from
#' @param enhancer enhancer to extract 
#' @param full whether to extract from full model
#' 
#' @return resampled observations
#' 
#' @export
resampleObs <- function(obj, enhancer=NULL, full=TRUE){
    dpar <- getDistrParam.DNA(obj, enhancer=enhancer, full=full)
    rpar <- getDistrParam.RNA(obj, enhancer=enhancer, full=full)
    if(obj@model=="gamma.pois") {
        dsample <- apply(dpar, 1, function(x) {
            rgamma(n = 1, shape = x["shape"], rate = x["rate"])
        })
        rsample <- apply(rpar, 1, function(x) {
            rnbinom(n = 1, size = x["size"], mu = x["mu"])
        })
    } else if(obj@model=="ln.nb") {
        dsample <- apply(dpar, 1, function(x) {
            rlnorm(n = 1, meanlog = x["meanlog"], sdlog = x["sdlog"])
        })
        rsample <- apply(rpar, 1, function(x) {
            rnbinom(n = 1, size = x["size"], mu = x["mu"])
        })
    } 
    return(data.frame(dna=dsample, rna=rsample))
}