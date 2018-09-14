#' Get DNA model-based estimates from an MpraObject (the expected values based 
#' on the model). These can be compared with the observed counts to assess 
#' goodness of fit.
#' 
#' @param obj MpraObject to extract from
#' @param enhancers which enhancers to get the fits for. Can be character 
#' vectors with enhancer names, logical or numeric enhancer indices, or NULL if 
#' all enhancers are to be extracted (default)
#' @param depth include depth correction in the model fitting (default TRUE)
#' @param full if LRT modeling was used, TRUE (default) would return the fits
#' of the full model, FALSEwould return the reduced model fits.
#' @param transition use the DNA->RNA transition matrix (deafult: FALSE). This 
#' is useful if the DNA observations need to be distributed to match the RNA
#' observations.
#' 
#' @return DNA fits (numeric, enhancers x samples)
#' 
#' @export
#' 
#' @examples
#' data <- simulateMPRA(tr = rep(2,5), da=NULL, nbatch=2, nbc=15)
#' obj <- MpraObject(dnaCounts = data$obs.dna, 
#'                   rnaCounts = data$obs.rna, 
#'                   colAnnot = data$annot)
#' obj <- estimateDepthFactors(obj, lib.factor = "batch", which.lib = "both")
#' obj <- analyzeQuantification(obj, dnaDesign = ~ batch + barcode, 
#'                               rnaDesign = ~1)
#' dna.fits <- getFits_DNA(obj)
getFits_DNA <- function(obj, enhancers=NULL, depth=TRUE, full=TRUE, 
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
    
    coef.mat <- fit$d.coef[enhancers,,drop=FALSE]
    coef.mat[is.na(coef.mat)] <- 0
    
    dmat <- obj@designs@dna
    if(transition) {
        dmat <- obj@designs@dna2rna
    }
    
    if(model(obj) == "gamma.pois") {
        # disp * (1/beta)
        disp <- rep(coef.mat[,1], NROW(dmat))
        beta.inv <- coef.mat[,-1] %*% t(dmat)
        dfit <- exp(disp + beta.inv)
        
    } else if(model(obj) == "ln.nb" | model(obj) == "ln.ln") {
        dfit <- exp(coef.mat[-1,,drop=FALSE] %*% t(dmat))
    }
    
    if(depth == TRUE){
        dfit <- dfit * rep(dnaDepth(obj), each = NROW(coef.mat))
    }
    
    rownames(dfit) <- enhancers
    if(transition) {
        colnames(dfit) <- rownames(rnaAnnot(obj))
    } else {
        colnames(dfit) <- rownames(dnaAnnot(obj))
    }
    return(dfit)
}

#' Get RNA model-based estimates from an MpraObject (the expected values based 
#' on the model). These can be compared with the observed counts to assess 
#' goodness of fit.
#' 
#' @param obj MpraObject to extract from
#' @param enhancers which enhancers to get the fits for. Can be character 
#' vectors with enhancer names, logical or numeric enhancer indices, or NULL if 
#' all enhancers are to be extracted (default)
#' @param depth include depth correction in the model fitting (default TRUE)
#' @param full if LRT modeling was used, TRUE (default) would return the fits
#' of the full model, FALSEwould return the reduced model fits.
#' @param rnascale if controls were used to correct the fitting (in comparative
#' analyses), use these factors to re-adjust the estimates back.
#' 
#' @return RNA fits (numeric, enhancers x samples)
#' 
#' @export
#' @examples
#' data <- simulateMPRA(tr = rep(2,5), da=NULL, nbatch=2, nbc=15)
#' obj <- MpraObject(dnaCounts = data$obs.dna, 
#'                   rnaCounts = data$obs.rna, 
#'                   colAnnot = data$annot)
#' obj <- estimateDepthFactors(obj, lib.factor = "batch", which.lib = "both")
#' obj <- analyzeQuantification(obj, dnaDesign = ~ batch + barcode, 
#'                               rnaDesign = ~1)
#' rna.fits <- getFits_RNA(obj)
getFits_RNA <- function(obj, enhancers=NULL, depth=TRUE, full=TRUE, 
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
    
    dfit <- getFits_DNA(obj=obj, enhancers=enhancers, depth=FALSE, 
                        full=full, transition=TRUE)
    
    joint.des.mat <- cbind(rdesign, rctrldesign)
    
    coef.mat <- t(fit$r.coef[enhancers,-1,drop=FALSE])
    if(!is.null(obj@rnaCtrlScale) & rnascale) {
        coef.mat <- rbind(coef.mat,
                        replicate(NCOL(coef.mat), obj@rnaCtrlScale))
    }
    rfit <- exp(joint.des.mat %*% coef.mat)
    
    if(depth == TRUE){
        rfit <- rfit * replicate(NCOL(rfit), rnaDepth(obj))
    }
    
    rfit <- t(rfit) * dfit
    
    rownames(rfit) <- enhancers
    colnames(rfit) <- rownames(rnaAnnot(obj))
    return(rfit)
}

#' extract the DNA model parameters
#' 
#' @rdname extractModelParameters
#' @aliases extractModelParameters_DNA
#' extractModelParameters_RNA
#' 
#' @param obj the MpraObject to extract the parameters from
#' @param features the features to extract the parameters from (by default, 
#' parameters will be returned for all features)
#' @param full if TRUE (default), return the parameters of the full model. 
#' Otherwise, return the parameters of the reduced model (only relevant for 
#' LRT-based analyses)
#' @return a data.frame of features (rows) by parameters (cols). By convension,
#' the first parameter is related to the second moment, and the interpretation 
#' of  it depends on the distributional model used (`alpha` for `gamma.pois`, 
#' variance  for `ln.nb` and `ln.ln`)
#' 
#' @examples
#' data <- simulateMPRA(tr = rep(2,5), da=NULL, nbatch=2, nbc=15)
#' obj <- MpraObject(dnaCounts = data$obs.dna, 
#'                   rnaCounts = data$obs.rna, 
#'                   colAnnot = data$annot)
#' obj <- estimateDepthFactors(obj, lib.factor = "batch", which.lib = "both")
#' obj <- analyzeQuantification(obj, dnaDesign = ~ batch + barcode, 
#'                               rnaDesign = ~1)
#' model.params.dna <- getModelParameters_DNA(obj)
#' model.params.rna <- getModelParameters_RNA(obj)
#' 
#' @export
#' @rdname extractModelParameters
getModelParameters_DNA <- function(obj, features=NULL, full=TRUE) {
    if(is.null(features)) {
        features <- seq_len(NROW(dnaCounts(obj)))
    } else if (is.character(features)) {
        features <- which(rownames(dnaCounts(obj)) %in% features)
    }
    
    if(is.null(obj@modelFits)){
        stop("can't extract model parameters before fitting a model. An \
             analysis function must be called first.")
    }
    if(full) {
        coef.mat <- obj@modelFits$d.coef[features,,drop=FALSE] 
    } else if (!full & !is.null(obj@modelFits.red)) {
        coef.mat <- obj@modelFits.red$d.coef[features,,drop=FALSE] 
    } else {
        stop("Parameters can't be extracted from reduced model, since analysis \
             did not include fitting a reduced model")
    }
    colnames(coef.mat) <- c("disp", colnames(obj@designs@dna))
    rownames(coef.mat) <- rownames(dnaCounts(obj))[features]
    return(as.data.frame(coef.mat))
}

#' @export
#' @rdname extractModelParameters
getModelParameters_RNA <- function(obj, features=NULL, full=TRUE) {
    if(is.null(obj@modelFits)){
        stop("can't extract model parameters before fitting a model")
    }
    if(is.null(features)) {
        features <- seq_len(NROW(dnaCounts(obj)))
    } else if (is.character(features)) {
        features <- which(rownames(dnaCounts(obj)) %in% features)
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
    
    rownames(coef.mat) <- rownames(dnaCounts(obj))[features]
    return(as.data.frame(coef.mat))
}

#' Get model distribution parameters from an MpraObject of a given candidate 
#' enhancer
#' 
#' @rdname getDistrParam
#' @aliases getDistrParam_DNA
#' getDistrParam_RNA
#' 
#' @param obj MpraObject to extract from
#' @param enhancer enhancer to extract 
#' @param full whether to extract from full model
#' 
#' @return fit parameters (numeric, samples x parameters)
#' 
#' @export
#' 
#' @examples
#' data <- simulateMPRA(tr = rep(2,5), da=NULL, nbatch=2, nbc=15)
#' obj <- MpraObject(dnaCounts = data$obs.dna, 
#'                   rnaCounts = data$obs.rna, 
#'                   colAnnot = data$annot)
#' obj <- estimateDepthFactors(obj, lib.factor = "batch", which.lib = "both")
#' obj <- analyzeQuantification(obj, dnaDesign = ~ batch + barcode, 
#'                               rnaDesign = ~1)
#' ## get distributional parameters of the first enhancer:
#' dist.params.dna <- getDistrParam_DNA(obj, 1)
#' dist.params.rna <- getDistrParam_RNA(obj, 1)
#' 
#' @rdname getDistrParam
getDistrParam_DNA <- function(obj, enhancer, full=TRUE){
    if(full == TRUE){
        fit <- obj@modelFits
    } else {
        fit <- obj@modelFits.red
    }
    
    coef.mat <- t(fit$d.coef[enhancer,,drop=FALSE])
    coef.mat[is.na(coef.mat)] <- 0
    
    if(model(obj) == "gamma.pois") {
        par.shape <- as.vector(exp(coef.mat[1,]))
        par.rate <- as.vector(exp(-obj@designs@dna %*% 
                                      coef.mat[-1,,drop=FALSE]))
        par <- data.frame(shape = par.shape, rate = par.rate)
    } else if(model(obj) == "ln.nb" | model(obj) == "ln.ln") {
        par.sdlog <- as.vector(exp(coef.mat[1,]))
        par.meanlog <- as.vector(obj@designs@dna %*% coef.mat[-1,,drop=FALSE])
        par <- data.frame(meanlog = par.meanlog, sdlog = par.sdlog)
    }
    
    rownames(par) <- rownames(dnaAnnot(obj))
    return(par)
}

#' @export
#' @rdname getDistrParam
getDistrParam_RNA <- function(obj, enhancer=NULL, full=TRUE){
    
    if(full == TRUE){
        fit <- obj@modelFits
    } else {
        fit <- obj@modelFits.red
    }
    rfit <- as.vector(getFits_RNA(obj, enhancers=enhancer, depth=FALSE, 
                                full=full, rnascale=TRUE))
    
    if(model(obj) == "gamma.pois") {
        par.size <- as.vector(exp(fit$r.coef[enhancer,1]))
        par.mu <- rfit
        par <- data.frame(size = par.size, mu = par.mu)
    } else if(model(obj) == "ln.nb") {
        par.size <- as.vector(exp(fit$r.coef[enhancer,1]))
        par.mu <- rfit
        par <- data.frame(size = par.size, mu = par.mu)
    } else if(model(obj) == "ln.ln") {
        par.sdlog <- as.vector(exp(fit$r.coef[enhancer,1]))
        par.meanlog <- log(rfit)
        par <- data.frame(meanlog = par.meanlog, sdlog = par.sdlog)
    }
    
    rownames(par) <- rownames(rnaAnnot(obj))
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
#' design. Will return the corresponding rates for all values of the given 
#' factor
#' @param full if true, return rate of the full model (default), otherwise of
#' the reduced model (only applies if an LRT-based analysis was used)
#' @export
#' @return the estimate for transcription rate as fitted by the model
#' 
#' @examples
#' data <- simulateMPRA(tr = rep(2,10), da=c(rep(2,5), rep(2.5,5)), 
#'                      nbatch=2, nbc=15)
#' obj <- MpraObject(dnaCounts = data$obs.dna, 
#'                   rnaCounts = data$obs.rna, 
#'                   colAnnot = data$annot)
#' obj <- estimateDepthFactors(obj, lib.factor = "batch", which.lib = "both")
#' obj <- analyzeComparative(obj, dnaDesign = ~ batch + barcode + condition, 
#'                               rnaDesign = ~ condition, reducedDesign = ~ 1)
#' ## get alpha estimate for the two conditions
#' alpha <- getAlpha(obj, by.factor="condition")
getAlpha <- function(obj, by.factor=NULL, full=TRUE) {
    des <- obj@designs@rnaFull
    if(!full) {
        des <- obj@designs@rnaRed
        if(is.null(des)) {
            stop("reduced model not available")
        }
    }
    
    if(is.null(by.factor)) {
        return(data.frame(alpha=getSingleAlpha(obj, full=full)))
    } else if (by.factor=="all") {
        ## extract all parametres except for the dispersion
        alpha.mat <- getModelParameters_RNA(obj, full=full)[,-1]
        if(checkForIntercept(des)) {
            ## add the intercept to all other columns
            alpha.mat[,-1] <- alpha.mat[,-1] + alpha.mat[,1]
        }
        return(as.data.frame(exp(alpha.mat)))
    } else {
        l <- levels(rnaAnnot(obj)[,by.factor])
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
        colnames(alpha.mat) <- levels(rnaAnnot(obj)[,by.factor])
        
        return(as.data.frame(alpha.mat))
    }
}