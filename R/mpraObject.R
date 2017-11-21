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

##TODO: document this
setClass("MpraObject", validity = validateMpraObject,
         slots = c(
             ## provided by user
             dnaCounts = "matrix",
             rnaCounts = "matrix",
             colAnnot = "data.frame",
             controls = "integerORNULL", #idx of negative controls
             controls.forfit = "integerORNULL", #idx of negative controls used in fitting of full and red
             
             dnaDepth = "numeric",
             rnaDepth = "numeric",
             rnaCtrlScale = "numericORNULL",
             
             mode = "character",
             model = "character",
             designs = "Designs",
             modelFits = "list",
             modelFits.red = "list", ##only used for LRT diff mode
             modelPreFits.dna.ctrl = "listORNULL", ##only used for LRT_iter diff mode 
             
             results = "data.frame",
             
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
                       BPPARAM=NULL, ncores=NULL) {
    if(is.null(BPPARAM) & is.null(ncores) | ncores == 1) {
        BPPARAM <- bpparam("SerialParam")
    } else if(is.null(BPPARAM) & !is.null(ncores)) {
        BPPARAM <- MulticoreParam(workers = ncores)
    } else if(!is.null(BPPARAM) & !is.null(ncores)) {
        stop("supply either BPPARAM or ncores but not both")
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
getDNAFits <- function(obj, enhancers, depth=FALSE, full=TRUE){
    if(full == TRUE){
        fit <- obj@modelFits
    } else {
        fit <- obj@modelFits.red
    }
    if(obj@model == "gamma.pois") {
        # mu_gamma = alpha / beta 
        #          = alpha * model_dna
        # shape_gamma = alpha
        # rate_gamma = beta = 1/model_dna
        dfit <- do.call(rbind, lapply(enhancers, function(i) {
            exp(fit[[i]]$d.coef[1] + fit[[i]]$d.coef[-1] %*% t(obj@designs@dna))
        }))
    } else if(obj@model == "ln.nb") {
        dfit <- do.call(rbind, lapply(enhancers, function(i) {
            exp(fit[[i]]$d.coef[-1] %*% t(obj@designs@dna))
        }))
    }
    if(depth == TRUE){
        dfit <- dfit * obj@dnaDepth
    }
    return(dfit)
}

#' Get RNA full model fits from an MpraObject
#' 
#' @param obj MpraObject to extract from
#' @param enhancers enhancer to extract 
#' @param depth include depth correction
#' @param full whether to extract from full model
#' 
#' @return RNA fits (numeric, enhancers x samples)
#' 
#' @export
getRNAFits <- function(obj, enhancers, depth=TRUE, full=TRUE){
    if(full == TRUE){
        fit <- obj@modelFits
        rdesign <- obj@designs@rnaFull
        rctrldesign <- obj@designs@rnaCtrlFull
    } else {
        fit <- obj@modelFits.red
        rdesign <- obj@designs@rnaRed
        rctrldesign <- obj@designs@rnaCtrlRed
    }
    dfit <- getDNAFits(obj=obj, enhancers=enhancers, 
                       depth=FALSE, full=full)
    if(obj@model == "gamma.pois") {
        # mu_NB = alpha / beta * rna_model
        #       = alpha * dna_model * rna_model
        #       = mu_rna * rna_model
        # size_alpha = alpha
        rfit <- do.call(rbind, lapply(enhancers, function(i) {
            dfit*exp( c(fit[[i]]$r.coef, obj@rnaCtrlScale) %*% t(cbind(rdesign, rctrldesign)) )
        }))
    } else if(obj@model == "ln.nb") {
        # mu_NB = rna_model
        # Note that the first parameter in r.coef is the variance parameter and that
        # the design matrix is padded with 0s at this position.
        rfit <- do.call(rbind, lapply(enhancers, function(i) {
            dfit*exp( c(fit[[i]]$r.coef, obj@rnaCtrlScale) %*% t(cbind(rdesign, rctrldesign)) )
        }))
    }
    if(depth == TRUE){
        rfit <- rfit * matrix(obj@rnaDepth, nrow = NROW(rfit), 
                              ncol = NCOL(rfit), byrow = TRUE)
    }
    return(rfit)
}
