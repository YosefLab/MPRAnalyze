
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
             
             lib.factor = "factor",
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
    if(is.null(BPPARAM) & is.null(ncores)) {
        BPPARAM <- SerialParam()
    } else if(is.null(BPPARAM) & !is.null(ncores)) {
        BPPARAM <- MulticoreParam(workers = ncores)
    } else if(!is.null(BPPARAM) & !is.null(ncores)) {
        stop("supply either BPPARAM or ncores but not both")
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
    if(obj@model == "gamma.pois") {
        # mu_gamma = alpha / beta 
        #          = alpha * model_dna
        # shape_gamma = alpha
        # rate_gamma = beta = 1/model_dna
        dfit <- do.call(cbind, lapply(enhancers, function(i) {
            exp(fit[[i]]$d.coef[1] + obj@designs@dna %*% fit[[i]]$d.coef[-1])
        }))
    } else if(obj@model == "ln.nb") {
        dfit <- do.call(cbind, lapply(enhancers, function(i) {
            exp(obj@designs@dna %*% fit[[i]]$d.coef[-1])
        }))
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
    dfit <- getDNAFits(obj=obj, enhancers=enhancers, 
                       depth=FALSE, full=full)
    
    joint.des.mat <- cbind(rdesign, rctrldesign)
    rfit <- do.call(cbind, lapply(enhancers, function(i) {
        exp(joint.des.mat %*% c(fit[[i]]$r.coef[-1], rctrlscale))
    }))
        
    if(depth == TRUE){
        rfit <- rfit * matrix(obj@rnaDepth, nrow = NROW(rfit), 
                              ncol = NCOL(rfit), byrow = TRUE)
    }
    
    rfit <- t(rfit) * dfit
    
    rownames(rfit) <- enhancers
    colnames(rfit) <- rownames(obj@colAnnot)
    return(rfit)
}
