setClassUnion('numericORNULL', members=c('numeric', 'NULL'))
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
    controls = "integer", #idx of negative controls

    dnaDepth = "numeric",
    rnaDepth = "numeric",
    rnaCtrlScale = "numericORNULL",
    
    mode = "character",
    model = "character",
    designs = "Designs",
    modelFits = "list",
    modelFits.red = "list", ##only used for LRT diff mode
    modelPreFits.dna.ctrl = "list", ##only used for LRT_iter diff mode 

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
                       BPPARAM=NULL) {
    if(is.null(BPPARAM)) {
        BPPARAM <- bpparam()
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
getDNAFits <- function(obj, enhancers, depth=TRUE, full=TRUE){
    if(full == TRUE){
        fit <- obj@modelFits
    } else {
        fit <- obj@modelFits.red
    }
    if(obj@model == "gamma.poisson") {
        dfit <- do.call(rbind, lapply(enhancers, function(i) {
            fit[[i]]$d.par[1] + 
                exp(sum(fit[[i]]$d.par[-1] * t(obj@designs$dna[i,]), na.rm = TRUE))
        }))
    } else if(obj@model == "ln.nb") {
        dfit <- do.call(rbind, lapply(enhancers, function(i) {
            exp(sum(fit[[i]]$d.par[-1] * t(obj@designs$dna[i,]), na.rm = TRUE))
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
getRNAFits <- function(obj, enhancers, depth=TRUE, full=TRUE){
    if(full == TRUE){
        fit <- obj@modelFits
        rdesign <- obj@designs$rnaFull
    } else {
        fit <- obj@modelFits.red
        rdesign <- obj@designs$rnaRed
    }
    dfit <- getDNAFits(obj=obj, enhancers=enhancers, 
                       depth=FALSE, full=full)
    if(obj@model == "gamma.poisson") {
        rfit <- dfit * do.call(rbind, lapply(enhancers, function(i) {
                exp(sum(fit[[i]]$r.par * t(rdesign[i,]), na.rm = TRUE))
        }))
    } else if(obj@model == "ln.nb") {
        rfit <- dfit * do.call(rbind, lapply(enhancers, function(j) {
            exp(sum(fit[[i]]$r.par[-1] * t(rdesign[i,]), na.rm = TRUE))
        }))
    }
    if(depth == TRUE){
        rfit <- rfit * obj@rnaDepth
    }
    return(rfit)
}
