
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
    rnaCtrlRed = "Design",
    
    dna2rna = "Design"
))

#' validate an MpraObject
#'
#' @param object the MpraObject instance to validate
#'
#' @return TRUE if object is valid, otherwise vector of character strings
#' explaining how it's invalid
#' @noRd
validateMpraObject <- function(object) {
    errors = character()
    if (NROW(dnaCounts(object)) != NROW(rnaCounts(object))) {
        errors <- c(errors, "DNA and RNA don't have the same number of Rows")
    }
    if (NCOL(dnaCounts(object)) != NROW(dnaAnnot(object))) {
        errors <- c(errors, "DNA observations and annotations don't match")
    }
    if (NCOL(rnaCounts(object)) != NROW(rnaAnnot(object))) {
        errors <- c(errors, "RNA observations and annotations don't match")
    }
    if (is.null(rownames(dnaCounts(object))) | 
        is.null(rownames(rnaCounts(object))) |
        any(rownames(dnaCounts(object)) != rownames(rnaCounts(object)))) {
        errors <- c(errors,
                    "RNA, DNA feature names either missing or don't match")
    }
    if (length(unique(rownames(dnaCounts(object)))) != 
        length(rownames(dnaCounts(object)))) {
        errors <- c(errors, "enhancer IDs must by unique")
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
#' @slot dnaAnnot DNA column annotations 
#' (info on condition, batch, barcode, etc)
#' @slot rnaAnnot RNA column annotations 
#' (info on condition, batch, barcode, etc)
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
#' @slot modelFits.red fitted reduced models, populated by an LRT analysis 
#' function
#' @slot modelPreFits.dna.ctrl fitted models for control DNA
#' @slot BPPARAM The BiocParallel parallelization backend to use throughout
setClass("MpraObject", validity = validateMpraObject,
         slots = c(
             ## provided by user
             dnaCounts = "matrix",
             rnaCounts = "matrix",
             dnaAnnot = "data.frame",
             rnaAnnot = "data.frame",
             controls = "integerORNULL", 
             controls.forfit = "integerORNULL", 
             
             lib.factor = "factor",
             dnaDepth = "numeric",
             rnaDepth = "numeric",
             rnaCtrlScale = "numericORNULL",
             
             model = "character",
             designs = "Designs",
             modelFits = "list",
             modelFits.red = "list", 
             modelPreFits.dna.ctrl = "listORNULL",
             
             BPPARAM = "BiocParallelParam"
         ))


#' Initialize a MpraObject object
#'
#' @import BiocParallel
#'
#' @param dnaCounts the DNA counts matrix
#' @param rnaCounts the RNA counts matrix
#' @param dnaAnnot column annotations for the DNA matrix
#' @param rnaAnnot column annotations for the RNA matrix
#' @param colAnnot (deprecated) column annotations - should only be used if DNA,
#' RNA annotations are the same.
#' @param controls a vector specifying which enhancers are negative controls
#' (scrambles)
#' @param BPPARAM the biocParalell backend to use for parallelization throughout
#' the analysis
#' @return an MpraObject
#' @export
#' @examples
#' data <- simulateMPRA(tr = rep(2,10), da=c(rep(2,5), rep(2.5,5)), 
#'                      nbatch=2, nbc=20)
#' ## use 3 of the non-active enhancers as controls
#' obj <- MpraObject(dnaCounts = data$obs.dna, 
#'                   rnaCounts = data$obs.rna, 
#'                   colAnnot = data$annot,
#'                   controls = as.integer(c(1,2,4)))
MpraObject <- function(dnaCounts, rnaCounts, dnaAnnot=NULL, rnaAnnot=NULL, 
                       colAnnot=NULL, controls=NA_integer_,
                       BPPARAM=NULL) {
    if(is.null(BPPARAM)) {
        BPPARAM <- SerialParam()
    }
    
    if(is.logical(controls)) {
        controls <- which(controls)
    } else if (is.character(controls)) {
        controls <- which(rownames(dnaCounts) %in% controls)
    }
    if((is.null(dnaAnnot) | is.null(rnaAnnot)) & !is.null(colAnnot)) {
        rnaAnnot <- dnaAnnot <- colAnnot
    }
    
    ## remove invalid enhancers: either all dna or all rna counts are 0.
    invalid <- union(which(apply(dnaCounts, 1, function(x) all(x==0))),
                     which(apply(rnaCounts, 1, function(x) all(x==0))))
    if(length(invalid) > 0) {
        warning(length(invalid), " enhancers were removed from the analysis")
        if(length(controls) > 1) {
            ctrl <- rep(FALSE, NROW(dnaCounts))
            ctrl[controls] <- TRUE
            controls <- which(ctrl[-invalid])
        }
        dnaCounts <- dnaCounts[-invalid,]
        rnaCounts <- rnaCounts[-invalid,]
    }
    
    obj <- new("MpraObject", dnaCounts=dnaCounts, rnaCounts=rnaCounts,
               dnaAnnot=dnaAnnot, rnaAnnot=rnaAnnot, controls=controls, 
               BPPARAM=BPPARAM)
    return(obj)
}

#' get properties from a given MpraObject
#' @param obj the MpraObject to extract the property from
#' @return The value of the property
#' @export
#' @rdname MpraObject_getters
#' @examples
#' data <- simulateMPRA(tr = rep(2,15), nbatch=2, nbc=15)
#' obj <- MpraObject(dnaCounts = data$obs.dna, 
#'                   rnaCounts = data$obs.rna, 
#'                   colAnnot = data$annot,
#'                   controls = 1:5)
#' dnaCounts <- dnaCounts(obj)
#' rnaCounts <- rnaCounts(obj)
#' dnaAnnot <- dnaAnnot(obj)
#' rnaAnnot <- rnaAnnot(obj)
#' controls <- controls(obj)
#' dnaDepth <- dnaDepth(obj)
#' rnaDepth <- rnaDepth(obj)
#' model <- model(obj)
setGeneric("dnaCounts", function(obj) standardGeneric("dnaCounts"))

#' @rdname MpraObject_getters
setMethod("dnaCounts", signature(obj="MpraObject"), function(obj) obj@dnaCounts)

#' @rdname MpraObject_getters
#' @export
setGeneric("rnaCounts", function(obj) standardGeneric("rnaCounts"))

#' @rdname MpraObject_getters
setMethod("rnaCounts", signature(obj="MpraObject"), function(obj) obj@rnaCounts)

#' @rdname MpraObject_getters
#' @export
setGeneric("dnaAnnot", function(obj) standardGeneric("dnaAnnot"))

#' @rdname MpraObject_getters
setMethod("dnaAnnot", signature(obj="MpraObject"), function(obj) obj@dnaAnnot)

#' @rdname MpraObject_getters
#' @export
setGeneric("rnaAnnot", function(obj) standardGeneric("rnaAnnot"))

#' @rdname MpraObject_getters
setMethod("rnaAnnot", signature(obj="MpraObject"), function(obj) obj@rnaAnnot)

#' @rdname MpraObject_getters
#' @export
setGeneric("controls", function(obj) standardGeneric("controls"))

#' @rdname MpraObject_getters
setMethod("controls", signature(obj="MpraObject"), function(obj) obj@controls)

#' @rdname MpraObject_getters
#' @export
setGeneric("dnaDepth", function(obj) standardGeneric("dnaDepth"))

#' @rdname MpraObject_getters
setMethod("dnaDepth", signature(obj="MpraObject"), function(obj) obj@dnaDepth)

#' @rdname MpraObject_getters
#' @export
setGeneric("rnaDepth", function(obj) standardGeneric("rnaDepth"))

#' @rdname MpraObject_getters
setMethod("rnaDepth", signature(obj="MpraObject"), function(obj) obj@rnaDepth)

#' @rdname MpraObject_getters
#' @export
setGeneric("model", function(obj) standardGeneric("model"))

#' @rdname MpraObject_getters
setMethod("model", signature(obj="MpraObject"), function(obj) obj@model)