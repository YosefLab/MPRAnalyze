
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

#' Full documentation of MpraObject structure and slots
#' 
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
#' 
#' @noRd
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

#' MpraObject
#' 
#' @description The main object MPRAnalyze works with, contains the input data,
#' associated annotations, model parameters and analysis results.
#' 
#' @export
#' @import BiocParallel
#' @importFrom SummarizedExperiment SummarizedExperiment assay colData
#' 
#' @param dnaCounts the DNA count matrix, or a SummarizedExperiment object
#' containing the DNA Counts and column annotations for the DNA data. If the
#' input is a SummarizedExperiment object, the dnaAnnot (or colAnnot) arguments
#' will be ignored
#' @param rnaCounts the RNA count matrix, or a SummarizedExperiment object
#' containing the RNA Counts and column annotations for the RNA data. If the
#' input is a SummarizedExperiment object, the rnaAnnot (or colAnnot) arguments
#' will be ignored
#' @param dnaAnnot data.frame with the DNA column (sample) annotations
#' @param rnaAnnot data.frame with the RNA column (sample) annotations
#' @param colAnnot if annotations for DNA and RNA are identical, they can be set
#' at the same time using colAnnot instead of using both rnaAnnot and dnaAnnot
#' @param controls IDs of the rows in the matrices that correspond to negative 
#' control enhancers. These are used to establish the null for quantification
#' purposes, and to correct systemic bias in comparative analyses. Can be a 
#' character vectors (corresponding to rownames in the data matrices), logical
#' or numeric indices.
#' @param BPPARAM a parallelization backend using the BiocParallel package, see
#' more details [here](http://bioconductor.org/packages/release/bioc/html/BiocParallel.html)
#' 
#' @param obj The MpraObject to extract properties from
#' 
#' @return an initialized MpraObject
#' 
#' @section Accessors:
#' MpraObject properties can be accessed using accessor functions
#' \describe{
#' \item{dnaCounts}{the DNA count matrix}
#' \item{rnaCounts}{the RNA count matrix}
#' \item{dnaAnnot}{data.frame with the DNA column (sample) annotations}
#' \item{ranAnnot}{data.frame with the RNA column (sample) annotations}
#' \item{model}{the distributional model used. the Gamma-Poisson convolutional 
#' model is used by default. see \code{\link{setModel}}}
#' \item{dnaDepth}{The library size correction factors computed for the DNA
#' libraries. These are computed by the \code{\link{estimateDepthFactors}} 
#' function and can be set manually using the \code{\link{setDepthFactors}}
#' function}
#' \item{rnaDepth}{The library size correction factors computed for the RNA
#' libraries These are computed by the \code{\link{estimateDepthFactors}} 
#' function and can be set manually using the \code{\link{setDepthFactors}}
#' function}
#' }
#' 
#' @examples
#' data <- simulateMPRA(tr = rep(2,10), da=c(rep(2,5), rep(2.5,5)), 
#'                      nbatch=2, nbc=20)
#' ## use 3 of the non-active enhancers as controls
#' obj <- MpraObject(dnaCounts = data$obs.dna, 
#'                   rnaCounts = data$obs.rna, 
#'                   colAnnot = data$annot,
#'                   controls = as.integer(c(1,2,4)))
#' ## alternatively, initialize the object with SummarizedExperiment objects:
#' \dontrun{
#' se.DNA <- SummarizedExperiment(list(data$obs.dna), colData=data$annot)
#' se.RNA <- SummarizedExperiment(list(data$obs.rna), colData=data$annot)
#' obj <- MpraObject(dnaCounts = se.DNA, rnaCounts = rna.se, 
#'                   controls = as.integer(c(1,2,4)))
#' }
#' dnaCounts <- dnaCounts(obj)
#' rnaCounts <- rnaCounts(obj)
#' dnaAnnot <- dnaAnnot(obj)
#' rnaAnnot <- rnaAnnot(obj)
#' controls <- controls(obj)
#' model <- model(obj)
#' 
#' obj <- estimateDepthFactors(obj, lib.factor=c("batch", "condition"))
#' dnaDepth <- dnaDepth(obj)
#' rnaDepth <- rnaDepth(obj)
#' 
#' @rdname MpraObject
setGeneric("MpraObject", 
           function(dnaCounts, rnaCounts, dnaAnnot=NULL, rnaAnnot=NULL, 
                    colAnnot=NULL, controls=NA_integer_,
                    BPPARAM=NULL) standardGeneric("MpraObject"))

#' @rdname MpraObject
#' @export
setMethod("MpraObject", signature = signature(dnaCounts = "matrix"),
          function(dnaCounts, rnaCounts, dnaAnnot=NULL, rnaAnnot=NULL, 
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
              
              ## remove invalid enhancers: all dna or all rna counts are 0.
              invalid <- union(which(apply(dnaCounts, 1, 
                                           function(x) all(x==0))),
                               which(apply(rnaCounts, 1, 
                                           function(x) all(x==0))))
              if(length(invalid) > 0) {
                  warning(length(invalid), 
                          " enhancers were removed from the analysis")
                  if(length(controls) > 1) {
                      ctrl <- rep(FALSE, NROW(dnaCounts))
                      ctrl[controls] <- TRUE
                      controls <- which(ctrl[-invalid])
                  }
                  dnaCounts <- dnaCounts[-invalid,]
                  rnaCounts <- rnaCounts[-invalid,]
              }
              
              obj <- new("MpraObject", dnaCounts=dnaCounts, 
                         rnaCounts=rnaCounts, dnaAnnot=dnaAnnot, 
                         rnaAnnot=rnaAnnot, controls=controls, 
                         BPPARAM=BPPARAM)
              return(obj)
              
})

#' @rdname MpraObject
#' @export
setMethod("MpraObject", 
          signature = signature(dnaCounts = "SummarizedExperiment"),
          function(dnaCounts, rnaCounts, dnaAnnot=NULL, rnaAnnot=NULL, 
                   colAnnot=NULL, controls=NA_integer_,
                   BPPARAM=NULL) {
              return(MpraObject(dnaCounts = assay(dnaCounts),
                                rnaCounts = assay(rnaCounts),
                                dnaAnnot = as.data.frame(colData(dnaCounts)),
                                rnaAnnot = as.data.frame(colData(rnaCounts)),
                                controls = controls,
                                BPPARAM = BPPARAM))
          })



#' @rdname MpraObject
setGeneric("dnaCounts", function(obj) standardGeneric("dnaCounts"))

#' @rdname MpraObject
#' @export
setMethod("dnaCounts", signature(obj="MpraObject"), function(obj) obj@dnaCounts)

#' @rdname MpraObject
setGeneric("rnaCounts", function(obj) standardGeneric("rnaCounts"))

#' @rdname MpraObject
#' @export
setMethod("rnaCounts", signature(obj="MpraObject"), function(obj) obj@rnaCounts)

#' @rdname MpraObject
setGeneric("dnaAnnot", function(obj) standardGeneric("dnaAnnot"))

#' @rdname MpraObject
#' @export
setMethod("dnaAnnot", signature(obj="MpraObject"), function(obj) obj@dnaAnnot)

#' @rdname MpraObject
setGeneric("rnaAnnot", function(obj) standardGeneric("rnaAnnot"))

#' @rdname MpraObject
#' @export
setMethod("rnaAnnot", signature(obj="MpraObject"), function(obj) obj@rnaAnnot)

#' @rdname MpraObject
setGeneric("controls", function(obj) standardGeneric("controls"))

#' @rdname MpraObject
#' @export
setMethod("controls", signature(obj="MpraObject"), function(obj) obj@controls)

#' @rdname MpraObject
setGeneric("dnaDepth", function(obj) standardGeneric("dnaDepth"))

#' @rdname MpraObject
#' @export
setMethod("dnaDepth", signature(obj="MpraObject"), function(obj) obj@dnaDepth)

#' @rdname MpraObject
setGeneric("rnaDepth", function(obj) standardGeneric("rnaDepth"))

#' @rdname MpraObject
#' @export
setMethod("rnaDepth", signature(obj="MpraObject"), function(obj) obj@rnaDepth)

#' @rdname MpraObject
setGeneric("model", function(obj) standardGeneric("model"))

#' @rdname MpraObject
#' @export
setMethod("model", signature(obj="MpraObject"), function(obj) obj@model)