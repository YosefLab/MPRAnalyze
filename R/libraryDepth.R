#' easy access container for depth estimation functions
#' @noRd
DEPTH_EST_FUNCTIONS = list(
    uq = function(x) {
        x[x==0] <- NA
        return(quantile(x, .75, na.rm=TRUE))
        },
    rle = function(x) exp(mean(log(x), na.rm=TRUE)),
    totsum = function(x) sum((x), na.rm = TRUE)
)


#' estimate library size correction factors
#'
#' @note since in most MPRA experiments multiple barcodes exist within a single
#' library, each column in the matrix is usually not a separate library. For
#' this reason, it is recommended to supply this function with the appropriate
#' partitioning of the data matrix columns into libraries, see lib.factor
#'
#' @param obj the MpraObject
#' @param which.lib which library to compute the depth factors for. Options are
#' "both" (default), "dna" or "rna". If the DNA and RNA counts have different 
#' library factors, this function should be called twice: once with "dna" and
#' once with "rna"
#' @param lib.factor the factor associating each sample to a library. Can be a
#' factor or the name of a column in the object's colAnnot. If not provided, the
#' data is assumed to have been generated from a single library, and constant
#' library depth is set.
#' @param depth.estimator a character indicating which depth estimation to use,
#' or a function to perform the estimation. Currently supported values are "uq" 
#' for upper quantile of non-zero values (default), "rle" for RLE (uses 
#' geometric mean, and is therefore not recommended if libraries have 0 counts),
#'  or "totsum" for total sum.
#' For a function input: function should take a numeric vector and return a 
#' single numeric, and preferably handle NA values. See examples.
#'
#' @return the MpraObject with estimated values for sequencing depth factors
#' 
#' @export
#' 
#' @examples
#' data <- simulateMPRA(tr = rep(2,10), da=NULL, nbatch=2, nbc=20)
#' obj <- MpraObject(dnaCounts = data$obs.dna, 
#'                   rnaCounts = data$obs.rna, 
#'                   colAnnot = data$annot)
#' obj <- estimateDepthFactors(obj, lib.factor = "batch", which.lib = "both")
#' ## Upper quantile, using a higher quantile than 0.75:
#' obj <- estimateDepthFactors(obj, lib.factor = "batch", which.lib = "both",
#'                             depth.estimator = function(x) quantile(x, .95, 
#'                                  na.rm=TRUE))
estimateDepthFactors <- function(obj, lib.factor=NULL, which.lib="both", 
                                depth.estimator="uq") {
    
    if(!(which.lib %in% c("dna", "rna", "both"))) {
        stop("which.lib must be 'dna', 'rna' or 'both' (default)")
    }
    if (which.lib %in% c("dna", "both")) {
        dDepth <- estimateFactors(counts=obj@dnaCounts, 
                                  annotations=obj@dnaAnnot, 
                                  lib.factor=lib.factor, 
                                  depth.estimator=depth.estimator)
        if (any(dDepth == 0)) {
            warning("Some DNA library size factors are 0.")
        }
        obj@dnaDepth <- dDepth
    }
    if (which.lib %in% c("rna", "both")) {
        rDepth <- estimateFactors(counts=obj@rnaCounts, 
                                  annotations=obj@rnaAnnot, 
                                  lib.factor=lib.factor, 
                                  depth.estimator=depth.estimator)
        if (any(rDepth == 0)) {
            warning("Some RNA library size factors are 0.")
        }
        obj@rnaDepth <- rDepth
    }
    return(obj)
}

#' estimate the library depth factors
#' @param counts the counts
#' @param annotations the annotations (data.frame)
#' @param lib.factor the name of the factors describing the libraries
#' @param depth.estimator a character indicating which depth estimation to use.
#' Currently supported values are "uq" for upper quartile (default), "totsum"
#' for total sum normalization, and "rle" for RLE (uses geometric mean, and is 
#' therefore not recommended if libraries have 0 counts)
#' @return computed library depth corection factors
#' @noRd
estimateFactors <- function(counts, annotations, lib.factor=NULL, 
                            depth.estimator='uq') {
    if(is.null(lib.factor)) {
        return(rep(1, NCOL(counts)))
    }
    
    if(is.character(lib.factor)){
        if(!all(lib.factor %in% colnames(annotations))) {
            stop("character input for lib.factor must be valid annotations")
        }
        if(length(lib.factor) == 1) {
            lib.factor <- annotations[,lib.factor]
        } else {
            lib.factor <- do.call(paste0, annotations[,lib.factor])
        }
    }
    if(is.function(depth.estimator)) {
        est <- depth.estimator
    } else if (is.character(depth.estimator)) {
        est <- DEPTH_EST_FUNCTIONS[[depth.estimator]]
        if(is.null(est)) {
            stop(depth.estimator, " depth estimation is not supported")
        }
    }
    
    # obj@lib.factor <- as.factor(lib.factor)
    return(compute.depth(data=counts, lib.factor=lib.factor, func=est))
}

#' compute depth factor
#'
#' @param data the data matrix to compute the depth factor for
#' @param lib.factor the partitioning of the data matrix columns into library.
#' A separate depth factor is computed per library.
#' @param func the function to use to compute the library depth
#'
#' @return a numeric of length NCOL(data) with the appropriate epth factor for
#' each column (these are not unique values)
#' @noRd
compute.depth <- function(data, lib.factor, func) {
    lib.factor <- as.factor(lib.factor)

    depth <- by(data = t(data), INDICES = lib.factor, FUN = function(x) {
        func(c(unlist(x)))
    })
    if(depth[1] > 0) {
        depth <- depth / depth[1]
    }

    return(as.vector(depth[lib.factor]))
}

#' Manually set library depth correction factors
#'
#' @param obj the MpraObject
#' @param dnaDepth library size factors for the DNA data, a numeric vector of 
#' length of the number of columns in the DNA data matrix
#' @param rnaDepth library size factors for the RNA data, a numeric vector of 
#' length of the number of columns in the RNA data matrix
#'
#' @return the MpraObject with library depth factors
#' 
#' @export
#' 
#' @examples
#' data <- simulateMPRA(tr = rep(2,10), da=NULL, nbatch=2, nbc=20)
#' obj <- MpraObject(dnaCounts = data$obs.dna, 
#'                   rnaCounts = data$obs.rna, 
#'                   colAnnot = data$annot)
#' ## set constant depth factors (no depth correction)
#' obj <- setDepthFactors(obj, dnaDepth = rep(1, NCOL(data$obs.dna)),
#'                               rnaDepth = rep(1, NCOL(data$obs.rna)))
setDepthFactors <- function(obj, dnaDepth, rnaDepth) {
    if(length(dnaDepth) != NCOL(dnaCounts(obj)) |
       length(rnaDepth) != NCOL(rnaCounts(obj))) {
        stop("invalid length of depth factors")
    }
    obj@dnaDepth <- dnaDepth
    obj@rnaDepth <- rnaDepth

    return(obj)
}
