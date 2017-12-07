
#' estimate library size correction factors
#'
#' @note since in most MPRA experiments multiple barcodes exist within a single
#' library, each column in the matrix is usually not a separate library. For this
#' reason, it is recommended to supply this function with the appropriate
#' partitioning of the data matrix columns into libraries, see lib.factor
#'
#' @param obj the MpraObject
#' @param lib.factor the factor associating each sample to a library. Can be a
#' factor or the name of a column in the object's colAnnot. If not provided, the
#' data is assumed to have been generated from a single library, and constant
#' library depth is set.
#' @param depthEstimator a character indicating which depth estimation to use.
#' Currently supported values are "uq" for upper quantile (default) and "rle"
#' for RLE (uses geometric mean, and is therefore not recommended if libraries
#' have 0 counts)
#'
#' @return the MpraObject with estimated values for sequencing depth factors
#' 
#' @export
estimateDepthFactors <- function(obj, lib.factor=NULL, depthEstimator='uq') {
    ## library size already set
    if(length(obj@dnaDepth) > 0) {
        return(obj)
    }

    if(is.null(lib.factor)) {
        ##library factor not provided, assume single library (=neutral depth)
        obj@lib.factor <- as.factor(rep(1, NCOL(obj@dnaCounts)))
        obj@dnaDepth <- rep(1, NCOL(obj@dnaCounts))
        obj@rnaDepth <- rep(1, NCOL(obj@rnaCounts))
        return(obj)
    }
    
    if(is.character(lib.factor)){
        if(!(lib.factor %in% colnames(obj@colAnnot))) {
            stop("character input for lib.factor must be name of a column in object's colAnnot")
        }
        lib.factor <- obj@colAnnot[,lib.factor]
    }

    est <- DEPTH_EST_FUNCTIONS[[depthEstimator]]
    if(is.null(est)) {
        stop(depthEstimator, " depth estimation is not supported")
    }
    
    obj@lib.factor <- as.factor(lib.factor)
    obj@dnaDepth <- compute.depth(data=obj@dnaCounts,
                                  lib.factor=lib.factor,
                                  func=est)
    obj@rnaDepth <- compute.depth(data=obj@rnaCounts,
                                  lib.factor=lib.factor,
                                  func=est)
    return(obj)
}

#' compute depth factor
#'
#' @param data the data matrix to compute the depth factor for
#' @param lib.factor the partitioning of the data matrix columns into library.
#' A separate depth factor is computer per library.
#'
#' @return a numeric of length NCOL(data) with the appropriate epth factor for
#' each column (these are not unique values)
compute.depth <- function(data, lib.factor, func) {
    lib.factor <- as.factor(lib.factor)

    depth <- by(data = t(data), INDICES = lib.factor, FUN = func)
    if(median(depth) > 0) {
        depth <- depth / median(depth)
    }

    return(as.vector(depth[lib.factor]))
}

#' easy access container for depth estimation functions
DEPTH_EST_FUNCTIONS = list(
    uq = function(x) quantile(c(unlist((x))), .75),
    rle = function(x) exp(mean(log(c(unlist((x))))))
)

#' Manually set library depth correction factors
#'
#' @param obj the MpraObject
#' @param dnaDepth library size factors for the DNA data, a numeric vector of length
#' of the number of columns in the DNA data matrix
#' @param rnaDepth library size factors for the RNA data, a numeric vector of length
#' of the number of columns in the RNA data matrix
#'
#' @return the MpraObject with library depth factors
#' 
#' @export
setDepthFactors <- function(obj, dnaDepth, rnaDepth) {
    if(length(dnaDepth) != NCOL(obj@dnaCounts) 
       | length(rnaDepth) != NCOL(obj@rnaCounts)) {
        stop("invalid length of depth factors")
    }
    obj@dnaDepth <- dnaDepth
    obj@rnaDepth <- rnaDepth

    return(obj)
}
