#' estimate library size correction factors
#'
#'@param obj an mpraObject
estimateDepthFactors <- function(obj, batchids=NULL) {
    if(is.null(batchids)) {
        ##TODO: now what?
    } else {
        ## TODO: compute size factor per batch
        ## TODO: redistribute factors to the appropriate vector
    }
}

setDepthFactors <- function(obj, dnaDepth, rnaDepth) {
    ##TODO: validate lengths, etc
    obj@dnaDepth <- dnaDepth
    obj@rnaDepth <- rnaDepth
}
