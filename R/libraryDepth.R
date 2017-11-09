#' estimate library size correction factors
#'
#'@param obj an mpraObject
estimateDepthFactors <- function(obj, batchids=NULL) {
    ## library size already set
    if(!is.null(obj@dnaDepth)) {
        return(obj)
    } else if(is.null(batchids)) {
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
