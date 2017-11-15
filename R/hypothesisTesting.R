
#' Calculate likelihood ratio test for the specific nested model
#'
#' @note Must be run after running fit.differential
#' TODO: adatapt this to only perform test and not fit reduced model
#'
#' @param obj the MpraObject containing the full model
#' @param dnaDesign.reduced the design of the DNA reduced model. If NULL, DNA is
#' assumed to follow same model as the full model - and raDesign.reduced must be
#' provided 
#' @param rnaDesign.reduced the design of the RNA reduced model. If NULL, RNA is
#' assumed to follow same model as the full model, and dnaDesign.reduced must be
#' provided
test.lrt <- function(obj, dnaDesign.reduced=NULL, rnaDesign.reduced) {
    ## TODO: check that mode is 'LRT', or just that it's not 'quant'?
    if (obj@model != "LRT") {
        stop("fit.differential must be run in 'LRT' mode prior to calling this function")
    } else if (is.null(dnaDesign.reduced) & is.null(rnaDesign.reduced)) {
        stop("DNA and RNA can't both be the same as the full model")
    } else if (is.null(dnaDesign.reduced)) {
        dnaDesign.reduced <- obj@designs@dnaFull
    } else if (is.null(rnaDesign.reduced)) {
        rnaDesign.reduced <- obj@designs@rnaFull
    }

    obj@designs@dnaRed <- getDesignMat(obj, dnaDesign.reduced)
    obj@designs@rnaRed <- getDesignMat(obj, rnaDesign.reduced)

    ##TODO: make sure the reduced model is nested within the full model?

    obj@ModelFits.reduced <- bplapply(rownames(obj@dnaCounts), function(rn) {
        return(fitfun(dcounts=obj@dnaCounts[rn,],
                      rcounts=obj@rnaCounts[rn,],
                      ddepth=obj@dnaDepth,
                      rdepth=obj@rnaDepth,
                      ddesign.mat=obj@designs@dnaRed,
                      rdesign.mat=obj@designs@rnaRed,
                      compute.hessian=compHess))
    }, BPPARAM = obj@BPPARAM)

    ll.full <- unlist(lapply(obj@modelFits,
                             function(x) x$ll), use.names = TRUE)
    ll.red <- unlist(lapply(obj@modelFits.reduced,
                            function(x) x$ll), use.names = TRUE)
    res <- merge(data.frame('full'=ll.full),
                 data.frame('reduced'=ll.red), by = 0, all=FALSE)
    colnames(res)[1] <- "Feature" ##TODO: what do we call these?

    res$lrt.statistic <- 2*(res$full - res$reduced)

    df.full <- NCOL(obj@designs@dnaFull) + NCOL(obj@designs@rnaFull)
    df.red <- NCOL(obj@designs@dnaRed) + NCOL(obj@designs@rnaRed)
    res$p.val <- pchisq(res$lrt.statistic, df = df.full-df.red)
    res$fdr <- p.adjust(res$p.val, 'BH')

    obj@hyptestResults <- res
    return(obj)
}
