#' Calculate likelihood ratio test for the specific nested model
#'
#' @note Must be run after running fit.differential
#' TODO: adatapt this to only perform test and not fit reduced model
#'
#' @param obj the MpraObject containing the full and reduced
#' @param dnaDesign.reduced the design of the DNA reduced model. If NULL, DNA is
#' assumed to follow same model as the full model - and raDesign.reduced must be
#' provided 
#' @param rnaDesign.reduced the design of the RNA reduced model. If NULL, RNA is
#' assumed to follow same model as the full model, and dnaDesign.reduced must be
#' provided
#' 
#' @return results data frame
test.lrt <- function(obj) {
    
    ll.full <- sapply(obj@modelFits, function(x) x$ll)
    ll.red <- sapply(obj@modelFits.red, function(x) x$ll)
    df.dna <- sapply(obj@modelFits, function(x) x$d.df )
    df.full.rna <- sapply(obj@modelFits, function(x) x$r.df )
    df.red.rna <- sapply(obj@modelFits.red, function(x) x$r.df)
    df.full.rna.ctrl <- sapply(obj@modelFits, function(x) x$r.ctrl.df )
    df.red.rna.ctrl <- sapply(obj@modelFits.red, function(x) x$r.ctrl.df )
    df.full <- sapply(obj@modelFits, function(x) x$d.df + x$r.df + x$r.ctrl.df )
    df.red <- sapply(obj@modelFits.red, function(x) x$d.df + x$r.df + x$r.ctrl.df )
    #df.full <- NCOL(obj@designs@dnaFull) + NCOL(obj@designs@rnaFull)
    #df.red <- NCOL(obj@designs@dnaRed) + NCOL(obj@designs@rnaRed)
    res <- data.frame(enhancer=rownames(obj@dnaCounts),
                      ll.full=ll.full,
                      ll.red=ll.red,
                      df.full=df.full,
                      df.red=df.red,
                      df.dna=df.dna,
                      df.full.rna=df.full.rna,
                      df.red.rna=df.red.rna,
                      df.full.rna.ctrl=df.full.rna.ctrl,
                      df.red.rna.ctrl=df.red.rna.ctrl,
                      row.names=rownames(obj@dnaCounts))
    
    res$lrt.statistic <- 2*(res$ll.full - res$ll.red)
    res$pval <- pchisq(res$lrt.statistic, df = res$df.full-res$df.red, 
                       lower.tail=FALSE)
    res$fdr.pval <- p.adjust(res$pval, 'BH')
    
    return(res)
}
