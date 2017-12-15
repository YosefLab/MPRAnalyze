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
    ##TODO: format these as a data.frame to begin with?
    ll.full <- sapply(obj@modelFits, function(x) x$ll)
    ll.red <- sapply(obj@modelFits.red, function(x) x$ll)
    df.dna <- sapply(obj@modelFits, function(x) x$d.df)
    df.full.rna <- sapply(obj@modelFits, function(x) x$r.df)
    df.red.rna <- sapply(obj@modelFits.red, function(x) x$r.df)
    df.full.rna.ctrl <- sapply(obj@modelFits, function(x) x$r.ctrl.df)
    df.red.rna.ctrl <- sapply(obj@modelFits.red, function(x) x$r.ctrl.df)
    df.full <- sapply(obj@modelFits, function(x) x$d.df + x$r.df + x$r.ctrl.df)
    df.red <- sapply(obj@modelFits.red, function(x) x$d.df + x$r.df + x$r.ctrl.df)
    
    lrt <- 2*(ll.full - ll.red)
    pval <- pchisq(lrt, df = df.full-df.red, lower.tail=FALSE)
    fdr <- p.adjust(pval, 'BH')
    
    return(data.frame(statistic=lrt, pval=pval, fdr=fdr))
}

test.ttest <- function(obj, condition) {
    ## check if there is an intercept or not
    
    ## mode: comparative.ttest.1ref (has intercept, coeff = offset)
    ## mode: comparative.ttest.mrefs (no intercept, all pairs possible)
    
    ## extract slopes, reformat for results data.frame
    ## TODO: test.ttest:
    ##    if .1ref: get one argmuent, compute ttest
    ##    if .mrefs: get two arguments, compute ttest
    ##    return results and append them to the results data.frame
}