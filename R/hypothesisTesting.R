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
    df.rna.full <- sapply(obj@modelFits, function(x) x$r.df)
    df.rna.red <- sapply(obj@modelFits.red, function(x) x$r.df)
    df.full <- sapply(obj@modelFits, function(x) x$d.df + x$r.df + x$r.ctrl.df)
    df.red <- sapply(obj@modelFits.red, function(x) x$d.df + x$r.df + x$r.ctrl.df)
    
    lrt <- 2*(ll.full - ll.red)
    df <- df.full-df.red
    pval <- pchisq(lrt, df=df, lower.tail=FALSE)
    fdr <- p.adjust(pval, 'BH')
    
    res <- data.frame(statistic=lrt, pval=pval, fdr=fdr, df.test=df,
                      df.dna=df.dna, df.rna.full=df.rna.full, 
                      df.rna.red=df.rna.red)
    
    ## if condition is single term, extract the corresponding coefficient as logFC
    condition.name <- colnames(obj@designs@rnaFull)[!(colnames(obj@designs@rnaFull) %in% 
                                 colnames(obj@designs@rnaRed))]
    if(length(condition.name) == 1) {
        ## single coefficient is the log Fold Change
        res$logFC <- extractModeParameters.RNA(obj)[,condition.name]
    }
    
    return(res)
}

#' Calculate the significance of a factor in the regression model
#' 
#' @param obj the MpraObject
#' @param factor the name of the factor to make the comparison on
#' @param contrast the character value of the factor to use as a contrast. See details.
#' 
#' @return a data.frame of the results
#' this include the test statistic, logFC, p-value and BH-corrected FDR.
test.coefficient <- function(obj, factor, contrast) {
    if(!(obj@mode == "comparative.coef")) {
        stop("function `analyze.comparative.coef` must be run first for this functionality to be available")
    } else if(!(factor %in% colnames(obj@colAnnot))) {
        stop("given factor: ", factor, " is not included in object annotations")
    } 
    ref <- levels(as.factor(obj@colAnnot[,factor]))[1]
    if (ref == contrast) {
        stop("given contrast ", contrast, " is the reference level of factor ", factor)
    }
    coef.id <- colnames(obj@designs@rnaFull) %in% paste0(factor, contrast)
    if(!any(coef.id)) {
        stop("no matching coefficient for given arguments")
    }
    
    message("comparing ", contrast, " to ", ref, " in factor ", factor, "...")
    coef.id <- 1 + which(coef.id)
    valids <- vapply(obj@modelFits, function(x) !is.null(x$r.se), TRUE)
    logFC <- se <- statistic <- pval <- fdr <- rep(NA, length(obj@modelFits))
    
    logFC[valids] <- vapply(obj@modelFits[valids], function(x) x$r.coef[coef.id], 0.0)
    se[valids] <- vapply(obj@modelFits[valids], function(x) x$r.se[coef.id], 0.0)
    statistic <- (logFC / se) ^ 2
    pval <- pchisq(q = statistic, df = 1, lower.tail = FALSE)
    fdr <- p.adjust(pval, 'BH')
    
    return(data.frame(logFC=logFC, statistic=statistic, pval=pval, fdr=fdr, row.names = names(obj@modelFits)))
}