#' Calculate likelihood ratio test for the specific nested model
#' @import stats
#'
#' @note Must be run after running an LRT-based analysis
#'
#' @param obj the MpraObject containing the full and reduced
#' 
#' @export 
#' @return results data frame
test.lrt <- function(obj) {
    if(length(obj@modelFits) == 0 | length(obj@modelFits.red) == 0) {
        stop("An LRT analysis must be performed before computing the test")
    }
    
    message("Performing Likelihood Ratio Test...")
    
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
        res$logFC <- extractModelParameters.RNA(obj)[,condition.name]
    }
    
    return(res)
}

#' Calculate the significance of a factor in the regression model
#' 
#' @import stats
#' 
#' @param obj the MpraObject
#' @param factor the name of the factor to make the comparison on
#' @param contrast the character value of the factor to use as a contrast. See 
#' details.
#' 
#' @export
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
    
    message("Testing for significance:  ", contrast, " vs. ", ref, ", in factor ", factor, "...")
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

#' test for significant activity (quantitative analysis) using various empirical
#' tests (see details)
#' 
#' @param obj the MpraObject, after running an analysis function
#' @param statistic if null [default], the intercept term is used as the score.
#' An alternate score can be provided by setting 'statistic'. Must be a numeric
#' vector.
#' 
#' @export
#' @return a data.frame of empirical summary statistics based on the model's 
#' estimate of slope, or the given statistic. These are:
#' \itemize{
#'     \item statistic: the statistic (either the provided, or extracted from the
#'     models)
#'     \item zscore: Z-score of the statistic (number of standard devisations 
#'     from the mean)
#'     \item mad.score: a median-baed equivalent of the Z-score, with less 
#'     sensitivity to outlier values
#'     \item zscore.ctrl: only available if negative controls are provided.
#'     a Z-score based on the controls distribution, instead of the distribution 
#'     of the complete set of observations
#'     \item mad.score.ctrl: only available if negative controls are provided.
#'     a MAD-score based on the controls distribution, instead of the distribution 
#'     of the complete set of observations
#'     \item epval: only available if negative controls are provided. empirical 
#'     P-value, using the control distribution as the null
#'     \item fdr: only available if negative controls are provided. BH adjusted-
#'     empricial p-values
#' }
test.empirical <- function(obj, statistic=NULL) {
    
    if(is.null(statistic)) {
        ## extract slope - second coefficient of the rna model
        statistic <- vapply(obj@modelFits, function(x) x$r.coef[2], 0.0)
        names(statistic) <- names(obj@modelFits)
    }
    
    zscore <- (statistic - mean(statistic, na.rm=TRUE)) / sd(statistic, 
                                                             na.rm=TRUE)
    mad.score <- (statistic - median(statistic, na.rm=TRUE)) / mad(statistic, 
                                                                   na.rm=TRUE)
    
    res <- data.frame(statistic=statistic,
                      zscore=zscore,
                      mad.score=mad.score)
    
    if(!is.null(obj@controls)) {
        ctrls <- statistic[obj@controls]
        res$zscore.ctrl <- ((statistic - mean(ctrls, na.rm=TRUE)) / 
                                sd(ctrls, na.rm=TRUE))
        res$mad.score.ctrl <- ((statistic - median(ctrls, na.rm=TRUE)) / 
                                   mad(ctrls, na.rm=TRUE))
        
        res$epval <- 1 - ecdf(ctrls)(statistic)
        res$fdr <- p.adjust(res$epval, "BH")
    }
    
    return(res)
}
    