#' Calculate likelihood ratio test for the specific nested model
#' @import stats
#'
#' @note Must be run after running an LRT-based analysis
#'
#' @param obj the MpraObject containing the full and reduced
#' 
#' @export 
#' @return results data frame
#' 
#' @examples
#' data <- simulateMPRA(tr = rep(2,5), da=c(rep(0,2), rep(1,3)), 
#'                      nbatch=2, nbc=15)
#' obj <- MpraObject(dnaCounts = data$obs.dna, 
#'                   rnaCounts = data$obs.rna, 
#'                   colAnnot = data$annot)
#' obj <- estimateDepthFactors(obj, lib.factor = "batch", which.lib = "both")
#' obj <- analyzeComparative(obj, dnaDesign = ~ batch + barcode + condition, 
#'                               rnaDesign = ~ condition, reducedDesign = ~ 1)
#' results <- testLrt(obj)
testLrt <- function(obj) {
    if ("mode" %in% slotNames(obj) & (obj@mode == "scale")) {
        return(testLrt.scale(obj))
    }
    if(length(obj@modelFits) == 0 | length(obj@modelFits.red) == 0) {
        stop("An LRT analysis must be performed before computing the test")
    }
    
    message("Performing Likelihood Ratio Test...")
    
    ll.full <- obj@modelFits$ll
    ll.red <- obj@modelFits.red$ll
    df.dna <- obj@modelFits$d.df
    df.rna.full <- obj@modelFits$r.df + obj@modelFits$r.ctrl.df + 
        length(obj@rnaCtrlScale)
    df.rna.red <- obj@modelFits.red$r.df + obj@modelFits.red$r.ctrl.df + 
        length(obj@rnaCtrlScale)
    df.full <- df.dna + df.rna.full 
    df.red <- df.dna + df.rna.red 
    
    lrt <- 2*(ll.full - ll.red)
    df <- df.full-df.red
    pval <- pchisq(lrt, df=df, lower.tail=FALSE)
    fdr <- p.adjust(pval, 'BH')
    
    res <- data.frame(statistic=lrt, pval=pval, fdr=fdr, df.test=df,
                    df.dna=df.dna, df.rna.full=df.rna.full, 
                    df.rna.red=df.rna.red)
    ## if condition is single term, extract the coefficient as logFC
    condition.name <- (colnames(obj@designs@rnaFull)
                       [!(colnames(obj@designs@rnaFull) %in% 
                              colnames(obj@designs@rnaRed))])
    if(length(condition.name) == 1) {
        ## single coefficient is the log Fold Change
        res$logFC <- getModelParameters_RNA(obj)[,condition.name]
    }
    
    return(res)
}

#' @rdname testLrt
#' @noRd
testLrt.scale <- function(obj) {
    if(length(obj@modelFits) == 0 | length(obj@modelFits.red) == 0) {
        stop("An LRT analysis must be performed before computing the test")
    }
    
    message("Performing Likelihood Ratio Test...")
    
    ll.full <- obj@modelFits$ll
    ll.red <- obj@modelFits.red$ll
    df.rna.full <- obj@modelFits$r.df + obj@modelFits$r.ctrl.df + 
        length(obj@rnaCtrlScale)
    df.rna.red <- obj@modelFits.red$r.df + obj@modelFits.red$r.ctrl.df + 
        length(obj@rnaCtrlScale)
    df.full <- df.rna.full 
    df.red <- df.rna.red 
    
    lrt <- 2*(ll.full - ll.red)
    df <- df.full-df.red
    pval <- pchisq(lrt, df=df, lower.tail=FALSE)
    fdr <- p.adjust(pval, 'BH')
    
    res <- data.frame(statistic=lrt, pval=pval, fdr=fdr, df.test=df,
                      df.rna.full=df.rna.full, 
                      df.rna.red=df.rna.red)
    ## if condition is single term, extract the coefficient as logFC
    condition.name <- (colnames(obj@designs@rnaFull)
                       [!(colnames(obj@designs@rnaFull) %in% 
                              colnames(obj@designs@rnaRed))])
    if(length(condition.name) == 1) {
        ## single coefficient is the log Fold Change
        res$logFC <- getModelParameters_RNA(obj)[,condition.name]
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
#' @examples
#' data <- simulateMPRA(tr = rep(2,5), da=c(rep(0,2), rep(1,3)), 
#'                      nbatch=2, nbc=15)
#' obj <- MpraObject(dnaCounts = data$obs.dna, 
#'                   rnaCounts = data$obs.rna, 
#'                   colAnnot = data$annot)
#' obj <- estimateDepthFactors(obj, lib.factor = "batch", which.lib = "both")
#' 
#' ## fit.se must be TRUE for coefficient based testing to work
#' obj <- analyzeComparative(obj, dnaDesign = ~ batch + barcode + condition, 
#'                               rnaDesign = ~ condition, fit.se = TRUE)
#' results <- testCoefficient(obj, "condition", "contrast")
testCoefficient <- function(obj, factor, contrast) {
    if(is.null(obj@modelFits$r.se)) {
        stop("Model fitting did not include standard error estimation.\
             Coefficient-based testing cannot be perfromed.")
    }
    
    if(!(factor %in% colnames(rnaAnnot(obj)))) {
        stop("given factor: ", factor, 
            " is not included in object annotations")
    } 
    ref <- levels(as.factor(rnaAnnot(obj)[,factor]))[1]
    if (ref == contrast) {
        stop("given contrast ", contrast, 
            " is the reference level of factor ", factor)
    }
    coef.id <- colnames(obj@designs@rnaFull) %in% paste0(factor, contrast)
    if(!any(coef.id)) {
        stop("no matching coefficient for given arguments")
    }
    
    message("Testing for significance:  ", contrast, " vs. ", ref, 
            ", in factor ", factor, "...")
    coef.id <- 1 + which(coef.id)
    
    # valids <- !is.null(obj@modelFits$r.se)
    valids <- apply(obj@modelFits$r.se, 1, function(x) !all(is.na(x)))
    logFC <- se <- statistic <- pval <- fdr <- rep(NA, NROW(dnaCounts(obj)))
    logFC[valids] <- obj@modelFits$r.coef[valids,coef.id]
    se[valids] <- obj@modelFits$r.se[valids,coef.id]
    statistic <- (logFC / se) ^ 2
    pval <- pchisq(q = statistic, df = 1, lower.tail = FALSE)
    fdr <- p.adjust(pval, 'BH')
    
    res <- data.frame(logFC=logFC, statistic=statistic, 
                      pval=pval, fdr=fdr, 
                      row.names = rownames(dnaCounts(obj)))
    return(res)
}

#' test for significant activity (quantitative analysis) using various empirical
#' tests (see details)
#' 
#' @param obj the MpraObject, after running an analysis function
#' @param statistic if null [default], the intercept term is used as the score.
#' An alternate score can be provided by setting 'statistic'. Must be a numeric
#' vector.
#' @param useControls is TRUE and controls are available, use the controls to
#' establish the background model and compare against. This allows for more
#' accurate zscores as well as empircal p-values.
#' @param twoSided should the p-value be from a two-sided test (default: FALSE,
#' right-side test)
#' @param subset only test a subset of the enhancers in the object (logical,
#' indices or names). Default is NULL, then all the enhancers are included.
#' @export
#' @return a data.frame of empirical summary statistics based on the model's 
#' estimate of slope, or the given statistic. These are:
#' \itemize{
#'     \item statistic: the statistic (either the provided, or extracted from 
#'     the models)
#'     \item zscore: Z-score of the statistic (number of standard devisations 
#'     from the mean). If controls are available, the score is based on their 
#'     distribution: so it's the number of control-sd from the control-mean
#'     \item mad.score: a median-baed equivalent of the Z-score, with less 
#'     sensitivity to outlier values. If controls are provided, it's based
#'     on their distribution.
#'     \item pval.zscore: a p-value based on the normal approximation of the
#'     Z-scores
#'     \item pval.empirical: only available if negative controls are provided. 
#'     empirical P-value, using the control distribution as the null
#' }
#' 
#' @examples
#' data <- simulateMPRA(tr = rep(2,10), da=NULL, nbatch=2, nbc=15)
#' obj <- MpraObject(dnaCounts = data$obs.dna, 
#'                   rnaCounts = data$obs.rna, 
#'                   colAnnot = data$annot)
#' obj <- estimateDepthFactors(obj, lib.factor = "batch", which.lib = "both")
#' obj <- analyzeQuantification(obj, dnaDesign = ~ batch + barcode, 
#'                               rnaDesign = ~1)
#' results <- testEmpirical(obj)
#' 
#' ## or test with a different statistic:
#' aggregated.ratio <- rowSums(data$obs.rna) / rowSums(data$obs.dna)
#' results <- testEmpirical(obj, aggregated.ratio)
testEmpirical <- function(obj, statistic=NULL, useControls=TRUE, twoSided=FALSE,
                          subset=NULL) {
    
    if(is.null(statistic)) {
        alpha <- getAlpha(obj)
        statistic <- alpha$alpha
        names(statistic) <- rownames(alpha)
    }
    
    if(!is.null(subset)) {
        if(is.character(subset)) {
            subset <- rownames(dnaCounts(obj)) %in% subset
        }
        if(is.logical(subset)) {
            subset <- which(subset)
        }
        statistic <- statistic[subset]
    }
    
    res <- data.frame(statistic=statistic)
    
    if(!any(controls(obj)) | !useControls) {
        ## No controls, use bottom of the distribution to establish the baseline
        
        #estimate mode of distribution
        dist.peak <- est.mode(statistic)
        base <- statistic[statistic < dist.peak]
        
        std.div <- sqrt(mean((base - dist.peak) ^ 2, na.rm = TRUE))
        mad.score <- 1.4826 * median(abs(base - dist.peak), na.rm = TRUE)
        
        res$zscore <- (statistic - dist.peak) / std.div
        res$mad.score <- (statistic - dist.peak) / mad.score
        
    } else {
        ctrl.idx <- controls(obj)
        if (!is.null(subset)) {
            ctrl.idx <- ctrl.idx[subset]
        }
        ctrls <- statistic[ctrl.idx]
        
        res$control <- FALSE
        res$control[ctrl.idx] <- TRUE
        
        res$zscore <- ((statistic - mean(ctrls, na.rm=TRUE)) / 
                                sd(ctrls, na.rm=TRUE))
        res$mad.score <- ((statistic - median(ctrls, na.rm=TRUE)) / 
                                mad(ctrls, na.rm=TRUE))
        res$pval.empirical <- 1 - ecdf(ctrls)(statistic)
    }
    
    if (twoSided) {
        res$pval.mad <- 2 * pnorm(abs(res$mad.score), lower.tail = FALSE)
        res$pval.zscore <- 2 * pnorm(abs(res$zscore), lower.tail = FALSE)
    } else {
        res$pval.mad <- pnorm(res$mad.score, lower.tail = FALSE)
        res$pval.zscore <- pnorm(res$zscore, lower.tail = FALSE)
    }
    
    return(res)
}

#' estimate mode of a sample
#' @param x the sample to estimate the mode of
#' @return the estimate mode of the sample distribution
#' @details the etimated mode is the center of the window that contains the most
#' observations. Window size is set to 2% of the range of values.
#' @noRd
est.mode <- function(x) {
    cdf <- ecdf(x)
    win.size <- (max(x) - min(x)) * 0.01
    f <- cdf(x + win.size) - cdf(x - win.size)
    m <- x[which.max(f)]
    return(m)
}