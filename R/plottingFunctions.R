#' Output observed ratios as boxplots and add fits
#' 
#' Aid intuition of model fits by visualisation of the fit 
#' to one enhancer.
#' 
#' @param obj (MpraObject)
#' MpraObject from which enhancer is to be visualised.
#' @param id (string)
#' ID of enhancer which is to be visualised. The control set is 
#' taken from obj if controls were used.
#' @param condition (string)
#' Column of colAnnot that contains condition which should be resolved on
#' x-axis. Typically, this would be the condition that was tested for 
#' via analyse.condition.*()
#' @param batch (string)
#' Column of colAnnot that contains batch allocation of confounding variable
#' that was corrected for during fitting.
#' @param full (bool)
#' Whether to plot full or reduced model fits.
#' @param show.outliers (bool)
#' Whether to show outliers in boxplot. This slows if many samples were observed.
#' 
#' @return (ggplot) Boxplot graphic.
#' 
#' @import ggplot2
#' @export
plotBoxplots <- function(obj, id, condition=NULL, batch=NULL, full=TRUE, 
                         show.outliers=FALSE){
    
    ## set outlier handling
    outlier.shape <- if(show.outliers) "x" else NA

    ## check whether test results available
    
    if(!is.null(obj@results)) {
        fdr <- round(log(obj@results[id,]$fdr)/log(10),2)
    } else {
        fdr <- NA
    }
    
    ## extract observations and format into data-frame of rna:dna ratios
    # case enhancer observations
    gplot.data.boxplot <- data.frame(
        ratio=(log(obj@rnaCounts[id,])-log(obj@rnaDepth))-
            (log(obj@dnaCounts[id,])-log(obj@dnaDepth)),
        cond=obj@dnaAnnot[,condition],
        batch=obj@dnaAnnot[,batch],
        enhancer=id,
        type="obs",
        stringsAsFactors=FALSE
    )
    # control enhancer observations
    if(!is.null(obj@controls)){
        gplot.data.obs.ctrl <- do.call(rbind, lapply(obj@controls, function(i) {
            data.frame(
                ratio=(log(obj@rnaCounts[i,])-log(obj@rnaDepth))-
                    (log(obj@dnaCounts[i,])-log(obj@dnaDepth)),
                cond=obj@dnaAnnot[,condition],
                batch=obj@dnaAnnot[,batch],
                enhancer="controls",
                type="obs",
                stringsAsFactors=FALSE
            )
        }))
        gplot.data.boxplot <- rbind(gplot.data.boxplot, gplot.data.obs.ctrl)
    }
    
    ## generate sample from fitted distribution
    # extract model fits for each observation
    obs.resampled <- resampleObs(obj, enhancer=id, full=full)
    gplot.data.fit <- data.frame(
        ratio=log(obs.resampled$rna)-log(obs.resampled$dna),
        cond=obj@dnaAnnot[,condition],
        batch=obj@dnaAnnot[,batch],
        enhancer=id,
        type="fit",
        stringsAsFactors=FALSE
    )
    gplot.data.fit <- gplot.data.fit[!duplicated(gplot.data.fit),]
    if(!is.null(obj@controls)){
        gplot.data.fit.ctrl <- do.call(rbind, lapply(obj@controls, function(i) {
            obs.resampled <- resampleObs(obj, enhancer=i, full=full)
            gplot.data.fit.ctrl.i <- data.frame(
                ratio=log(obs.resampled$rna)-log(obs.resampled$dna),
                cond=obj@dnaAnnot[,condition],
                batch=obj@dnaAnnot[,batch],
                enhancer="controls",
                type="fit",
                stringsAsFactors=FALSE
            )
            return(gplot.data.fit.ctrl.i)
        }))
        gplot.data.fit.ctrl <- gplot.data.fit.ctrl[!duplicated(gplot.data.fit.ctrl),]
        gplot.data.fit <- rbind(gplot.data.fit, gplot.data.fit.ctrl)
    }
    gplot.data.boxplot <- rbind(gplot.data.boxplot, gplot.data.fit)
    
    ## create ggplot
    # colour-blind palette
    cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                    "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    if(!is.null(condition)) {
        if(!is.null(batch)) {
            gplot.boxplot <- ggplot() + geom_boxplot(
                data=gplot.data.boxplot, aes_(
                    x = ~cond, y = ~ratio, fill = ~enhancer, 
                    alpha = ~batch, linetype = ~type),
                outlier.shape = outlier.shape) +
                scale_shape_discrete(solid = FALSE) +
                labs(title=paste0(id, " log10 fdr-corrected p-value: ", fdr))
        } else {
            gplot.boxplot <- ggplot() + geom_boxplot(
                data=gplot.data.boxplot, aes_(
                    x = ~cond, y = ~ratio, fill = ~enhancer, linetype = ~type),
                outlier.shape = outlier.shape) +
                scale_shape_discrete(solid = FALSE) +
                labs(title=paste0(id, " log10 fdr-corrected p-value: ", fdr))
        }
    } else {
        if(!is.null(batch)) {
            gplot.boxplot <- ggplot() + geom_boxplot(
                data=gplot.data.boxplot, aes_(
                    x = ~enhancer, y = ~ratio, fill = ~enhancer, alpha = ~batch, 
                    linetype = ~type),
                outlier.shape = outlier.shape) +
                scale_shape_discrete(solid = FALSE) +
                labs(title=paste0(id, " log10 fdr-corrected p-value: ", fdr))
        } else {
            gplot.boxplot <- ggplot() + geom_boxplot(
                data=gplot.data.boxplot, aes_(
                    x = ~enhancer, y = ~ratio, 
                    fill = ~enhancer, linetype = ~type),
                outlier.shape = outlier.shape) +
                scale_shape_discrete(solid = FALSE) +
                labs(title=paste0(id, " log10 fdr-corrected p-value: ", fdr))
        }
    }
    if(length(unique(gplot.data.boxplot$batch)) <= length(cbbPalette)) {
        gplot.boxplot <- gplot.boxplot + 
            scale_fill_manual(values = cbbPalette) + 
            scale_colour_manual(values = cbbPalette)
    }
    
    return(gplot.boxplot)
}

#' create a scatter plot of the model alpha (trascription rate) vs. the naive
#' estimator (RNA / DNA).
#' @param obj the MpraObject
#' @param condition if NULL (default) create a single plot of the model baseline
#' (intercept) vs the total data. Else, condition must be a valid factor in the 
#' object column annotations that was included in the design. In that case, a 
#' plot is generated for each level of the factor
#' @param log plot in log scale
#' @export
plotAlphaRatio <- function(obj, condition = NULL, log=FALSE, categories=NULL) {
    if(is.null(condition)) {
        ratio <- rowMeans(obj@rnaCounts) / rowMeans(obj@dnaCounts)
        alpha <- getAlpha(obj)
        res <- ggplot(data = data.frame(ratio = ratio, 
                                        alpha = alpha, 
                                        category=categories)) + 
            geom_point(mapping = aes(x = ratio, y = alpha, color=category)) + 
            geom_abline(intercept=0, slope=1) + 
            xlab("RNA / DNA") + ylab("alpha")
        if(log) {
            res <- res + scale_x_log10() + scale_y_log10() +
                xlab("log(RNA / DNA)") + ylab("log(alpha)")
        }
    } else {
        first <- TRUE
        res <- lapply(levels(as.factor(obj@dnaAnnot[,condition])), function(l) {
            if(first) {
                alpha <- getAlpha(obj)
                first <- FALSE
            } else {
                alpha <- getAlpha(obj, condition, l)
            }
            if(condition %in% colnames(obj@dnaAnnot)) {
                idx.dna <- obj@dnaAnnot[,condition] == l
            } else {
                idx.dna <- 1:NCOL(obj@dnaCounts)
            }
            
            idx.rna <- obj@rnaAnnot[,condition] == l
            ratio <- (rowMeans(obj@rnaCounts[,idx.rna,drop=FALSE]) / 
                     rowMeans(obj@dnaCounts[,idx.dna,drop=FALSE]))
            
            res <- ggplot(data = data.frame(ratio = ratio, 
                                            alpha = alpha,
                                            category = categories)) + 
                geom_point(mapping = aes(x = ratio, y = alpha, color=categories)) + 
                geom_abline(intercept=0, slope=1) + 
                xlab("RNA / DNA") + ylab("alpha") + ggtitle(paste(condition, l))
            if(log) {
                res <- res + scale_x_log10() + scale_y_log10() +
                    xlab("log(RNA / DNA)") + ylab("log(alpha)")
            }
            return(res)
        })
    }
    return(res)
}


#' plot the CDF of the pvalues
#' @param p the p-values to plot
#' @param categories a factor (same length as p) dividing the p-values into 
#' categories to be plotted separately
#' @param adjusted if TRUE, adjust the p-values using the BH FDR method.
plotPvalCDF <- function(p, categories, adjusted=FALSE) {
    if(adjusted) {
        p <- p.adjust(p, "BH")
    }
    df <- data.frame(p=p, cat=categories)
    ggplot(df, aes_(~p, group=~cat, color=~cat)) + stat_ecdf(geom="step") + 
        xlab("P-value") + ylab("CDF")
    
}

#' Plot the observed and expected distributions of a given enhancer
#' @param obj MpraObject after model fitting
#' @param enhancer the id of the enhancer to plot (index or name)
#' @param rna if TRUE, plot the RNA distribution. Otherwise plot the DNA.
#' @import ggplot2
#' @export
plotObsExpDistributions <- function(obj, enhancer, bins=NULL) {
    if(length(enhancer) > 1) {
        stop("please supply a single enhancer (index or name)")
    }
    df <- data.frame(obs.r = obj@rnaCounts[enhancer,],
                     exd.r = as.numeric(getFits.RNA(obj, enhancer)),
                     obs.d = obj@dnaCounts[enhancer,], 
                     exd.d = as.numeric(getFits.DNA(obj, enhancer)))
    df$obs.r[df$obs.r <= 0] <- NA
    df$obs.d[df$obs.d <= 0] <- NA
    
    dp <- ggplot(df) + 
        geom_histogram(aes_(x = ~obs.d, ~..density.., fill="Observed")) + 
        geom_density(aes_(x= ~exd.d, color="Expected"), size=2) + 
        scale_fill_manual(name=element_blank(), values = c("Observed"='grey33')) + 
        scale_colour_manual(name=element_blank(), values = c('Expected'='black')) + 
        theme(text=element_text(size=20), legend.position = "none") +
        xlab("counts")
    rp <- ggplot(df) + 
        geom_histogram(aes_(x = ~obs.r, ~..density.., fill="Observed")) + 
        geom_density(aes_(x= ~exd.r, color="Expected"), size=2) + 
        scale_fill_manual(name=element_blank(), values = c("Observed"='grey33')) + 
        scale_colour_manual(name=element_blank(), values = c('Expected'='black')) + 
        theme(text=element_text(size=20), legend.position = "none") +
        xlab("counts")
    return(list(dna=dp, rna=rp))
}

#' plot the relationship between the mean and the variance of the DNA and RNA 
#' distributions. Plots are in log scale.
#' @param obj the MpraObject to plot
#' @return a list of two ggplot objects: 'dna' and 'rna'.
#' @import ggplot2
#' @export
plotMeanVariance <- function(obj) {
    df.d <- data.frame(mean = rowMeans(obj@dnaCounts),
                       var = apply(obj@dnaCounts, 1, var))
                     
    df.r <- data.frame(mean = rowMeans(obj@rnaCounts),
                       var = apply(obj@rnaCounts, 1, var))
    
    pd <- ggplot(data = df.d, aes(log(mean), log(var))) + 
        geom_point() +
        theme(text=element_text(size=20), legend.position = "none")
    pr <- ggplot(data = df.r, aes(log(mean), log(var))) + 
        geom_point() + 
        theme(text=element_text(size=20), legend.position = "none")
    
    return(list(dna = pd, rna = pr))
}
