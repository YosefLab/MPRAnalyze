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
#' @return (gplot) Boxplot graphic.
#' 
#' @import ggplot2
#' 
#' @export
plotBoxplots <- function(obj, id, condition=NULL, batch=NULL, full=TRUE, show.outliers=FALSE){
    
    ## set outlier handling
    if(show.outliers) {
        outlier.shape = "x"
    } else {
        outlier.shape = NA
    } 
    
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
                data=gplot.data.boxplot, aes(
                    x=cond, y=ratio, fill=enhancer, alpha=batch, linetype=type),
                outlier.shape = outlier.shape ) +
                scale_shape_discrete(solid = FALSE) +
                labs(title=paste0(id, " log10 fdr-corrected p-value: ", fdr ))
        } else {
            gplot.boxplot <- ggplot() + geom_boxplot(
                data=gplot.data.boxplot, aes(
                    x=cond, y=ratio, fill=enhancer, linetype=type),
                outlier.shape = outlier.shape ) +
                scale_shape_discrete(solid = FALSE) +
                labs(title=paste0(id, " log10 fdr-corrected p-value: ", fdr ))
        }
    } else {
        if(!is.null(batch)) {
            gplot.boxplot <- ggplot() + geom_boxplot(
                data=gplot.data.boxplot, aes(
                    x=enhancer, y=ratio, fill=enhancer, alpha=batch, linetype=type),
                outlier.shape = outlier.shape ) +
                scale_shape_discrete(solid = FALSE) +
                labs(title=paste0(id, " log10 fdr-corrected p-value: ", fdr ))
        } else {
            gplot.boxplot <- ggplot() + geom_boxplot(
                data=gplot.data.boxplot, aes(
                    x=enhancer, y=ratio, fill=enhancer, linetype=type),
                outlier.shape = outlier.shape ) +
                scale_shape_discrete(solid = FALSE) +
                labs(title=paste0(id, " log10 fdr-corrected p-value: ", fdr ))
        }
    }
    if(length(unique(gplot.data.boxplot$batch)) <= length(cbbPalette)) {
        gplot.boxplot <- gplot.boxplot + scale_fill_manual(values = cbbPalette) + 
            scale_colour_manual(values = cbbPalette)
    }
    
    return(gplot.boxplot)
}

#' Scatter plot observed log-fold changes and p-values for condition test:
#' volcano plot
#' 
#' @param obj (MpraObject)
#' MpraObject from which enhancer is to be visualised.
#' 
#' @return (ggplot) scatter plot
#' 
#' @import ggplot2
#' 
#' @export
plotVolcano <- function(obj){
    
    ## extract model fits
    #TODO only extract coeffient
    stop()
    dfit <- getDNAFits(obj, enhancers=id, depth=FALSE, full=full)
    rfit <- getRNAFits(obj, enhancers=id, depth=FALSE, full=full)
    gplot.data.fit <- data.frame(
        ratio=log(dfit)-log(rfit),
        cond=obj@dnaAnnot[,condition],
        batch=obj@dnaAnnot[,batch],
        enhancer="model_fit",
        stringsAsFactors=FALSE
    )
    
    ## create ggplot
    gplot.volcano <- ggplot() + geom_boxplot(
        data=gplot.data, aes(x=lfc, y=padj) ) +
        labs(title="volcano plot" )
    
    return(gplot.volcano)
}

#' create a scatter plot of the model alpha (trascription rate) vs. the naive
#' estimator (RNA / DNA).
#' @param obj the MpraObject
#' @param condition if NULL (default) create a single plot of the model baseline
#' (intercept) vs the total data. Else, condition must be a valid factor in the 
#' object column annotations that was included in the design. In that case, a 
#' plot is generated for each level of the factor
#' @export
plotAlphaRatio <- function(obj, condition = NULL, logScale=TRUE) {
    if(is.null(condition)) {
        ratio <- rowMeans(obj@rnaCounts) / rowMeans(obj@dnaCounts)
        alpha <- getAlpha(obj)
        res <- ggplot(data = data.frame(ratio = ratio, alpha = alpha)) + 
            geom_point(mapping = aes(x = ratio, y = alpha)) + 
            geom_abline(intercept=0, slope=1) + 
            xlab("RNA / DNA") + ylab("alpha")
    } else {
        first <- TRUE
        res <- lapply(levels(as.factor(obj@dnaAnnot[,condition])), function(l) {
            if(first) {
                alpha <- getAlpha(obj)
                first <- FALSE
            } else {
                alpha <- getAlpha(obj, condition, l)
            }
            idx.dna <- obj@dnaAnnot[,condition] == l
            idx.rna <- obj@rnaAnnot[,condition] == l
            ratio <- (rowMeans(obj@rnaCounts[,idx.rna,drop=FALSE]) / 
                     rowMeans(obj@dnaCounts[,idx.dna,drop=FALSE]))
            
            ggplot(data = data.frame(ratio = ratio, alpha = alpha)) + 
                geom_point(mapping = aes(x = ratio, y = alpha)) + 
                geom_abline(intercept=0, slope=1) + 
                xlab("RNA / DNA") + ylab("alpha") + ggtitle(paste(condition, l))
        })
    }
    return(res)
}

plotObsExpDistributions <- function(obj, enhancer, RNA=TRUE) {
    if(length(enhancer) > 1) {
        stop("plase supply a single enhancer (index or name)")
    }
    if(RNA) {
        df <- data.frame(obs=obj@rnaCounts[enhancer,], 
                         expctd=getRNAFits(obj, enhancer))
    } else {
        df <- data.frame(obs=obj@dnaCounts[enhancer,], 
                         expctd=getDNAFits(obj, enhancer))
    }
    print(summary(df))
    ggplot(df) + 
        geom_histogram(aes(x = obs, color="blue")) + 
        geom_density(aes(x=expctd, fill="red"), alpha=0.1)
}

plotPvalCDF <- function(p, categories, adjusted=FALSE) {
    if(adjusted) {
        p <- p.adjust(p, "BH")
    }
    df <- data.frame(p=p, cat=categories)
    ggplot(df, aes(p, group=cat, color=cat)) + stat_ecdf(geom="step") + 
        xlab("P-value") + ylab("CDF")
    
}

plotObsExpDistributions <- function(obj, enhancer, rna=TRUE, KS = TRUE) {
    if(length(enhancer) > 1) {
        stop("plase supply a single enhancer (index or name)")
    }
    if(rna) {
        df <- data.frame(obs = obj@rnaCounts[enhancer,], 
                         exd = as.numeric(getRNAFits(obj, enhancer)))
    } else {
        df <- data.frame(obs = obj@dnaCounts[enhancer,], 
                         exd = as.numeric(getDNAFits(obj, enhancer)))
    }
    
    df <- df[df$obs > 0,]
    print(ks.test(df$obs, df$exd)$statistic)
    ggplot(df) + 
        geom_histogram(aes(obs, ..density.., fill="Observed")) + 
        geom_density(aes(exd, color="Expected"), size=2) + 
        scale_fill_manual(name=element_blank(), values = c("Observed"='grey33')) + 
        scale_colour_manual(name=element_blank(), values = c('Expected'='black')) + 
        theme(text=element_text(size=20), legend.position = "none") +
        # theme(text=element_text(size=20), legend.position = c(0.8, 0.8)) +
        xlab("counts")
    
    
    # df <- data.frame(counts=c(df$exd,df$obs), 
    #                  source=as.factor(c(rep(1,NROW(df)), rep(2,NROW(df)))))
    # levels(df$source) <- c("Expected", "Observed")
    # ggplot(df) + 
    #     geom_density(aes(counts, group=source, col=source, fill=source), 
    #                  size=1, alpha=0.25) +
    #     theme(text=element_text(size=20),
    #           legend.direction = "horizontal", legend.position = "bottom",
    #           legend.title = element_blank()) +
    #     scale_fill_manual(values = "")
}
