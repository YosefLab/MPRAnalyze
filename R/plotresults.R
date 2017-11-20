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
#' 
#' @return (gplot) Boxplot graphic.
#' 
#' @import ggplot2
#' 
#' @export
plotBoxplots <- function(obj, id, condition=NULL, batch=NULL, full=TRUE){
    
    ## extract observations and format into data-frame of rna:dna ratios
    # case enhancer observations
    gplot.data.obs <- data.frame(
        ratio=(log(obj@rnaCounts[id,])-log(obj@rnaDepth))-
            (log(obj@dnaCounts[id,])-log(obj@dnaDepth)),
        cond=obj@colAnnot[,condition],
        batch=obj@colAnnot[,batch],
        enhancer=id,
        stringsAsFactors=FALSE
    )
    # control enhancer observations
    if(!is.null(obj@controls)){
        gplot.data.obs.ctrl <- data.frame(
            ratio=(log(obj@rnaCounts[obj@controls,])-log(obj@rnaDepth))-
                (log(obj@dnaCounts[obj@controls,])-log(obj@dnaDepth)),
            cond=obj@colAnnot[,condition],
            batch=obj@colAnnot[,batch],
            enhancer="control",
            stringsAsFactors=FALSE
        )
        gplot.data.obs <- rbind(gplot.data.obs, gplot.data.obs.ctrl)
    }
    
    ## extract model fits
    dfit <- as.vector(getDNAFits(obj, enhancers=id, depth=FALSE, full=full))
    rfit <- as.vector(getRNAFits(obj, enhancers=id, depth=FALSE, full=full))
    gplot.data.fit <- data.frame(
        ratio=log(rfit)-log(dfit),
        cond=obj@colAnnot[,condition],
        batch=obj@colAnnot[,batch],
        enhancer="model_fit",
        stringsAsFactors=FALSE
    )
    
    ## create ggplot
    # colour-blind palette
    cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                    "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    gplot.data <- rbind(gplot.data.obs, gplot.data.fit)
    gplot.boxplot <- ggplot() + geom_boxplot(
        data=gplot.data.obs, aes(x=cond, y=ratio, fill=batch, alpha=batch),
        outlier.shape = "x" ) +
        geom_point(data = gplot.data.fit, aes(x=cond, y=ratio, shape=batch)) +
        labs(title=paste0(id, " log10 fdr-corrected p-value: ", 
                          round(log(obj@results[id,]$fdr.pval)/log(10),2)) )
    if(length(unique(gplot.data.obs$batch <= length(cbbPalette)))) {
        gplot.boxplot <- gplot.boxplot + scale_fill_manual(values = cbbPalette)
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
        cond=obj@colAnnot[,condition],
        batch=obj@colAnnot[,batch],
        enhancer="model_fit",
        stringsAsFactors=FALSE
    )
    
    ## create ggplot
    gplot.volcano <- ggplot() + geom_boxplot(
        data=gplot.data, aes(x=lfc, y=padj) ) +
        labs(title="volcano plot" )

    return(gplot.volcano)
}