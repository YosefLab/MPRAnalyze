#' export
plotBoxplots <- function(obj, strID){
    
    vecboolObs <- !is.na(obj@matDNACountsProc[strID,]) & !is.na(obj@matRNACountsProc[strID,])
    dfgplotObs <- data.frame(
        ratio=log(obj@matRNACountsProc[strID,]/obj@vecRNADepth)-
            log(obj@matDNACountsProc[strID,]/obj@vecDNADepth),
        timecateg=obj@dfAnnotationProc$timeCateg,
        batch=obj@dfAnnotationProc$batch,
        enhancer="case"
    )
    if(!is.null(obj@vecCtrlIDs)){
        dfgplotObsCtrl <- data.frame(
            ratio=do.call(c, lapply(obj@vecCtrlIDs, function(strIDctrl) {
                log(obj@matRNACountsProc[strIDctrl,]/obj@vecRNADepth)-
                    log(obj@matDNACountsProc[strIDctrl,]/obj@vecDNADepth) 
            })),
            timecateg=rep(obj@dfAnnotationProc$timeCateg, length(obj@vecCtrlIDs)),
            batch=rep(obj@dfAnnotationProc$batch, length(obj@vecCtrlIDs)),
            enhancer="scrambled"
        )
        dfgplotObs <- rbind(dfgplotObs, dfgplotObsCtrl)
    }
    dfgplotFit <- data.frame(
        ratio=log(obj@lsModelFitsFull[[strID]]$vecFitRNA[vecboolObs]/
                      obj@lsModelFitsFull[[strID]]$vecFitDNA[vecboolObs]),
        timecateg=obj@dfAnnotationProc$timeCateg[vecboolObs],
        batch=obj@dfAnnotationProc$batch[vecboolObs],
        enhancer="fit"
    )
    dfgplotObs <- rbind(dfgplotObs, dfgplotFit)
    gplotBoxplot <- ggplot() + geom_boxplot(
        data=dfgplotObs, aes(x=timecateg, y=ratio, fill=enhancer) ) +
        labs(title=paste0(strID, " log10 q-value: ", 
                          round(log(obj$dfMPRAnalyzeResults[strID,]$padj)/log(10),2)) )
    
    return(gplotBoxplot)
}