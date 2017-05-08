runDEAnalysis <- function(obj){
    
    # Extract fit reporters
    vecLLfull <- sapply(obj@lsModelFitsFull, function(fit) fit$scaLL )
    vecLLred <- sapply(obj@lsModelFitsRed, function(fit) fit$scaLL )
    
    vecDFfull <- sapply(obj@lsModelFitsFull, function(fit) fit$scaDF )
    vecDFred <- sapply(obj@lsModelFitsRed, function(fit) fit$scaDF )
    
    vecCONVfull <- sapply(obj@lsModelFitsFull, function(fit) fit$scaConvergence )
    vecCONVred <- sapply(obj@lsModelFitsRed, function(fit) fit$scaConvergence )
    
    # Compute P-value from chi-square distributed deviance
    scaDeltaDF <- vecDFfull[1] - vecDFred[1]
    vecDeviance <- 2*(vecLLfull - vecLLred)
    vecPvalue <- pchisq(vecDeviance,df=scaDeltaDF,lower.tail=FALSE)
    # Benjamini-Hochberg multiple testing correction
    vecPvalueBH <- p.adjust(vecPvalue, method = "BH")
    
    if(is.null(obj@vecCtrlIDs)){
        vecidxCase <- seq(1, dim(obj@matDNACountsProc)[1])
    } else {
        vecidxComb <- seq(1, dim(obj@matDNACountsProc)[1])
        vecidxCtrl <- match(obj@vecCtrlIDs, rownames(obj@matDNACountsProc))
        vecidxCase <- setdiff(vecidxComb, vecidxCtrl)
    }
    
    # Prepare output table
    dfDEAnalysis<- data.frame(
        region=row.names(obj@matRNACountsProc[vecidxCase,,drop=FALSE]),
        p=vecPvalue,
        padj=vecPvalueBH,
        loglik_full=vecLLfull,
        loglik_red=vecLLred,
        df_full=vecDFfull,
        df_red=vecDFred,
        meanRNA=apply(obj@matRNACountsProc[vecidxCase,,drop=FALSE], 1, function(g) mean(g, na.rm=TRUE) ),
        meanDNA=apply(obj@matDNACountsProc[vecidxCase,,drop=FALSE], 1, function(g) mean(g, na.rm=TRUE) ),
        converge_full=vecCONVfull,
        converge_red=vecCONVred,
        stringsAsFactors = FALSE,
        row.names = NULL)
    
    dfDEAnalysis <- dfDEAnalysis[match(obj@vecAllIDs,dfDEAnalysis$region),]
    rownames(dfDEAnalysis) <- obj@vecAllIDs
    dfDEAnalysis$allZero <- !(obj@vecAllIDs %in% rownames(obj@matRNACountsProc))
    
    obj@dfMPRAnalyzeResults <- dfDEAnalysis
    return(obj)
}