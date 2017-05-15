#' Run hypothesis test
#' @author David Sebastian Fischer
#' @export 
doHypothesisTests <- function(obj){
    
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
    if(!is.null(obj@vecCtrlIDs)){
        # Compute naive reference statistic assuming no sample structure,
        # ie. ignoring meta data.
        # T-test if there are multiple observations,
        # Z-test if there is one observation per sample.
        if(dim(obj@matRNACountsProc)[2] > 1) {
            dfDEAnalysis$ttest_ctrl <- calcTTest(
                obj, vecCaseIDs=row.names(obj@matRNACountsProc)[vecidxCase], 
                vecCtrlIDs=obj@vecCtrlIDs)
        } else {
            dfDEAnalysis$ztest_ctrl <- calcZTest(
                obj, vecCaseIDs=row.names(obj@matRNACountsProc)[vecidxCase], 
                vecCtrlIDs=obj@vecCtrlIDs)
        }
    }
    
    dfDEAnalysis <- dfDEAnalysis[match(obj@vecAllIDs,dfDEAnalysis$region),]
    rownames(dfDEAnalysis) <- obj@vecAllIDs
    dfDEAnalysis$allZero <- !(obj@vecAllIDs %in% rownames(obj@matRNACountsProc))
    
    obj@dfMPRAnalyzeResults <- dfDEAnalysis
    return(obj)
}

#' T-test on RNA/DNA ratios of case enhancer sequence versus control set
#' One-sided test for case ration being larger than control.
#' Appropriate to compare a set of ratios of a given enhancer with a set of control ratios,
#' i.e. multiple observations per enhancer.
#' @seealso calcZTest
#' @author David Sebastian Fischer
calcTTest <- function(obj, vecCaseIDs, vecCtrlIDs) {
    vecCtrlRatios <- do.call(c, lapply(vecCtrlIDs, function(id) {
        obj@matRNACountsProc[id,]/obj@matDNACountsProc[id,]
    })) 
    vecPvalTTest <- sapply(vecCaseIDs, function(id) {
        t.test(x=obj@matRNACountsProc[id,]/obj@matDNACountsProc[id,],
               y=vecCtrlRatios,
               alternative = "greater")$p.value # one sided test for whether case ratio is larger than enhancers
    })
    return(vecPvalTTest)
}

#' Z-score on ratios of cumulative RNA/DNA counts per case enhancer sequence and control set.
#' One-sided test for case ration being larger than control.
#' Appropriate to compare one ration of a given enhancer with a set of control ratios,
#' i.e. one observation per enhancer.
#' @seealso calcTTest
#' @author David Sebastian Fischer
calcZTest <- function(obj, vecCaseIDs, vecCtrlIDs) {
    vecCtrlRatios <- obj@matRNACountsProc[vecCtrlIDs,]/obj@matDNACountsProc[vecCtrlIDs,]
    vecPvalTTest <- sapply(vecCaseIDs, function(id) {
        pnorm(q=obj@matRNACountsProc[id,]/obj@matDNACountsProc[id,],
              mean = mean(vecCtrlRatios, na.rm=TRUE),
              sd = sd(vecCtrlRatios, na.rm=TRUE),
              lower.tail = FALSE, # one sided test for whether case ratio is larger than enhancers
              log.p = FALSE)
    })
    return(vecPvalTTest)
}