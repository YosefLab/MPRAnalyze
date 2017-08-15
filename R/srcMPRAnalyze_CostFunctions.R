#' @importFrom compiler cmpfun
#' @author David Sebastian Fischer
evalLogLikRNA_pointDNAnbRNA <- function(
    vecTheta,
    matRNACounts,
    matDNAEst,
    scaDisp, 
    vecRNADepth,
    boolBaselineCtrl,
    lsvecidxBatchRNA,
    lsvecidxBatchRNACtrl ){
    
    scaNParamUsed <- 0
    if(is.null(scaDisp)){
        scaDisp <- exp(vecTheta[scaNParamUsed+1])
        scaNParamUsed <- scaNParamUsed + 1
        if(scaDisp < 10^(-10)) scaDisp <- 10^(-10)
        if(scaDisp > 10^(10)) scaDisp <- 10^(10)
    }
    scaSlopeRNAvsDNA <- vecTheta[scaNParamUsed+1]
    if(scaSlopeRNAvsDNA < -23) scaSlopeRNAvsDNA <- -23 
    if(scaSlopeRNAvsDNA > 23) scaSlopeRNAvsDNA <- 23 
    scaNParamUsed <- scaNParamUsed + 1
    
    vecBatchFactors <- array(0, length(vecRNADepth))
    if(!is.null(lsvecidxBatchRNA)){
        for(vecidxConfounder in lsvecidxBatchRNA){
            scaNBatchFactors <- max(vecidxConfounder)-1 # Batches are counted from 1
            # Factor of first batch is one (constant), the remaining
            # factors scale based on the first batch.
            vecBatchFactorsConfounder <- c(0, vecTheta[(scaNParamUsed+1):(scaNParamUsed+scaNBatchFactors)])[vecidxConfounder]
            scaNParamUsed <- scaNParamUsed+scaNBatchFactors
            # Prevent batch factor shrinkage and explosion:
            vecBatchFactorsConfounder[vecBatchFactorsConfounder < -23] <- -23
            vecBatchFactorsConfounder[vecBatchFactorsConfounder > 23] <- 23
            
            vecBatchFactors <- vecBatchFactors+vecBatchFactorsConfounder
        }
    }
    
    if(boolBaselineCtrl){
        scaSlopeRNAvsDNACtrl <- vecTheta[scaNParamUsed+1]
        if(scaSlopeRNAvsDNACtrl < -23) scaSlopeRNAvsDNACtrl <- -23 
        if(scaSlopeRNAvsDNACtrl > 23) scaSlopeRNAvsDNACtrl <- 23 
        scaNParamUsed <- scaNParamUsed + 1
    } else {
        scaSlopeRNAvsDNACtrl <- 0
    }
    vecBatchFactorsCtrl <- array(0, length(vecRNADepth))
    if(!is.null(lsvecidxBatchRNACtrl)){
        for(vecidxConfounder in lsvecidxBatchRNACtrl){
            scaNBatchFactors <- max(vecidxConfounder)-1 # Batches are counted from 1
            # Factor of first batch is one (constant), the remaining
            # factors scale based on the first batch.
            vecBatchFactorsConfounder <- c(0, vecTheta[(scaNParamUsed+1):(scaNParamUsed+scaNBatchFactors)])[vecidxConfounder]
            scaNParamUsed <- scaNParamUsed+scaNBatchFactors
            # Prevent batch factor shrinkage and explosion:
            vecBatchFactorsConfounder[vecBatchFactorsConfounder < -23] <- -23
            vecBatchFactorsConfounder[vecBatchFactorsConfounder > 23] <- 23
            
            vecBatchFactorsCtrl <- vecBatchFactorsCtrl+vecBatchFactorsConfounder
        }
    }
    
    # Compute log likelihood under constant model by
    # adding log likelihood of model at each timepoint.
    scaLogLik <- sum(sapply(seq(1,dim(matRNACounts)[1]), function(i){
        vecboolObserved <- !is.na(matRNACounts[i,]) & !is.na(matDNAEst[i,])
        if(i==1) {
            scaCtrlCorrection <- 0
        } else {
            scaCtrlCorrection <- scaSlopeRNAvsDNACtrl+vecBatchFactorsCtrl[vecboolObserved]
        }
        scaLogLik <- sum(dnbinom(
            matRNACounts[i,vecboolObserved], 
            mu=matDNAEst[i,vecboolObserved]*
                exp(scaSlopeRNAvsDNA+vecBatchFactors[vecboolObserved]+
                scaCtrlCorrection) *
                vecRNADepth[vecboolObserved], 
            size=scaDisp, 
            log=TRUE))
        return(scaLogLik)
    }))
    
    # Maximise log likelihood: Return likelihood as value to optimisation routine
    return(scaLogLik)
}
#' @author David Sebastian Fischer
evalLogLikRNA_pointDNAnbRNA_comp <- cmpfun(evalLogLikRNA_pointDNAnbRNA)

#' @author David Sebastian Fischer
evalLogLikDNA_lnDNA <- function(
    vecTheta,
    vecDNACounts,
    vecboolObs,
    vecLogDNADepth,
    lsvecidxBatchDNA){
    
    scaNParamUsed <- 0
    scaMuDNA <- vecTheta[scaNParamUsed+1]
    if(scaMuDNA < -23) scaMuDNA <- -23 
    if(scaMuDNA > 23) scaMuDNA <- 23
    scaNParamUsed <- scaNParamUsed + 1
    scaSdDNA <- exp(vecTheta[scaNParamUsed+1])
    if(scaSdDNA < -23) scaSdDNA <- -23 
    if(scaSdDNA > 23) scaSdDNA <- 23
    scaNParamUsed <- scaNParamUsed + 1
    
    vecDNAFit <- array(scaMuDNA, length(vecLogDNADepth))
    if(!is.null(lsvecidxBatchDNA)){
        for(vecidxConfounder in lsvecidxBatchDNA){
            scaNBatchFactors <- max(vecidxConfounder)-1 # Batches are counted from 1
            # Factor of first batch is one (constant), the remaining
            # factors scale based on the first batch.
            vecBatchFacConf <- c(0, vecTheta[(scaNParamUsed+1):(scaNParamUsed+scaNBatchFactors)])[vecidxConfounder]
            scaNParamUsed <- scaNParamUsed+scaNBatchFactors
            # Prevent batch factor shrinkage and explosion:
            vecBatchFacConf[vecBatchFacConf < -23] <- -23
            vecBatchFacConf[vecBatchFacConf > 23] <- 23
            
            vecDNAFit <- vecDNAFit+vecBatchFacConf
        }
    }
    vecLogDNAHat <- vecDNAFit+vecLogDNADepth
    
    # Compute log likelihood under constant model by
    # adding log likelihood of model at each timepoint.
    scaLogLik <- sum(dlnorm(
        x=vecDNACounts[vecboolObs],
        meanlog=vecLogDNAHat[vecboolObs],
        sdlog=scaSdDNA,
        log=TRUE), na.rm=TRUE)
    
    # Maximise log likelihood: Return likelihood as value to optimisation routine
    return(scaLogLik)
}
#' @author David Sebastian Fischer
evalLogLikDNA_lnDNA_comp <- cmpfun(evalLogLikDNA_lnDNA)

#' Objective for fittin DNA and RNA model under gammaDNApoisRNA framework
#' 
#' Used if DNA model was not pre-fit.
#' 
#' @author David Sebastian Fischer
evalLogLikDNARNA_gammaDNApoisRNA <- function(
    vecTheta,
    matDNACounts,
    matRNACounts,
    vecDNADepth,
    vecRNADepth,
    boolBaselineCtrl,
    lsvecidxBatchRNA,
    lsvecidxBatchRNACtrl,
    lsvecidxBatchDNA,
    lsDNAModelFitsCtrl=NULL){
    
    if(is.null(lsDNAModelFitsCtrl)) {
        vecDNAModelsToFit <- seq(1, dim(matDNACounts)[1])
    } else {
        vecDNAModelsToFit <- 1
    }
    scaNParamUsed <- 0
    # DNA model
    lsvecDNAModel <- list()
    lsvecCumulBatchFacDNA <- list()
    for(i in vecDNAModelsToFit){
        vecDNAModel <- exp(vecTheta[(scaNParamUsed+1):(scaNParamUsed+2)])
        vecDNAModel[vecDNAModel < 10^(-10)] <- 10^(-10) 
        vecDNAModel[vecDNAModel > 10^(10)] <- 10^(10)
        scaNParamUsed <- scaNParamUsed + 2
        
        vecBatchFacDNA <- array(1, length(vecRNADepth))
        if(!is.null(lsvecidxBatchDNA)){
            for(vecidxConfounder in lsvecidxBatchDNA){
                scaNBatchFactors <- max(vecidxConfounder)-1 # Batches are counted from 1
                # Factor of first batch is one (constant), the remaining
                # factors scale based on the first batch.
                vecBatchFacConf <- c(1, exp(vecTheta[(scaNParamUsed+1):(scaNParamUsed+scaNBatchFactors)]))[vecidxConfounder]
                scaNParamUsed <- scaNParamUsed+scaNBatchFactors
                # Prevent batch factor shrinkage and explosion:
                vecBatchFacConf[vecBatchFacConf < 10^(-10)] <- 10^(-10)
                vecBatchFacConf[vecBatchFacConf > 10^(10)] <- 10^(10)
                
                vecBatchFacDNA <- vecBatchFacDNA*vecBatchFacConf
            }
        }
        lsvecDNAModel[[i]] <- vecDNAModel
        lsvecCumulBatchFacDNA[[i]] <- vecBatchFacDNA
    }
    if(!is.null(lsDNAModelFitsCtrl)) {
        for(i in seq(1, length(lsDNAModelFitsCtrl), by=1)){
            lsvecDNAModel[[i+1]] <- lsDNAModelFitsCtrl[[i]]$vecDNAModel
            lsvecCumulBatchFacDNA[[i+1]] <- lsDNAModelFitsCtrl[[i]]$vecCumulBatchFacDNA
        }
    }
    
    scaRNAModel <- exp(vecTheta[(scaNParamUsed+1)])
    scaRNAModel[scaRNAModel < 10^(-10)] <- 10^(-10) 
    scaRNAModel[scaRNAModel > 10^(10)] <- 10^(10)
    scaNParamUsed <- scaNParamUsed + 1
    
    # RNA model
    vecRNAFit <- array(scaRNAModel, length(vecRNADepth))
    if(!is.null(lsvecidxBatchRNA)){
        for(vecidxConfounder in lsvecidxBatchRNA){
            scaNBatchFactors <- max(vecidxConfounder)-1 # Batches are counted from 1
            # Factor of first batch is one (constant), the remaining
            # factors scale based on the first batch.
            vecBatchFacConf <- c(1, exp(vecTheta[(scaNParamUsed+1):(scaNParamUsed+scaNBatchFactors)]))[vecidxConfounder]
            scaNParamUsed <- scaNParamUsed+scaNBatchFactors
            # Prevent batch factor shrinkage and explosion:
            vecBatchFacConf[vecBatchFacConf < 10^(-10)] <- 10^(-10)
            vecBatchFacConf[vecBatchFacConf > 10^(10)] <- 10^(10)
            
            vecRNAFit <- vecRNAFit*vecBatchFacConf
        }
    }
    if(boolBaselineCtrl){
        scaSlopeRNAvsDNACtrl <- exp(vecTheta[scaNParamUsed+1])
        if(scaSlopeRNAvsDNACtrl < 10^(-10)) scaSlopeRNAvsDNACtrl <- 10^(-10) 
        if(scaSlopeRNAvsDNACtrl > 10^(10)) scaSlopeRNAvsDNACtrl <- 10^(10) 
        scaNParamUsed <- scaNParamUsed + 1
    } else {
        scaSlopeRNAvsDNACtrl <- 1
    }
    vecBatchFactorsCtrl <- array(1, length(vecRNADepth))
    if(!is.null(lsvecidxBatchRNACtrl)){
        for(vecidxConfounder in lsvecidxBatchRNACtrl){
            scaNBatchFactors <- max(vecidxConfounder)-1 # Batches are counted from 1
            # Factor of first batch is one (constant), the remaining
            # factors scale based on the first batch.
            vecBatchFacConf <- c(1, exp(vecTheta[(scaNParamUsed+1):(scaNParamUsed+scaNBatchFactors)]))[vecidxConfounder]
            scaNParamUsed <- scaNParamUsed+scaNBatchFactors
            # Prevent batch factor shrinkage and explosion:
            vecBatchFacConf[vecBatchFacConf < 10^(-10)] <- 10^(-10)
            vecBatchFacConf[vecBatchFacConf > 10^(10)] <- 10^(10)
            
            vecBatchFactorsCtrl <- vecBatchFactorsCtrl*vecBatchFacConf
        }
    }
    
    scaLogLik <- sum(sapply(seq(1,dim(matDNACounts)[1]), function(i){
        vecboolObsBoth <- !is.na(matDNACounts[i,]) & !is.na(matRNACounts[i,]) 
        if(i==1) {
            vecCtrlCorrection <- rep(1, sum(vecboolObsBoth))
        } else {
            vecCtrlCorrection <- scaSlopeRNAvsDNACtrl*vecBatchFactorsCtrl[vecboolObsBoth]
        }
        vecLogLikNB <- dnbinom(
            x=matRNACounts[i,vecboolObsBoth], 
            mu=lsvecDNAModel[[i]][1]/lsvecDNAModel[[i]][2]*
                lsvecCumulBatchFacDNA[[i]][vecboolObsBoth]*
                vecRNAFit[vecboolObsBoth]*vecCtrlCorrection*
                vecRNADepth[vecboolObsBoth], 
            size=lsvecDNAModel[[i]][1], #1/a?
            log=TRUE)
        vecboolObsDNA <- !is.na(matDNACounts[i,])
        vecLogLikGamma <- dgamma(
            x=matDNACounts[i,vecboolObsDNA],
            shape=lsvecDNAModel[[i]][1], 
            rate=lsvecDNAModel[[i]][2]/
                (lsvecCumulBatchFacDNA[[i]][vecboolObsDNA]*
                     vecDNADepth[vecboolObsDNA]),
            log=TRUE)
        scaLogLik <- sum(vecLogLikNB) + sum(vecLogLikGamma)
        return(scaLogLik)
    }))
    return(scaLogLik)
}
#' @author David Sebastian Fischer
evalLogLikDNARNA_gammaDNApoisRNA_comp <- cmpfun(evalLogLikDNARNA_gammaDNApoisRNA)

#' Objective for fittin DNA model under gammaDNApoisRNA framework
#' 
#' Used to pre-fit DNA model or to fit DNA-model only in iterative estimation.
#' 
#' @author David Sebastian Fischer
evalLogLikDNA_gammaDNApoisRNA <- function(
    vecTheta,
    vecDNACounts,
    vecRNACounts,
    vecboolObsDNA,
    vecboolObsBoth,
    vecDNADepth,
    vecRNADepth,
    vecRNAModelFit,
    lsvecidxBatchDNA){
    
    scaNParamUsed <- 0
    # DNA model
    vecDNAModel <- exp(vecTheta[(scaNParamUsed+1):(scaNParamUsed+2)])
    vecDNAModel[vecDNAModel < 10^(-10)] <- 10^(-10) 
    vecDNAModel[vecDNAModel > 10^(10)] <- 10^(10)
    scaNParamUsed <- scaNParamUsed + 2
    
    vecBatchFacDNA <- array(1, length(vecRNADepth))
    if(!is.null(lsvecidxBatchDNA)){
        for(vecidxConfounder in lsvecidxBatchDNA){
            scaNBatchFactors <- max(vecidxConfounder)-1 # Batches are counted from 1
            # Factor of first batch is one (constant), the remaining
            # factors scale based on the first batch.
            vecBatchFacConf <- c(1, exp(vecTheta[(scaNParamUsed+1):(scaNParamUsed+scaNBatchFactors)]))[vecidxConfounder]
            scaNParamUsed <- scaNParamUsed+scaNBatchFactors
            # Prevent batch factor shrinkage and explosion:
            vecBatchFacConf[vecBatchFacConf < 10^(-10)] <- 10^(-10)
            vecBatchFacConf[vecBatchFacConf > 10^(10)] <- 10^(10)
            
            vecBatchFacDNA <- vecBatchFacDNA*vecBatchFacConf
        }
    }
    vecLogLikNB <- dnbinom(
        x=vecRNACounts[vecboolObsBoth], 
        mu=vecDNAModel[1]/vecDNAModel[2]*
            vecBatchFacDNA[vecboolObsBoth]*
            vecRNAModelFit[vecboolObsBoth]*vecRNADepth[vecboolObsBoth], 
        size=vecDNAModel[1], #1/a?
        log=TRUE)
    vecLogLikGamma <- dgamma(
        x=vecDNACounts[vecboolObsDNA],
        shape=vecDNAModel[1], 
        rate=vecDNAModel[2]/(vecBatchFacDNA[vecboolObsDNA]*
                                 vecDNADepth[vecboolObsDNA]),
        log=TRUE)
    scaLogLik <- sum(vecLogLikNB) +
        sum(vecLogLikGamma)
    return(scaLogLik)
}
#' @author David Sebastian Fischer
evalLogLikDNA_gammaDNApoisRNA_comp <- cmpfun(evalLogLikDNA_gammaDNApoisRNA)

#' Objective for fittin RNA model under gammaDNApoisRNA framework
#' 
#' Used if DNA model was pre-fit or to fit RNA model only in iterative estimation.
#' 
#' @author David Sebastian Fischer
evalLogLikRNA_gammaDNApoisRNA <- function(
    vecTheta,
    lsvecDNAModel,
    lsvecCumulBatchFacDNA,
    matDNACounts,
    matRNACounts,
    vecDNADepth,
    vecRNADepth,
    lsvecidxBatchRNA,
    boolBaselineCtrl,
    lsvecidxBatchRNACtrl){
    
    scaNParamUsed <- 0
    scaRNAModel <- exp(vecTheta[(scaNParamUsed+1)])
    scaRNAModel[scaRNAModel < 10^(-10)] <- 10^(-10) 
    scaRNAModel[scaRNAModel > 10^(10)] <- 10^(10)
    scaNParamUsed <- scaNParamUsed + 1
    
    # RNA model
    vecRNAModelFit <- array(scaRNAModel, length(vecRNADepth))
    if(!is.null(lsvecidxBatchRNA)){
        for(vecidxConfounder in lsvecidxBatchRNA){
            scaNBatchFactors <- max(vecidxConfounder)-1 # Batches are counted from 1
            # Factor of first batch is one (constant), the remaining
            # factors scale based on the first batch.
            vecBatchFacConf <- c(1, exp(vecTheta[(scaNParamUsed+1):(scaNParamUsed+scaNBatchFactors)]))[vecidxConfounder]
            scaNParamUsed <- scaNParamUsed+scaNBatchFactors
            # Prevent batch factor shrinkage and explosion:
            vecBatchFacConf[vecBatchFacConf < 10^(-10)] <- 10^(-10)
            vecBatchFacConf[vecBatchFacConf > 10^(10)] <- 10^(10)
            
            vecRNAModelFit <- vecRNAModelFit*vecBatchFacConf
        }
    }
    vecRNAModelFitCtrl <- rep(1, length(vecRNADepth))
    if(boolBaselineCtrl){
        scaSlopeRNAvsDNACtrl <- exp(vecTheta[scaNParamUsed+1])
        if(scaSlopeRNAvsDNACtrl < 10^(-10)) scaSlopeRNAvsDNACtrl <- 10^(-10) 
        if(scaSlopeRNAvsDNACtrl > 10^(10)) scaSlopeRNAvsDNACtrl <- 10^(10) 
        scaNParamUsed <- scaNParamUsed + 1
        vecRNAModelFitCtrl <- vecRNAModelFitCtrl*scaSlopeRNAvsDNACtrl
    } else {
        scaSlopeRNAvsDNACtrl <- 1
    }
    if(!is.null(lsvecidxBatchRNACtrl)){
        for(vecidxConfounder in lsvecidxBatchRNACtrl){
            scaNBatchFactors <- max(vecidxConfounder)-1 # Batches are counted from 1
            # Factor of first batch is one (constant), the remaining
            # factors scale based on the first batch.
            vecBatchFacConf <- c(1, exp(vecTheta[(scaNParamUsed+1):(scaNParamUsed+scaNBatchFactors)]))[vecidxConfounder]
            scaNParamUsed <- scaNParamUsed+scaNBatchFactors
            # Prevent batch factor shrinkage and explosion:
            vecBatchFacConf[vecBatchFacConf < 10^(-10)] <- 10^(-10)
            vecBatchFacConf[vecBatchFacConf > 10^(10)] <- 10^(10)
            
            vecRNAModelFitCtrl <- vecRNAModelFitCtrl*vecBatchFacConf
        }
    }
    
    scaLogLik <- sum(sapply(seq(1,length(lsvecDNAModel)), function(i){
        vecboolObsBoth <- !is.na(matDNACounts[i,]) & !is.na(matRNACounts[i,]) 
        if(i==1) {
            vecCtrlCorrection <- rep(1, sum(vecboolObsBoth))
        } else {
            vecCtrlCorrection <- vecRNAModelFitCtrl[vecboolObsBoth]
        }
        vecLogLikNB <- dnbinom(
            x=matRNACounts[i,vecboolObsBoth], 
            mu=lsvecDNAModel[[i]][1]/lsvecDNAModel[[i]][2]*
                lsvecCumulBatchFacDNA[[i]][vecboolObsBoth]*
                vecRNAModelFit[vecboolObsBoth]*
                vecCtrlCorrection*
                vecRNADepth[vecboolObsBoth], 
            size=lsvecDNAModel[[i]][1], #1/a?
            log=TRUE)
        scaLogLik <- sum(vecLogLikNB)
        return(scaLogLik)
    }))
    return(scaLogLik)
}
#' @author David Sebastian Fischer
evalLogLikRNA_gammaDNApoisRNA_comp <- cmpfun(evalLogLikRNA_gammaDNApoisRNA)

#' @author David Sebastian Fischer
evalLogLikDNARNA_gammaDNApoisRNA_direct <- function(
    vecTheta,
    matDNACounts,
    matRNACounts,
    vecDNADepth,
    vecRNADepth,
    lsvecDNAModel,
    lsvecCumulBatchFacDNA,
    vecRNAModelFit,
    vecRNAModelFitCtrl){
    
    vecLogLik <- sapply(seq(1,dim(matDNACounts)[1]), function(i){
        vecboolObsBoth <- !is.na(matDNACounts[i,]) & !is.na(matRNACounts[i,]) 
        if(i==1) {
            vecCtrlCorrection <- rep(1, sum(vecboolObsBoth))
        } else {
            vecCtrlCorrection <- vecRNAModelFitCtrl[vecboolObsBoth]
        }
        vecLogLikNB <- dnbinom(
            x=matRNACounts[i,vecboolObsBoth], 
            mu=lsvecDNAModel[[i]][1]/lsvecDNAModel[[i]][2]*
                lsvecCumulBatchFacDNA[[i]][vecboolObsBoth]*
                vecRNAModelFit[vecboolObsBoth]*vecCtrlCorrection*
                vecRNADepth[vecboolObsBoth], 
            size=lsvecDNAModel[[i]][1], #1/a?
            log=TRUE)
        vecboolObsDNA <- !is.na(matDNACounts[i,])
        vecLogLikGamma <- dgamma(
            x=matDNACounts[i,vecboolObsDNA],
            shape=lsvecDNAModel[[i]][1], 
            rate=lsvecDNAModel[[i]][2]/
                (lsvecCumulBatchFacDNA[[i]][vecboolObsDNA]*
                     vecDNADepth[vecboolObsDNA]),
            log=TRUE)
        scaLogLik <- sum(vecLogLikNB) + 
            sum(vecLogLikGamma)
        return(scaLogLik)
    })
    return(sum(vecLogLik))
}
#' @author David Sebastian Fischer
evalLogLikDNARNA_gammaDNApoisRNA_direct_comp <- cmpfun(evalLogLikDNARNA_gammaDNApoisRNA_direct)
