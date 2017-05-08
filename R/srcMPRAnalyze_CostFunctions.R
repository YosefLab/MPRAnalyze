#' Cost function for constant model
#' 
#' Log likelihood cost function for numerical optimisation of
#' constant model.
#' Implements log linker function for the constant mean parameter
#' and the batch correction factors.
#' Implements lower sensitivity bound of likelihood with respect
#' to constant mean parameter.
#' Implements upper and lower sensitivity bound of likelihood with respect
#' to batch correction factors.
#' 
#' @seealso Compiled version: \link{evalLogLik_comp}
#' 
#' @param vecTheta (numeric vector number of parameters to be estimated) 
#'    Constant model parameter and batch correction factor estimates.
#' @param vecCounts (numeric vector number of samples)
#'    Read count data.
#' @param scaDisp (scalar) Gene-wise 
#'    negative binomial dispersion hyper-parameter. 
#'    Dispersion parameter is taken from vecTheta 
#'    if this is NULL (numerical fitting of dispersion).
#' @param vecSizeFactors (numeric vector number of samples) 
#'    Model scaling factors for each sample which take
#'    sequencing depth into account (size factors).
#' @param lsvecidxBatch (list length number of confounding variables)
#' 		List of index vectors. 
#' 		One vector per confounding variable.
#' 		Each vector has one entry per sample with the index batch
#' 		within the given confounding variable of the given sample.
#' 		Batches are enumerated from 1 to number of batches.
#' @param vecboolObserved (bool vector number of samples)
#'    Whether sample is observed (finite and not NA).
#'     
#' @return scaLogLik (scalar) Value of cost function (loglikelihood) for given gene.
#' 
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
        if(scaDisp < 10^(-10)) scaDisp <- 10^(-10) 
    }
    scaSlopeRNAvsDNA <- exp(vecTheta[scaNParamUsed+1])
    if(scaSlopeRNAvsDNA < 10^(-10)) scaSlopeRNAvsDNA <- 10^(-10) 
    if(scaSlopeRNAvsDNA < 10^(-10)) scaSlopeRNAvsDNA <- 10^(-10) 
    scaNParamUsed <- scaNParamUsed + 1
    
    vecBatchFactors <- array(1, length(vecRNADepth))
    if(!is.null(lsvecidxBatchRNA)){
        for(vecidxConfounder in lsvecidxBatchRNA){
            scaNBatchFactors <- max(vecidxConfounder)-1 # Batches are counted from 1
            # Factor of first batch is one (constant), the remaining
            # factors scale based on the first batch.
            vecBatchFactorsConfounder <- c(1, exp(vecTheta[(scaNParamUsed+1):(scaNParamUsed+scaNBatchFactors)]))[vecidxConfounder]
            scaNParamUsed <- scaNParamUsed+scaNBatchFactors
            # Prevent batch factor shrinkage and explosion:
            vecBatchFactorsConfounder[vecBatchFactorsConfounder < 10^(-10)] <- 10^(-10)
            vecBatchFactorsConfounder[vecBatchFactorsConfounder > 10^(10)] <- 10^(10)
            
            vecBatchFactors <- vecBatchFactors*vecBatchFactorsConfounder
        }
    }
    
    if(boolBaselineCtrl){
        scaSlopeRNAvsDNACtrl <- exp(vecTheta[scaNParamUsed+1])
        if(scaSlopeRNAvsDNACtrl < 10^(-10)) scaSlopeRNAvsDNACtrl <- 10^(-10) 
        if(scaSlopeRNAvsDNACtrl < 10^(-10)) scaSlopeRNAvsDNACtrl <- 10^(-10) 
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
            vecBatchFactorsConfounder <- c(1, exp(vecTheta[(scaNParamUsed+1):(scaNParamUsed+scaNBatchFactors)]))[vecidxConfounder]
            scaNParamUsed <- scaNParamUsed+scaNBatchFactors
            # Prevent batch factor shrinkage and explosion:
            vecBatchFactorsConfounder[vecBatchFactorsConfounder < 10^(-10)] <- 10^(-10)
            vecBatchFactorsConfounder[vecBatchFactorsConfounder > 10^(10)] <- 10^(10)
            
            vecBatchFactorsCtrl <- vecBatchFactorsCtrl*vecBatchFactorsConfounder
        }
    }
    
    # Compute log likelihood under constant model by
    # adding log likelihood of model at each timepoint.
    scaLogLik <- sum(sapply(seq(1,dim(matRNACounts)[1]), function(i){
        vecboolObserved <- !is.na(matRNACounts[i,]) & !is.na(matDNAEst[i,])
        if(i==1) {
            scaCtrlCorrection <- 1
        } else {
            scaCtrlCorrection <- scaSlopeRNAvsDNACtrl*vecBatchFactorsCtrl[vecboolObserved]
        }
        scaLogLik <- sum(dnbinom(
            matRNACounts[i,vecboolObserved], 
            mu=matDNAEst[i,vecboolObserved]*
                scaSlopeRNAvsDNA*vecBatchFactors[vecboolObserved]*
                scaCtrlCorrection*
                vecRNADepth[vecboolObserved], 
            size=scaDisp, 
            log=TRUE))
        return(scaLogLik)
    }))
    
    # Maximise log likelihood: Return likelihood as value to optimisation routine
    return(scaLogLik)
}
evalLogLikRNA_pointDNAnbRNA_comp <- cmpfun(evalLogLikRNA_pointDNAnbRNA)

evalLogLikDNA_lnDNA <- function(
    vecTheta,
    vecDNACounts,
    vecboolObs,
    vecDNADepth,
    lsvecidxBatchDNA){
    
    scaNParamUsed <- 0
    scaMuDNA <- exp(vecTheta[scaNParamUsed+1])
    if(scaMuDNA < 10^(-10)) scaMuDNA <- 10^(-10) 
    if(scaMuDNA > 10^(10)) scaMuDNA <- 10^(10)
    scaNParamUsed <- scaNParamUsed + 1
    scaSdDNA <- exp(vecTheta[scaNParamUsed+1])
    if(scaSdDNA < 10^(-10)) scaSdDNA <- 10^(-10) 
    if(scaSdDNA > 10^(10)) scaSdDNA <- 10^(10)
    scaNParamUsed <- scaNParamUsed + 1
    
    vecDNAFit <- array(scaMuDNA, length(vecDNADepth))
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
            
            vecDNAFit <- vecDNAFit*vecBatchFacConf
        }
    }
    vecDNAHat <- vecDNAFit*vecDNADepth
    
    # Compute log likelihood under constant model by
    # adding log likelihood of model at each timepoint.
    vecLogLik <- dlnorm(
        x=vecDNACounts[vecboolObs],
        meanlog=log(vecDNAHat[vecboolObs]),
        sdlog=scaSdDNA,
        log=TRUE)
    scaLogLik <- sum(vecLogLik, na.rm=TRUE)
    
    # Maximise log likelihood: Return likelihood as value to optimisation routine
    return(scaLogLik)
}
evalLogLikDNA_lnDNA_comp <- cmpfun(evalLogLikDNA_lnDNA)

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
        if(scaSlopeRNAvsDNACtrl < 10^(-10)) scaSlopeRNAvsDNACtrl <- 10^(-10) 
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
evalLogLikDNARNA_gammaDNApoisRNA_comp <- cmpfun(evalLogLikDNARNA_gammaDNApoisRNA)

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
evalLogLikDNA_gammaDNApoisRNA_comp <- cmpfun(evalLogLikDNA_gammaDNApoisRNA)

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
        if(scaSlopeRNAvsDNACtrl < 10^(-10)) scaSlopeRNAvsDNACtrl <- 10^(-10) 
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
evalLogLikRNA_gammaDNApoisRNA_comp <- cmpfun(evalLogLikRNA_gammaDNApoisRNA)

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
evalLogLikDNARNA_gammaDNApoisRNA_direct_comp <- cmpfun(evalLogLikDNARNA_gammaDNApoisRNA_direct)
