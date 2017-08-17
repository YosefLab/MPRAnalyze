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

#' Objective for fitting DNA and RNA model under gammaDNApoisRNA framework
#' 
#' Used if DNA model was or was not pre-fit. Whether or not the control
#' DNA models are fit or taken as given is selected within the function.
#' 
#' @author David Sebastian Fischer
evalLogLikDNARNA_gammaDNApoisRNA <- function(
    vecTheta,
    matDNACounts,
    matRNACounts,
    lsvecidxDNARNAObs,
    lsvecidxDNAObs,
    vecDNADepth,
    vecRNADepth,
    boolBaselineCtrl,
    lsvecidxBatchRNA,
    lsvecidxBatchRNACtrl,
    lsvecidxBatchDNA,
    lsDNAModelFitsCtrl=NULL,
    lsObjPars){
    
    # commented code out to test run time of preproc and LL eval here
    #t1 <- system.time({for(jj in 1:100) {
    scaNParamUsed <- 0
    # DNA model
    matDNAModel <- matrix(NA, nrow = lsObjPars$scaNEnhancers,
                          ncol = 2)
    matCumulBatchFacDNA <- matrix(NA, nrow = lsObjPars$scaNEnhancers,
                                  ncol = dim(matDNACounts)[2] )
    for(i in lsObjPars$vecidxDNAModelsToFit){
        vecDNAModel <- exp(vecTheta[(scaNParamUsed+1):(scaNParamUsed+2)])
        scaNParamUsed <- scaNParamUsed + 2
        vecDNAModel[vecDNAModel < 10^(-10)] <- 10^(-10) 
        vecDNAModel[vecDNAModel > 10^(10)] <- 10^(10)
        
        vecBatchFacDNA <- rep(1, length(vecRNADepth))
        if(!is.null(lsvecidxBatchDNA)){
            for(j in seq(1, length(lsObjPars$vecNBatchFacDNA), by=1)){
                scaNBatchFactors <- lsObjPars$vecNBatchFacDNA[j]
                # Factor of first batch is one (constant), the remaining
                # factors scale based on the first batch.
                vecBatchFacConf <- c(1, exp(vecTheta[(scaNParamUsed+1):(scaNParamUsed+scaNBatchFactors)]))
                scaNParamUsed <- scaNParamUsed+scaNBatchFactors
                # Prevent batch factor shrinkage and explosion:
                vecBatchFacConf[vecBatchFacConf < 10^(-10)] <- 10^(-10)
                vecBatchFacConf[vecBatchFacConf > 10^(10)] <- 10^(10)
                
                vecBatchFacDNA <- vecBatchFacDNA*vecBatchFacConf[lsvecidxBatchDNA[[j]]]
            }
        }
        matDNAModel[i,] <- vecDNAModel
        matCumulBatchFacDNA[i,] <- vecBatchFacDNA
    }
    if(!is.null(lsDNAModelFitsCtrl)) {
        matDNAModel[seq(2, length(lsDNAModelFitsCtrl)+1, by=1),] <- do.call(rbind, lapply(lsDNAModelFitsCtrl, function(x) x$vecDNAModel ))
        matCumulBatchFacDNA[seq(2, length(lsDNAModelFitsCtrl)+1, by=1),] <- do.call(rbind, lapply(lsDNAModelFitsCtrl, function(x) x$vecCumulBatchFacDNA ))
    }
    
    scaRNAModel <- exp(vecTheta[(scaNParamUsed+1)])
    scaNParamUsed <- scaNParamUsed + 1
    scaRNAModel[scaRNAModel < 10^(-10)] <- 10^(-10) 
    scaRNAModel[scaRNAModel > 10^(10)] <- 10^(10)
    
    # RNA model
    vecRNAFit <- rep(scaRNAModel, lsObjPars$scaNSamples)
    if(!is.null(lsvecidxBatchRNA)){
        for(j in seq(1, length(lsObjPars$vecNBatchFacRNA), by=1)){
            scaNBatchFactors <- lsObjPars$vecNBatchFacRNA[j]
            # Factor of first batch is one (constant), the remaining
            # factors scale based on the first batch.
            vecBatchFacConf <- c(1, exp(vecTheta[(scaNParamUsed+1):(scaNParamUsed+scaNBatchFactors)]))
            scaNParamUsed <- scaNParamUsed+scaNBatchFactors
            # Prevent batch factor shrinkage and explosion:
            vecBatchFacConf[vecBatchFacConf < 10^(-10)] <- 10^(-10)
            vecBatchFacConf[vecBatchFacConf > 10^(10)] <- 10^(10)
            
            vecRNAFit <- vecRNAFit*vecBatchFacConf[lsvecidxBatchRNA[[j]]]
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
    vecBatchFactorsCtrl <- rep(1, lsObjPars$scaNSamples)
    if(!is.null(lsvecidxBatchRNACtrl)){
        for(j in seq(1, length(lsObjPars$vecNBatchFacRNACtrl), by=1)){
            scaNBatchFactors <- lsObjPars$vecNBatchFacRNACtrl[j]
            # Factor of first batch is one (constant), the remaining
            # factors scale based on the first batch.
            vecBatchFacConf <- c(1, exp(vecTheta[(scaNParamUsed+1):(scaNParamUsed+scaNBatchFactors)]))
            scaNParamUsed <- scaNParamUsed+scaNBatchFactors
            # Prevent batch factor shrinkage and explosion:
            vecBatchFacConf[vecBatchFacConf < 10^(-10)] <- 10^(-10)
            vecBatchFacConf[vecBatchFacConf > 10^(10)] <- 10^(10)
            
            vecBatchFactorsCtrl <- vecBatchFactorsCtrl*vecBatchFacConf[lsvecidxBatchRNACtrl[[j]]]
        }
    }
    #} })["elapsed"]
    #print(t1)
    
    matRNAFits <- rbind(c(vecRNAFit * vecRNADepth), # case
                        c(vecRNAFit * scaSlopeRNAvsDNACtrl * vecBatchFactorsCtrl *vecRNADepth)) # control
    #t2 <- system.time({for(jj in 1:100) {
    scaLogLikRNA <- sum(sapply(seq(1,lsObjPars$scaNEnhancers,by=1), function(i){
        return(sum(
            dnbinom(
                x=matRNACounts[i,lsvecidxDNARNAObs[[i]]], 
                mu=matDNAModel[i,1]/matDNAModel[i,2]*
                    matCumulBatchFacDNA[i,][lsvecidxDNARNAObs[[i]]]*
                    matRNAFits[as.numeric(i>1)+1,lsvecidxDNARNAObs[[i]]],
                size=matDNAModel[i,1],
                log=TRUE)
        ))
        # avoiding dnbinom via precompute gamma fun/factorial term
        #vecMu <- matDNAModel[i,1]/matDNAModel[i,2]*
        #    matCumulBatchFacDNA[i,][vecboolObsBoth]*
        #    vecRNAFit[vecboolObsBoth]*vecCtrlCorrection*
        #    vecRNADepth[vecboolObsBoth]
        #vecLogLikNB <- 1.45 + #c1
        #    matRNACounts[i,vecboolObsBoth] * (log(vecMu) - log(matDNAModel[i,1]+vecMu)) +
        #    matDNAModel[i,1] * (log(matDNAModel[i,1]) - log(matDNAModel[i,1]+vecMu))
    }))
    #} })["elapsed"]
    #print(t2)
    
    #t3 <- system.time({for(jj in 1:100){
    # only evaluate on the models that are actually re-estimated:
    # the likelihood of control enhancer DNA observations with pre-fit
    # models does not change here!
    scaLogLikDNA <- sum(sapply(lsObjPars$vecidxDNAModelsToFit, function(i){
        return(sum( 
            dgamma(
                x = matDNACounts[i,lsvecidxDNAObs[[i]]],
                shape = matDNAModel[i,1], 
                rate = matDNAModel[i,2] / 
                    (matCumulBatchFacDNA[i,][lsvecidxDNAObs[[i]]]* 
                         vecDNADepth[lsvecidxDNAObs[[i]]]),
                log=TRUE)
        ))
    }))
    #} })["elapsed"]
    #print(t3)
    #stop()
    
    return(scaLogLikRNA + scaLogLikDNA)
}
#' @author David Sebastian Fischer
evalLogLikDNARNA_gammaDNApoisRNA_comp <- cmpfun(evalLogLikDNARNA_gammaDNApoisRNA)

#' Objective for fitting DNA model under gammaDNApoisRNA framework
#' 
#' Used to pre-fit DNA model or to fit DNA-model only in iterative estimation.
#' 
#' @author David Sebastian Fischer
evalLogLikDNA_gammaDNApoisRNA <- function(
    vecTheta,
    vecDNACounts,
    vecRNACounts,
    vecidxDNARNAObs,
    vecidxDNAObs,
    vecDNADepth,
    vecRNAModelFitXRNADepth,
    vecRNAModelFit,
    lsvecidxBatchDNA,
    lsObjPars){
    
    # DNA model
    vecDNAModel <- exp(vecTheta[1:2])
    vecDNAModel[vecDNAModel < 10^(-10)] <- 10^(-10) 
    vecDNAModel[vecDNAModel > 10^(10)] <- 10^(10)
    scaNParamUsed <- 2
    
    vecBatchFacDNA <- rep(1, lsObjPars$scaNSamples)
    if(!is.null(lsvecidxBatchDNA)){
        for(j in seq(1, length(lsObjPars$vecNBatchFacDNA))){
            scaNBatchFactors <- lsObjPars$vecNBatchFacDNA[j]
            # Factor of first batch is one (constant), the remaining
            # factors scale based on the first batch.
            vecBatchFacConf <- c(1, exp(vecTheta[(scaNParamUsed+1):(scaNParamUsed+scaNBatchFactors)]))
            scaNParamUsed <- scaNParamUsed+scaNBatchFactors
            # Prevent batch factor shrinkage and explosion:
            vecBatchFacConf[vecBatchFacConf < 10^(-10)] <- 10^(-10)
            vecBatchFacConf[vecBatchFacConf > 10^(10)] <- 10^(10)
            
            vecBatchFacDNA <- vecBatchFacDNA*vecBatchFacConf[lsvecidxBatchDNA[[j]]]
        }
    }
    return(sum(
        sum(dnbinom(
            x=vecRNACounts[vecidxDNARNAObs], 
            mu=vecDNAModel[1]/vecDNAModel[2]*
                vecBatchFacDNA[vecidxDNARNAObs]*
                vecRNAModelFitXRNADepth[vecidxDNARNAObs], 
            size=vecDNAModel[1], 
            log=TRUE)) +
            sum(dgamma(
                x=vecDNACounts[vecidxDNAObs],
                shape=vecDNAModel[1], 
                rate=vecDNAModel[2] / 
                    (vecBatchFacDNA[vecidxDNAObs] * vecDNADepth[vecidxDNAObs]),
                log=TRUE))
    ))
}
#' @author David Sebastian Fischer
evalLogLikDNA_gammaDNApoisRNA_comp <- cmpfun(evalLogLikDNA_gammaDNApoisRNA)

#' Objective for fitting RNA model under gammaDNApoisRNA framework
#' 
#' Used to fit RNA model only in iterative estimation used for control DNA.
#' 
#' lsvecidxBatchRNACtrl depprecated for control DNA fitting only
#' 
#' @author David Sebastian Fischer
evalLogLikRNA_gammaDNApoisRNA <- function(
    vecTheta,
    matDNAModel,
    matCumulBatchFacDNA,
    matRNACounts,
    lsvecidxDNARNAObs,
    vecRNADepth,
    lsvecidxBatchRNA,
    boolBaselineCtrl=FALSE,
    lsvecidxBatchRNACtrl=NULL,
    lsObjPars){
    
    scaRNAModel <- exp(vecTheta[1])
    scaRNAModel[scaRNAModel < 10^(-10)] <- 10^(-10) 
    scaRNAModel[scaRNAModel > 10^(10)] <- 10^(10)
    scaNParamUsed <- 1
    
    # RNA model
    vecRNAModelFit <- rep(scaRNAModel, lsObjPars$scaNSamples)
    if(!is.null(lsvecidxBatchRNA)){
        for(j in seq(1,length(lsObjPars$vecNBatchFacRNA),by=1)){
            scaNBatchFactors <- lsObjPars$vecNBatchFacRNA[j]
            # Factor of first batch is one (constant), the remaining
            # factors scale based on the first batch.
            vecBatchFacConf <- c(1, exp(vecTheta[(scaNParamUsed+1):(scaNParamUsed+scaNBatchFactors)]))
            scaNParamUsed <- scaNParamUsed+scaNBatchFactors
            # Prevent batch factor shrinkage and explosion:
            vecBatchFacConf[vecBatchFacConf < 10^(-10)] <- 10^(-10)
            vecBatchFacConf[vecBatchFacConf > 10^(10)] <- 10^(10)
            
            vecRNAModelFit <- vecRNAModelFit*vecBatchFacConf[lsvecidxBatchRNA[[j]]]
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
    vecBatchFactorsCtrl <- rep(1, lsObjPars$scaNSamples)
    if(!is.null(lsvecidxBatchRNACtrl)){
        for(j in seq(1, length(lsObjPars$vecNBatchFacRNACtrl), by=1)){
            scaNBatchFactors <- lsObjPars$vecNBatchFacRNACtrl[j]
            # Factor of first batch is one (constant), the remaining
            # factors scale based on the first batch.
            vecBatchFacConf <- c(1, exp(vecTheta[(scaNParamUsed+1):(scaNParamUsed+scaNBatchFactors)]))
            scaNParamUsed <- scaNParamUsed+scaNBatchFactors
            # Prevent batch factor shrinkage and explosion:
            vecBatchFacConf[vecBatchFacConf < 10^(-10)] <- 10^(-10)
            vecBatchFacConf[vecBatchFacConf > 10^(10)] <- 10^(10)
            
            vecBatchFactorsCtrl <- vecBatchFactorsCtrl*vecBatchFacConf[lsvecidxBatchRNACtrl[[j]]]
        }
    }
    
    matRNAFits <- rbind(c(vecRNAModelFit * vecRNADepth), # case
                        c(vecRNAModelFit * scaSlopeRNAvsDNACtrl * vecBatchFactorsCtrl *vecRNADepth)) # control
    return(sum(
        sapply(seq(1, lsObjPars$scaNEnhancers, by=1), function(i){
            return(sum(
                vecLogLikNB <- dnbinom(
                    x=matRNACounts[i,lsvecidxDNARNAObs[[i]]], 
                    mu=matDNAModel[i,1]/matDNAModel[i,2]*
                        matCumulBatchFacDNA[i,lsvecidxDNARNAObs[[i]]]*
                        matRNAFits[as.numeric(i>1)+1,lsvecidxDNARNAObs[[i]]], 
                    size=matDNAModel[i,1], 
                    log=TRUE)
            ))
        })
    ))
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
