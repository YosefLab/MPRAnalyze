#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++     Cost Functions    +++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#' Optimisation wrapper for RNA model estimation under pointDNAnbRNA framework
#' 
#' This wrapper runs and handles input to and output from optim().
#' 
#' @param vecCounts (numeric vector number of samples)
#'    Read count data.
#' @param scaDisp (scalar) Gene-wise 
#'    negative binomial dispersion hyper-parameter.
#' @param vecSizeFactors (numeric vector number of samples) 
#'    Model scaling factors for each sample which take
#'    sequencing depth into account (size factors).
#' @param lsvecidxBatch (list length number of confounding variables)
#' 		List of index vectors. 
#' 		One vector per confounding variable.
#' 		Each vector has one entry per sample with the index batch
#' 		within the given confounding variable of the given sample.
#' 		Batches are enumerated from 1 to number of batches.
#' @param MAXIT (scalar) [Default 1000] 
#'    Maximum number of BFGS iterations for model fitting with \link{optim}.
#' @param RELTOL (scalar) [Default 10^(-8)]
#'    Maximum relative change in loglikelihood to reach convergence in
#'    numerical optimisation by BFGS in \link{optim}.
#' @param trace (scalar) [Defaul 0]
#'    Reporting parameter of \link{optim}.
#' @param REPORT (scalar) [Default 10]
#'    Reporting parameter of \link{optim}.
#'     
#' @return (list) List of constant fit parameters and results.
#'    \itemize{
#'      \item scaMu (scalar) Maximum likelihood estimator of
#'      negative binomial mean parameter.
#'      \item lsvecBatchFactors (list length number of confounders)
#'      List of vectors of scalar batch correction factors for each sample.
#'      These are also maximum likelihood estimators.
#'      NULL if no confounders given.
#'      \item scaDispParam (scalar) Dispersion parameter estimate
#'      used in fitting (hyper-parameter).
#'      \item scaLL (scalar) Loglikelihood of data under maximum likelihood
#'      estimator model.
#'      \item scaConvergence (scalar) 
#'      Convergence status of optim on constant model.
#'    }
#'    
#' @author David Sebastian Fischer
fitRNA_pointDNAnbRNA <- function(
    matRNACounts,
    matDNAEst,
    scaDisp,
    vecRNADepth,
    lsvecidxBatchRNA,
    boolBaselineCtrl=FALSE,
    lsvecidxBatchRNACtrl=NULL,
    MAXIT=1000,
    RELTOL=10^(-8)){
    
    # Initialise parameters
    # Mean and batch correction model
    scaSlopeRNAvsDNA <- median(matRNACounts/matDNAEst, na.rm=TRUE)
    vecParamGuess <- scaSlopeRNAvsDNA
    if(!is.null(lsvecidxBatchRNA)){
        for(vecidxConfounder in lsvecidxBatchRNA){
            vecParamGuess <- c(vecParamGuess, rep(0, length(unique(vecidxConfounder))-1))
        }
    }
    if(boolBaselineCtrl) vecParamGuess <- c(vecParamGuess, 0)
    if(!is.null(lsvecidxBatchRNACtrl)){
        for(vecidxConfounder in lsvecidxBatchRNACtrl){
            vecParamGuess <- c(vecParamGuess, rep(0, length(unique(vecidxConfounder))-1))
        }
    }
    # Dispersion model
    if(is.null(scaDisp)){
        scaDispGuess <- 1
        vecParamGuess <- c(log(scaDispGuess), vecParamGuess)
    }
    
    lsFit <- tryCatch({
        optim(
            par=vecParamGuess,
            fn=evalLogLikRNA_pointDNAnbRNA_comp,
            matRNACounts=matRNACounts,
            matDNAEst=matDNAEst,
            scaDisp=scaDisp, 
            vecRNADepth=vecRNADepth,
            lsvecidxBatchRNA=lsvecidxBatchRNA,
            boolBaselineCtrl=boolBaselineCtrl,
            lsvecidxBatchRNACtrl=lsvecidxBatchRNACtrl,
            method="BFGS",
            control=list(maxit=MAXIT,
                         reltol=RELTOL,
                         fnscale=-1)
        )[c("par","value","convergence")]
    }, error=function(strErrorMsg){
        print(paste0("ERROR: Fitting null model: fitRNA_pointDNAnbRNA()."))
        print(paste0("vecParamGuess ", paste(vecParamGuess,collapse=" ")))
        print(paste0("matRNACounts ", paste(matRNACounts,collapse=" ")))
        print(paste0("matDNAEst ", paste(matDNAEst,collapse=" ")))
        print(paste0("scaDisp ", paste(scaDisp,collapse=" ")))
        print(paste0("vecRNADepth ", paste(vecRNADepth,collapse=" ")))
        print(paste0("lsvecidxBatchRNA ", paste(lsvecidxBatchRNA,collapse=" ")))
        print(paste0("MAXIT ", MAXIT))
        print(strErrorMsg)
        stop(strErrorMsg)
    })
    
    # Extract parameter estimates
    scaNParamUsed <- 0
    # Extract dispersion
    if(is.null(scaDisp)){
        scaDisp <- exp(lsFit$par[scaNParamUsed+1])
        scaNParamUsed <- scaNParamUsed + 1
        if(scaDisp < 10^(-10)) scaDisp <- 10^(-10) 
        if(scaDisp > 10^(10)) scaDisp <- 10^(10) 
    }
    # Extract mean
    scaSlopeRNAvsDNA <- exp(lsFit$par[scaNParamUsed+1])
    if(scaSlopeRNAvsDNA < 10^(-10)) scaSlopeRNAvsDNA <- 10^(-10) 
    if(scaSlopeRNAvsDNA > 10^(10)) scaSlopeRNAvsDNA <- 10^(10) 
    scaNParamUsed <- scaNParamUsed + 1
    # - Extract batch correction factors and
    # - Compute fitted values per observation
    vecFitRNA <- matDNAEst[1,]*vecRNADepth*scaSlopeRNAvsDNA
    if(!is.null(lsvecidxBatchRNA)){
        lsvecBatchFactorsRNA <- list()
        for(i in seq(1, length(lsvecidxBatchRNA))){
            scaNBatchFactors <- max(lsvecidxBatchRNA[[i]])-1 # Batches are counted from 1
            # Factor of first batch is one (constant), the remaining
            # factors scale based on the first batch.
            vecBatchFactors <- c(1, exp(lsFit$par[(scaNParamUsed+1):(scaNParamUsed+scaNBatchFactors)]))
            scaNParamUsed <- scaNParamUsed+scaNBatchFactors
            # Catch boundary of likelihood domain on batch factor space:
            vecBatchFactors[vecBatchFactors < 10^(-10)] <- 10^(-10)
            vecBatchFactors[vecBatchFactors > 10^(10)] <- 10^(10)
            lsvecBatchFactorsRNA[[i]] <- vecBatchFactors
            
            vecFitRNA <- vecFitRNA*vecBatchFactors[lsvecidxBatchRNA[[i]]]
        }
    } else { 
        lsvecBatchFactorsRNA <- NULL 
    }
    
    if(boolBaselineCtrl) {
        scaSlopeRNAvsDNACtrl <- exp(lsFit$par[scaNParamUsed+1])
        if(scaSlopeRNAvsDNACtrl < 10^(-10)) scaSlopeRNAvsDNACtrl <- 10^(-10) 
        if(scaSlopeRNAvsDNACtrl > 10^(10)) scaSlopeRNAvsDNACtrl <- 10^(10) 
        scaNParamUsed <- scaNParamUsed + 1
    } else { 
        scaSlopeRNAvsDNACtrl <- NULL 
    }
    if(!is.null(lsvecidxBatchRNACtrl)){
        lsvecBatchFactorsCtrl <- list()
        for(i in seq(1, length(lsvecidxBatchRNACtrl))){
            scaNBatchFactors <- max(lsvecidxBatchRNACtrl[[i]])-1 # Batches are counted from 1
            # Factor of first batch is one (constant), the remaining
            # factors scale based on the first batch.
            vecBatchFactors <- c(1, exp(lsFit$par[(scaNParamUsed+1):(scaNParamUsed+scaNBatchFactors)]))
            scaNParamUsed <- scaNParamUsed+scaNBatchFactors
            # Catch boundary of likelihood domain on batch factor space:
            vecBatchFactors[vecBatchFactors < 10^(-10)] <- 10^(-10)
            vecBatchFactors[vecBatchFactors > 10^(10)] <- 10^(10)
            lsvecBatchFactorsCtrl[[i]] <- vecBatchFactors
        }
    } else { 
        lsvecBatchFactorsCtrl <- NULL 
    }
    
    return(list(vecExprModel=scaSlopeRNAvsDNA,
                vecExprModelCtrl=scaSlopeRNAvsDNACtrl,
                lsvecBatchFactorsRNA=lsvecBatchFactorsRNA,
                lsvecBatchFactorsRNACtrl=lsvecBatchFactorsCtrl,
                scaDisp=scaDisp,
                vecFitRNA=vecFitRNA,
                scaDF=length(lsFit$par),
                scaLL=lsFit$value,
                scaConvergence=lsFit$convergence))
}

#' Optimisation wrapper for DNA model estimation under pointDNA*RNA framework
#' 
#' This wrapper runs and handles input to and output from optim().
#' 
#' @author David Sebastian Fischer
fitDNA_lnDNA <- function(
    vecDNACounts,
    vecDNADepth,
    lsvecidxBatchDNA,
    MAXIT=1000,
    RELTOL=10^(-8)){
    
    scaMuDNA <- mean(vecDNACounts, na.rm=TRUE)
    scaSdDNA <- 1
    vecParamGuess <- c(log(scaMuDNA), log(scaSdDNA))
    # DNA model
    if(!is.null(lsvecidxBatchDNA)){ # DNA batch model
        for(vecidxConfounder in lsvecidxBatchDNA){
            #vecMeanEst <- tapply(vecDNACounts/vecDNADepth, vecidxConfounder, mean, na.rm=TRUE)
            #vecFacEst <- vecMeanEst/scaMuDNA
            #vecFacEst <- (vecFacEst/vecFacEst[1])[2:length(vecidxConfounder)]
            #vecParamGuess <- c(vecParamGuess, vecFacEst)
            vecParamGuess <- c(vecParamGuess, rep(0, length(unique(vecidxConfounder))-1))
        }
    }
    
    lsFit <- tryCatch({
        optim(
            par=vecParamGuess,
            fn=evalLogLikDNA_lnDNA_comp,
            vecDNACounts=vecDNACounts,
            vecboolObs=!is.na(vecDNACounts),
            vecLogDNADepth=log(vecDNADepth),
            lsvecidxBatchDNA=lsvecidxBatchDNA,
            method="BFGS",
            control=list(maxit=1000,
                         reltol=RELTOL,
                         fnscale=-1)
        )[c("par","value","convergence")]
    }, error=function(strErrorMsg){
        print(paste0("ERROR: Fitting null model: fitDNA_lnDNA()."))
        print(paste0("vecParamGuess ", paste(vecParamGuess,collapse=" ")))
        print(paste0("matDNACounts ", paste(matDNACounts,collapse=" ")))
        print(paste0("vecDNADepth ", paste(vecDNADepth,collapse=" ")))
        print(paste0("lsvecidxBatchDNA ", paste(lsvecidxBatchDNA,collapse=" ")))
        print(paste0("MAXIT ", MAXIT))
        print(strErrorMsg)
    })
    # Extract parameter estimates
    scaNParamUsed <- 0
    scaMuDNA <- exp(lsFit$par[scaNParamUsed+1])
    if(scaMuDNA < 10^(-10)) scaMuDNA <- 10^(-10)
    if(scaMuDNA > 10^(10)) scaMuDNA <- 10^(10)
    scaNParamUsed <- scaNParamUsed+1
    scaSdDNA <- exp(lsFit$par[scaNParamUsed+1])
    if(scaSdDNA < 10^(-10)) scaSdDNA <- 10^(-10)
    if(scaSdDNA > 10^(10)) scaSdDNA <- 10^(10)
    scaNParamUsed <- scaNParamUsed+1
    
    # Extract batch correction factors and
    # Compute fitted values per observation
    vecFitDNA <- rep(scaMuDNA, length(vecDNACounts))
    if(!is.null(lsvecidxBatchDNA)){
        lsvecBatchFactorsDNA <- list()
        for(i in seq(1, length(lsvecidxBatchDNA))){
            scaNBatchFactors <- max(lsvecidxBatchDNA[[i]])-1 # Batches are counted from 1
            # Factor of first batch is one (constant), the remaining
            # factors scale based on the first batch.
            vecBatchFactors <- c(1, exp(lsFit$par[(scaNParamUsed+1):(scaNParamUsed+scaNBatchFactors)]))
            scaNParamUsed <- scaNParamUsed+scaNBatchFactors
            # Catch boundary of likelihood domain on batch factor space:
            vecBatchFactors[vecBatchFactors < 10^(-10)] <- 10^(-10)
            vecBatchFactors[vecBatchFactors > 10^(10)] <- 10^(10)
            lsvecBatchFactorsDNA[[i]] <- vecBatchFactors
            
            vecFitDNA <- vecFitDNA*vecBatchFactors[lsvecidxBatchDNA[[i]]]
        }
    } else { 
        lsvecBatchFactorsDNA <- NULL
    }
    
    return(list(scaMuDNA=scaMuDNA,
                scaSdDNA=scaSdDNA,
                lsvecBatchFactorsDNA=lsvecBatchFactorsDNA,
                vecFitDNA=vecFitDNA*vecDNADepth,
                vecFitDNAHat=vecFitDNA,
                scaDF=length(lsFit$par),
                scaLL=lsFit$value,
                scaConvergence=lsFit$convergence))
}

#' Optimisation wrapper for DNA and RNA model estimation under gammaDNApoisRNA framework
#' 
#' This wrapper runs and handles input to and output from optim().
#' 
#' @seealso fitDNARNA_gammaDNApoisRNA_coordascent_coordascent() for the
#' same statistical framework but coordinate ascent based parameter estimation
#' which makes fitting in high dimensional parameter spaces feasible.
#' 
#' @author David Sebastian Fischer
fitDNARNA_gammaDNApoisRNA <- function(
    matDNACounts,
    matRNACounts,
    vecDNADepth,
    vecRNADepth,
    lsvecidxBatchDNA,
    lsvecidxBatchRNA,
    boolBaselineCtrl=FALSE,
    lsvecidxBatchRNACtrl=NULL,
    lsDNAModelFitsCtrl=NULL,
    MAXIT=1000,
    RELTOL=10^(-8)){
    
    # DNA model
    if(is.null(lsDNAModelFitsCtrl)) {
        vecidxDNAModelsToFit <- seq(1, dim(matDNACounts)[1])
    } else {
        vecidxDNAModelsToFit <- 1
    }
    
    vecParamGuess <- c()
    for(i in vecidxDNAModelsToFit) {
        # Gamma distr parameters, mean right and disp 1
        vecDNAModel <- c(1,1/mean(matDNACounts[i,], na.rm=TRUE)) 
        vecParamGuess <- c(vecParamGuess, log(vecDNAModel))
        if(!is.null(lsvecidxBatchDNA)){ # DNA batch model
            for(vecidxConfounder in lsvecidxBatchDNA){
                vecParamGuess <- c(vecParamGuess, rep(0, length(unique(vecidxConfounder))-1))
            }
        }
    }
    scaRNAModel <- 1
    vecParamGuess <- c(vecParamGuess, log(scaRNAModel)) # baseline #RNA per plasmid
    if(!is.null(lsvecidxBatchRNA)){ # DNA batch model
        for(vecidxConfounder in lsvecidxBatchRNA){
            vecParamGuess <- c(vecParamGuess, rep(0, length(unique(vecidxConfounder))-1))
        }
    }
    if(boolBaselineCtrl) vecParamGuess <- c(vecParamGuess, 0)
    if(!is.null(lsvecidxBatchRNACtrl)){
        for(vecidxConfounder in lsvecidxBatchRNACtrl){
            vecParamGuess <- c(vecParamGuess, rep(0, length(unique(vecidxConfounder))-1))
        }
    }
    
    lsFit <- tryCatch({
        optim(
            par=vecParamGuess,
            fn=evalLogLikDNARNA_gammaDNApoisRNA_comp,
            matDNACounts=matDNACounts,
            matRNACounts=matRNACounts,
            vecDNADepth=vecDNADepth,
            vecRNADepth=vecRNADepth,
            lsvecidxBatchDNA=lsvecidxBatchDNA,
            lsvecidxBatchRNA=lsvecidxBatchRNA,
            boolBaselineCtrl=boolBaselineCtrl,
            lsvecidxBatchRNACtrl=lsvecidxBatchRNACtrl,
            lsDNAModelFitsCtrl=lsDNAModelFitsCtrl,
            method="BFGS",
            control=list(maxit=1000,
                         reltol=RELTOL,
                         fnscale=-1) #, trace=TRUE, REPORT=1)
        )[c("par","value","convergence")]
    }, error=function(strErrorMsg){
        print(paste0("ERROR: Fitting null model: fitConstModel()."))
        print(paste0("vecParamGuess ", paste(vecParamGuess,collapse=" ")))
        print(paste0("matDNACounts ", paste(matDNACounts,collapse=" ")))
        print(paste0("vecDNADepth ", paste(vecDNADepth,collapse=" ")))
        print(paste0("matRNACounts ", paste(matRNACounts,collapse=" ")))
        print(paste0("vecRNADepth ", paste(vecRNADepth,collapse=" ")))
        print(paste0("lsvecidxBatchRNA ", paste(lsvecidxBatchRNA,collapse=" ")))
        print(paste0("lsvecidxBatchDNA ", paste(lsvecidxBatchDNA,collapse=" ")))
        print(paste0("MAXIT ", MAXIT))
        print(strErrorMsg)
    })
    
    # Extract parameter estimates
    scaNParamUsed <- 0
    lsvecDNAModel <- list()
    lslsvecBatchFacDNA <- list()
    lsvecCumulBatchFacDNA <- list()
    for(i in vecidxDNAModelsToFit) {
        vecDNAModel <- exp(lsFit$par[(scaNParamUsed+1):(scaNParamUsed+2)])
        scaNParamUsed <- scaNParamUsed + 2
        vecDNAModel[vecDNAModel < 10^(-10)] <- 10^(-10)
        vecDNAModel[vecDNAModel > 10^(10)] <- 10^(10)
        
        # DNA model
        vecBatchFacDNA <- rep(1, length(vecDNADepth))
        if(!is.null(lsvecidxBatchDNA)){
            lsvecBatchFacDNA <- list()
            for(j in seq(1, length(lsvecidxBatchDNA))){
                scaNBatchFactors <- max(lsvecidxBatchDNA[[j]])-1 # Batches are counted from 1
                # Factor of first batch is one (constant), the remaining
                # factors scale based on the first batch.
                vecBatchFactors <- c(1, exp(lsFit$par[(scaNParamUsed+1):(scaNParamUsed+scaNBatchFactors)]))
                scaNParamUsed <- scaNParamUsed+scaNBatchFactors
                # Catch boundary of likelihood domain on batch factor space:
                vecBatchFactors[vecBatchFactors < 10^(-10)] <- 10^(-10)
                vecBatchFactors[vecBatchFactors > 10^(10)] <- 10^(10)
                lsvecBatchFacDNA[[j]] <- vecBatchFactors
                
                vecBatchFacDNA <- vecBatchFacDNA*vecBatchFactors[lsvecidxBatchDNA[[j]]]
            }
        } else { 
            lsvecBatchFacDNA <- NULL
        }
        lsvecDNAModel[[i]] <- vecDNAModel
        lslsvecBatchFacDNA[[i]] <- lsvecBatchFacDNA
        lsvecCumulBatchFacDNA[[i]] <- vecBatchFacDNA
    }
    vecDNAModelCase <- lsvecDNAModel[[1]]
    lsvecBatchFacDNACase <- lslsvecBatchFacDNA[[1]]
    vecBatchFacDNACase <- lsvecCumulBatchFacDNA[[1]]
    vecFitDNACase <- vecDNAModelCase[1]/vecDNAModelCase[2]*vecBatchFacDNACase*vecDNADepth
    
    scaSlopeRNAvsDNA <- exp(lsFit$par[(scaNParamUsed+1)])
    scaNParamUsed <- scaNParamUsed + 1
    scaSlopeRNAvsDNA[scaSlopeRNAvsDNA < 10^(-10)] <- 10^(-10)
    scaSlopeRNAvsDNA[scaSlopeRNAvsDNA > 10^(10)] <- 10^(10)
    
    # RNA model
    vecFitRNAModel <- rep(scaSlopeRNAvsDNA, length(vecRNADepth))
    if(!is.null(lsvecidxBatchRNA)){
        lsvecBatchFactorsRNA <- list()
        for(i in seq(1, length(lsvecidxBatchRNA))){
            scaNBatchFactors <- max(lsvecidxBatchRNA[[i]])-1 # Batches are counted from 1
            # Factor of first batch is one (constant), the remaining
            # factors scale based on the first batch.
            vecBatchFactors <- c(1, exp(lsFit$par[(scaNParamUsed+1):(scaNParamUsed+scaNBatchFactors)]))
            scaNParamUsed <- scaNParamUsed+scaNBatchFactors
            # Catch boundary of likelihood domain on batch factor space:
            vecBatchFactors[vecBatchFactors < 10^(-10)] <- 10^(-10)
            vecBatchFactors[vecBatchFactors > 10^(10)] <- 10^(10)
            lsvecBatchFactorsRNA[[i]] <- vecBatchFactors
            
            vecFitRNAModel <- vecFitRNAModel*vecBatchFactors[lsvecidxBatchRNA[[i]]]
        }
    } else { 
        lsvecBatchFactorsRNA <- NULL
    }
    if(boolBaselineCtrl) {
        scaSlopeRNAvsDNACtrl <- exp(lsFit$par[scaNParamUsed+1])
        if(scaSlopeRNAvsDNACtrl < 10^(-10)) scaSlopeRNAvsDNACtrl <- 10^(-10) 
        if(scaSlopeRNAvsDNACtrl > 10^(10)) scaSlopeRNAvsDNACtrl <- 10^(10) 
        scaNParamUsed <- scaNParamUsed + 1
    } else { 
        scaSlopeRNAvsDNACtrl <- NULL 
    }
    if(!is.null(lsvecidxBatchRNACtrl)){
        lsvecBatchFactorsRNACtrl <- list()
        for(i in seq(1, length(lsvecidxBatchRNACtrl))){
            scaNBatchFactors <- max(lsvecidxBatchRNACtrl[[i]])-1 # Batches are counted from 1
            # Factor of first batch is one (constant), the remaining
            # factors scale based on the first batch.
            vecBatchFactors <- c(1, exp(lsFit$par[(scaNParamUsed+1):(scaNParamUsed+scaNBatchFactors)]))
            scaNParamUsed <- scaNParamUsed+scaNBatchFactors
            # Catch boundary of likelihood domain on batch factor space:
            vecBatchFactors[vecBatchFactors < 10^(-10)] <- 10^(-10)
            vecBatchFactors[vecBatchFactors > 10^(10)] <- 10^(10)
            lsvecBatchFactorsRNACtrl[[i]] <- vecBatchFactors
        }
    } else { 
        lsvecBatchFactorsRNACtrl <- NULL 
    }
    
    vecFitRNA <- vecDNAModelCase[1]/vecDNAModelCase[2] *
        vecBatchFacDNACase * vecFitRNAModel * vecRNADepth
    
    return(list( vecDNAModel=vecDNAModelCase,
                 lsvecBatchFacDNA=lsvecBatchFacDNACase,
                 vecExprModel=scaSlopeRNAvsDNA,
                 vecExprModelCtrl=scaSlopeRNAvsDNACtrl,
                 lsvecBatchFactorsRNA=lsvecBatchFactorsRNA,
                 lsvecBatchFactorsRNACtrl=lsvecBatchFactorsRNACtrl,
                 vecFitDNA=vecFitDNACase,
                 vecFitRNA=vecFitRNA,
                 scaDisp=vecDNAModelCase[1],
                 scaDF=length(lsFit$par),
                 scaLL=lsFit$value,
                 scaConvergence=lsFit$convergence))
}

#' Optimisation wrapper for DNA and RNA model estimation under gammaDNApoisRNA framework
#' with coordinate ascent estimation scheme
#' 
#' This wrapper runs coordinate ascent as an iteration over BFGS-based estimation
#' of parameter groups exploiting conditionaly independencies.
#' Within each iteration, this wrapper runs and handles input to and output from optim().
#' Coordinate ascent allows parameter estimation in large parameter spaces,
#' e.g. if many control sequences are given.
#' 
#' @seealso fitDNARNA_gammaDNApoisRNA_coordascent() the same statistical
#' framework without coordinate ascent for more exact and faster
#' estimation in low dimensional parameter spaces.
#' 
#' @author David Sebastian Fischer
fitDNARNA_gammaDNApoisRNA_coordascent <- function(
    matDNACounts,
    matRNACounts,
    vecDNADepth,
    vecRNADepth,
    lsvecidxBatchDNA,
    lsvecidxBatchRNA,
    boolBaselineCtrl=FALSE,
    lsvecidxBatchRNACtrl=NULL,
    lsModelsDNA=NULL,
    lsModelRNA=NULL,
    MAXIT=1000,
    RELTOL=10^(-6) ){
    
    # Initialise
    # DNA models
    if(is.null(lsModelsDNA)){
        lsModelsDNA <- list()
        # DNA batch model
        if(!is.null(lsvecidxBatchDNA)){ 
            lsvecBatchFacDNA <- list()
            for(i in seq(1, length(lsvecidxBatchDNA))){
                lsvecBatchFacDNA[[i]] <- rep(1, length(unique(lsvecidxBatchDNA[[i]])))
            }
        } else{
            lsvecBatchFacDNA <- NULL
        }
        for(i in seq(1, dim(matDNACounts)[1])){
            # Gamma distr parameters, mean right and disp 1
            vecDNAModel <- c(1,1/mean(matDNACounts[i,], na.rm=TRUE))
            lsModelsDNA[[i]] <- list()
            lsModelsDNA[[i]]$vecDNAModel <- vecDNAModel
            lsModelsDNA[[i]]$lsvecBatchFacDNA <- lsvecBatchFacDNA
        }
    }
    
    if(is.null(lsModelRNA)){
        # RNA model
        scaSlopeRNAvsDNA <- mean(matRNACounts/matDNACounts, na.rm=TRUE)
        if(!is.null(lsvecidxBatchRNA)){ # RNA batch model
            lsvecBatchFacRNA <- list()
            for(i in seq(1, length(lsvecidxBatchRNA))){
                lsvecBatchFacRNA[[i]] <- rep(1, length(unique(lsvecidxBatchRNA[[i]])))
            }
        } else{
            lsvecBatchFacRNA <- NULL
        }
        scaSlopeRNAvsDNACtrl <- 1
        if(!is.null(lsvecidxBatchRNACtrl)){ # RNA batch model
            lsvecBatchFacRNACtrl <- list()
            for(i in seq(1, length(lsvecidxBatchRNACtrl))){
                lsvecBatchFacRNACtrl[[i]] <- rep(1, length(unique(lsvecidxBatchRNACtrl[[i]])))
            }
        } else{
            lsvecBatchFacRNACtrl <- NULL
        }
        lsModelRNA <- list(
            vecExprModel=scaSlopeRNAvsDNA,
            vecExprModelCtrl=scaSlopeRNAvsDNACtrl,
            lsvecBatchFactorsRNA=lsvecBatchFacRNA,
            lsvecBatchFactorsRNACtrl=lsvecBatchFacRNACtrl,
            vecRNAModelFit=rep(scaSlopeRNAvsDNA, length(vecRNADepth)),
            vecRNAModelFitCtrl=rep(1, length(vecRNADepth)) )
    }
    
    # Coordinate ascent: group parameters into batches
    # of conditionally independent parameters and estiamte
    # each batch in turns with BFGS until LL converges over all.
    scaLLNew <- -Inf
    scaLLOld <- -Inf
    scaPrec <- 10^(-6)
    scaIter <- 0
    
    print("Enter itertative coordinate ascent")
    while(scaLLNew+(scaPrec*scaLLNew) > scaLLOld | scaIter==0){
        scaIter <- scaIter + 1
        scaLLOld <- scaLLNew
        
        # DNA model
        lsModelsDNA <- lapply(seq(1, dim(matDNACounts)[1]), function(i){
            if(i == 1) { 
                vecModel <- lsModelRNA$vecRNAModelFit
            } else {
                vecModel <- lsModelRNA$vecRNAModelFit*
                    lsModelRNA$vecRNAModelFitCtrl
            }
            
            # Gamma distr parameters, mean right and disp 1
            vecParamGuessDNA <- log(lsModelsDNA[[i]]$vecDNAModel)
            if(!is.null(lsvecidxBatchDNA)){ # DNA batch model
                for(j in seq(1, length(lsvecidxBatchDNA))){
                    vecBatchFac <- lsModelsDNA[[i]]$lsvecBatchFacDNA[[j]]
                    vecParamGuessDNA <- 
                        c(vecParamGuessDNA, log(vecBatchFac[2:length(vecBatchFac)]))
                }
            }
            lsFitDNA <- tryCatch({
                optim(
                    par=vecParamGuessDNA,
                    fn=evalLogLikDNA_gammaDNApoisRNA_comp,
                    vecDNACounts=matDNACounts[i,],
                    vecRNACounts=matRNACounts[i,],
                    vecboolObsDNA=!is.na(matDNACounts[i,]),
                    vecboolObsBoth=!is.na(matDNACounts[i,]) & !is.na(matRNACounts[i,]),
                    vecDNADepth=vecDNADepth,
                    vecRNADepth=vecRNADepth,
                    vecRNAModelFit=vecModel,
                    lsvecidxBatchDNA=lsvecidxBatchDNA,
                    method="BFGS",
                    control=list(maxit=1000,
                                 reltol=RELTOL,
                                 fnscale=-1)
                )[c("par","value","convergence")]
            }, error=function(strErrorMsg){
                print(paste0("ERROR: Fitting DNA model: fitDNARNA_gammaDNApoisRNA_coordascent()."))
                print(paste0("vecParamGuessDNA ", paste(vecParamGuessDNA,collapse=" ")))
                print(paste0("matDNACounts ", paste(matDNACounts,collapse=" ")))
                print(paste0("vecDNADepth ", paste(vecDNADepth,collapse=" ")))
                print(paste0("matRNACounts ", paste(matRNACounts,collapse=" ")))
                print(paste0("vecRNADepth ", paste(vecRNADepth,collapse=" ")))
                print(paste0("lsvecidxBatchRNA ", paste(lsvecidxBatchRNA,collapse=" ")))
                print(paste0("lsvecidxBatchDNA ", paste(lsvecidxBatchDNA,collapse=" ")))
                print(paste0("MAXIT ", MAXIT))
                print(strErrorMsg)
            })
            
            # Extract parameter estimates
            scaNParamUsed <- 0
            vecDNAModel <- exp(lsFitDNA$par[(scaNParamUsed+1):(scaNParamUsed+2)])
            scaNParamUsed <- scaNParamUsed + 2
            vecDNAModel[vecDNAModel < 10^(-10)] <- 10^(-10)
            vecDNAModel[vecDNAModel > 10^(10)] <- 10^(10)
            
            vecBatchFacDNA <- rep(1, length(vecDNADepth))
            if(!is.null(lsvecidxBatchDNA)){
                lsvecBatchFacDNA <- list()
                for(i in seq(1, length(lsvecidxBatchDNA))){
                    scaNBatchFactors <- max(lsvecidxBatchDNA[[i]])-1 # Batches are counted from 1
                    # Factor of first batch is one (constant), the remaining
                    # factors scale based on the first batch.
                    vecBatchFactors <- c(1, exp(lsFitDNA$par[(scaNParamUsed+1):(scaNParamUsed+scaNBatchFactors)]))
                    scaNParamUsed <- scaNParamUsed+scaNBatchFactors
                    # Catch boundary of likelihood domain on batch factor space:
                    vecBatchFactors[vecBatchFactors < 10^(-10)] <- 10^(-10)
                    vecBatchFactors[vecBatchFactors > 10^(10)] <- 10^(10)
                    lsvecBatchFacDNA[[i]] <- vecBatchFactors
                    
                    vecBatchFacDNA <- vecBatchFacDNA*vecBatchFactors[lsvecidxBatchDNA[[i]]]
                }
            } else { 
                lsvecBatchFacDNA <- NULL
            }
            vecFitDNAHat <- vecDNAModel[1]/vecDNAModel[2]*vecBatchFacDNA
            return(list(
                vecDNAModel=vecDNAModel,
                lsvecBatchFacDNA=lsvecBatchFacDNA,
                vecCumulBatchFacDNA=vecBatchFacDNA,
                vecFitDNAHat=vecFitDNAHat,
                scaDFDNA=length(lsFitDNA$par),
                scaLL=lsFitDNA$value,
                scaConvergence=lsFitDNA$convergence
            ))
        })
        vecFitDNA <- lsModelsDNA[[1]]$vecFitDNAHat*vecDNADepth
        
        scaLLNew <- evalLogLikDNARNA_gammaDNApoisRNA_direct(
            matDNACounts=matDNACounts,
            matRNACounts=matRNACounts,
            vecDNADepth=vecDNADepth,
            vecRNADepth=vecRNADepth,
            lsvecDNAModel=lapply(lsModelsDNA, function(f) f$vecDNAModel ),
            lsvecCumulBatchFacDNA=lapply(lsModelsDNA, function(f) f$vecCumulBatchFacDNA),
            vecRNAModelFit=lsModelRNA$vecRNAModelFit,
            vecRNAModelFitCtrl=lsModelRNA$vecRNAModelFitCtrl)
        print(paste0(scaIter, ". iteration DNA done with LL=", round(scaLLNew,5)))
        
        # RNA model
        vecParamGuessRNA <- log(lsModelRNA$vecExprModel) # baseline #RNA per plasmid
        if(!is.null(lsvecidxBatchRNA)){ # DNA batch model
            for(j in seq(1, length(lsvecidxBatchRNA))){
                vecBatchFac <- lsModelRNA$lsvecBatchFactorsRNA[[j]]
                vecParamGuessRNA <- 
                    c(vecParamGuessRNA, log(vecBatchFac[2:length(vecBatchFac)]))
            }
        }
        if(boolBaselineCtrl) {
            vecParamGuessRNA <- c(vecParamGuessRNA, 
                                  log(lsModelRNA$vecExprModelCtrl))
        }
        if(!is.null(lsvecidxBatchRNACtrl)){
            for(j in seq(1, length(lsvecidxBatchRNACtrl))){
                vecBatchFac <- lsModelRNA$lsvecBatchFactorsRNACtrl[[j]]
                vecParamGuessRNA <- 
                    c(vecParamGuessRNA, log(vecBatchFac[2:length(vecBatchFac)]))
            }
        }
        lsFitRNA <- tryCatch({
            optim(
                par=vecParamGuessRNA,
                fn=evalLogLikRNA_gammaDNApoisRNA_comp,
                lsvecDNAModel=lapply(lsModelsDNA, function(f) f$vecDNAModel ),
                lsvecCumulBatchFacDNA=lapply(lsModelsDNA, function(f) f$vecCumulBatchFacDNA),
                matDNACounts=matDNACounts,
                matRNACounts=matRNACounts,
                vecDNADepth=vecDNADepth,
                vecRNADepth=vecRNADepth,
                lsvecidxBatchRNA=lsvecidxBatchRNA,
                boolBaselineCtrl=boolBaselineCtrl,
                lsvecidxBatchRNACtrl=lsvecidxBatchRNACtrl,
                method="BFGS",
                control=list(maxit=1000,
                             reltol=RELTOL,
                             fnscale=-1)
            )[c("par","value","convergence")]
        }, error=function(strErrorMsg){
            print(paste0("ERROR: Fitting DNA model: fitDNARNA_gammaDNApoisRNA_coordascent()."))
            print(paste0("vecParamGuessRNA ", paste(vecParamGuessRNA,collapse=" ")))
            print(paste0("lsvecidxBatchRNA ", paste(lsvecidxBatchRNA,collapse=" ")))
            print(paste0("lsvecidxBatchDNA ", paste(lsvecidxBatchDNA,collapse=" ")))
            print(paste0("MAXIT ", MAXIT))
            print(strErrorMsg)
        })
        
        scaNParamUsed <- 0
        scaSlopeRNAvsDNA <- exp(lsFitRNA$par[(scaNParamUsed+1)])
        scaNParamUsed <- scaNParamUsed + 1
        scaSlopeRNAvsDNA[scaSlopeRNAvsDNA < 10^(-10)] <- 10^(-10)
        scaSlopeRNAvsDNA[scaSlopeRNAvsDNA > 10^(10)] <- 10^(10)
        
        vecRNAModelFit <- rep(scaSlopeRNAvsDNA, length(vecRNADepth))
        if(!is.null(lsvecidxBatchRNA)){
            lsvecBatchFactorsRNA <- list()
            for(i in seq(1, length(lsvecidxBatchRNA))){
                scaNBatchFactors <- max(lsvecidxBatchRNA[[i]])-1 # Batches are counted from 1
                # Factor of first batch is one (constant), the remaining
                # factors scale based on the first batch.
                vecBatchFactors <- c(1, exp(lsFitRNA$par[(scaNParamUsed+1):(scaNParamUsed+scaNBatchFactors)]))
                scaNParamUsed <- scaNParamUsed+scaNBatchFactors
                # Catch boundary of likelihood domain on batch factor space:
                vecBatchFactors[vecBatchFactors < 10^(-10)] <- 10^(-10)
                vecBatchFactors[vecBatchFactors > 10^(10)] <- 10^(10)
                lsvecBatchFactorsRNA[[i]] <- vecBatchFactors
                
                vecRNAModelFit <- vecRNAModelFit*vecBatchFactors[lsvecidxBatchRNA[[i]]]
            }
        } else { 
            lsvecBatchFactorsRNA <- NULL
        }
        vecRNAModelFitCtrl <- rep(1, length(vecRNADepth))
        if(boolBaselineCtrl) {
            scaSlopeRNAvsDNACtrl <- exp(lsFitRNA$par[scaNParamUsed+1])
            if(scaSlopeRNAvsDNACtrl < 10^(-10)) scaSlopeRNAvsDNACtrl <- 10^(-10) 
            if(scaSlopeRNAvsDNACtrl > 10^(10)) scaSlopeRNAvsDNACtrl <- 10^(10) 
            scaNParamUsed <- scaNParamUsed + 1
            vecRNAModelFitCtrl <- vecRNAModelFitCtrl*scaSlopeRNAvsDNACtrl
        } else { 
            scaSlopeRNAvsDNACtrl <- NULL 
        }
        if(!is.null(lsvecidxBatchRNACtrl)){
            lsvecBatchFactorsRNACtrl <- list()
            for(i in seq(1, length(lsvecidxBatchRNACtrl))){
                scaNBatchFactors <- max(lsvecidxBatchRNACtrl[[i]])-1 # Batches are counted from 1
                # Factor of first batch is one (constant), the remaining
                # factors scale based on the first batch.
                vecBatchFactors <- c(1, exp(lsFitRNA$par[(scaNParamUsed+1):(scaNParamUsed+scaNBatchFactors)]))
                scaNParamUsed <- scaNParamUsed+scaNBatchFactors
                # Catch boundary of likelihood domain on batch factor space:
                vecBatchFactors[vecBatchFactors < 10^(-10)] <- 10^(-10)
                vecBatchFactors[vecBatchFactors > 10^(10)] <- 10^(10)
                lsvecBatchFactorsRNACtrl[[i]] <- vecBatchFactors
                
                vecRNAModelFitCtrl <- vecRNAModelFitCtrl*vecBatchFactors[lsvecidxBatchRNACtrl[[i]]]
            }
        } else { 
            lsvecBatchFactorsRNACtrl <- NULL 
        }
        
        vecFitRNA <- lsModelsDNA[[1]]$vecFitDNAHat*vecRNAModelFit*vecRNADepth
        lsModelRNA <- list(
            vecExprModel=scaSlopeRNAvsDNA,
            vecExprModelCtrl=scaSlopeRNAvsDNACtrl,
            lsvecBatchFactorsRNA=lsvecBatchFactorsRNA,
            lsvecBatchFactorsRNACtrl=lsvecBatchFactorsRNACtrl,
            vecFitRNA=vecFitRNA,
            vecRNAModelFit=vecRNAModelFit,
            vecRNAModelFitCtrl=vecRNAModelFitCtrl,
            scaLL=lsFitRNA$value,
            scaConvergence=lsFitRNA$convergence)
        
        scaLLNew <- evalLogLikDNARNA_gammaDNApoisRNA_direct(
            matDNACounts=matDNACounts,
            matRNACounts=matRNACounts,
            vecDNADepth=vecDNADepth,
            vecRNADepth=vecRNADepth,
            lsvecDNAModel=lapply(lsModelsDNA, function(f) f$vecDNAModel ),
            lsvecCumulBatchFacDNA=lapply(lsModelsDNA, function(f) f$vecCumulBatchFacDNA),
            vecRNAModelFit=lsModelRNA$vecRNAModelFit,
            vecRNAModelFitCtrl=lsModelRNA$vecRNAModelFitCtrl)
        print(paste0(scaIter, ". iteration RNA done with LL=", round(scaLLNew,5)))
    }
    print("Coordinate ascent done.")
    
    scaDF <- sum(sapply(lsModelsDNA, function(f) length(f$par) )) + length(lsFitRNA$par)
    
    return(list( vecDNAModel=lsModelsDNA[[1]]$vecDNAModel,
                 lsvecBatchFacDNA=lsModelsDNA[[1]]$lsvecBatchFacDNA,
                 vecExprModel=lsModelRNA$vecExprModel,
                 vecExprModelCtrl=lsModelRNA$vecExprModelCtrl,
                 lsvecBatchFactorsRNA=lsModelRNA$lsvecBatchFactorsRNA,
                 lsvecBatchFactorsRNACtrl=lsModelRNA$lsvecBatchFactorsRNACtrl,
                 vecFitDNA=vecFitDNA,
                 vecFitRNA=vecFitRNA,
                 scaDisp=lsModelsDNA[[1]]$vecDNAModel,
                 scaDF=scaDF,
                 scaLL=scaLLNew,
                 scaConvergence=lsFitRNA$convergence))
}

#' Fits models to an MPRA dataset
#'  
#' @param obj (object class obj)
#'    Object to be fit.
#' @param vecModelFac (vector of strings number of discrete model variables)
#'    
#' @author David Sebastian Fischer
fitModels <- function(
    obj, vecModelFacRNA, vecModelFacRNACtrl, vecModelFacDNA,
    boolFitDNA = FALSE,
    MAXIT=1000, RELTOL=10^(-8), boolVerbose=TRUE ){
    
    scaNGenes <- dim(obj@matRNACountsProc)[1]
    # Get batch assignments of samples
    if(any(vecModelFacRNA != "1")){
        vecModelFacBatch <- vecModelFacRNA[vecModelFacRNA != "1"]
        lsvecBatchesRNA <- lapply(vecModelFacBatch, function(m){
            vecBatches <- obj@dfAnnotationProc[,m]
            names(vecBatches) <- obj@dfAnnotationProc$Sample
            return(vecBatches)
        })
        names(lsvecBatchesRNA) <- vecModelFacBatch
    } else { 
        lsvecBatchesRNA <- NULL 
    }
    
    # Get batch assignments
    if(length(lsvecBatchesRNA)>0){
        lsvecidxBatchRNA <- list()
        for(batchfactor in lsvecBatchesRNA){
            lsvecidxBatchRNA[[length(lsvecidxBatchRNA)+1]] <- match(batchfactor, unique(batchfactor))
        }
        names(lsvecidxBatchRNA) <- names(lsvecBatchesRNA)
    } else { 
        lsvecidxBatchRNA <- NULL
    }
    
    # Control batches
    if(any(vecModelFacRNACtrl != "1")){
        vecModelFacRNACtrlBatch <- vecModelFacRNACtrl[vecModelFacRNACtrl != "1"]
        lsvecBatchesRNACtrl <- lapply(vecModelFacRNACtrlBatch, function(m){
            vecBatches <- obj@dfAnnotationProc[,m]
            names(vecBatches) <- obj@dfAnnotationProc$Sample
            return(vecBatches)
        })
        names(lsvecBatchesRNACtrl) <- vecModelFacRNACtrlBatch
    } else { 
        lsvecBatchesRNACtrl <- NULL 
    }
    
    # Get batch assignments
    if(length(lsvecBatchesRNACtrl)>0){
        lsvecidxBatchRNACtrl <- list()
        for(batchfactor in lsvecBatchesRNACtrl){
            lsvecidxBatchRNACtrl[[length(lsvecidxBatchRNACtrl)+1]] <- match(batchfactor, unique(batchfactor))
        }
        names(lsvecidxBatchRNACtrl) <- names(lsvecBatchesRNACtrl)
    } else { 
        lsvecidxBatchRNACtrl <- NULL
    }
    
    # DNA batches 
    # Get batch assignments of samples
    if(any(vecModelFacDNA != "1")){
        vecModelFacDNABatch <- vecModelFacDNA[vecModelFacDNA != "1"]
        lsvecBatchesDNA <- lapply(vecModelFacDNABatch, function(m){
            vecBatches <- obj@dfAnnotationProc[,m]
            names(vecBatches) <- obj@dfAnnotationProc$Sample
            return(vecBatches)
        })
        names(lsvecBatchesDNA) <- vecModelFacDNABatch
    } else { 
        lsvecBatchesDNA <- NULL 
    }
    
    # Get batch assignments
    if(length(lsvecBatchesDNA)>0){
        lsvecidxBatchDNA <- list()
        for(batchfactor in lsvecBatchesDNA){
            lsvecidxBatchDNA[[length(lsvecidxBatchDNA)+1]] <- match(batchfactor, unique(batchfactor))
        }
        names(lsvecidxBatchDNA) <- names(lsvecBatchesDNA)
    } else { 
        lsvecidxBatchDNA <- NULL
    }
    
    # Find control genes
    if(is.null(obj@vecCtrlIDs)){
        vecidxCtrl <- NULL
        vecidxCase <- seq(1, dim(obj@matDNACountsProc)[1])
    } else {
        vecidxAllIDs <- seq(1, dim(obj@matDNACountsProc)[1])
        vecidxCtrl <- match(obj@vecCtrlIDs, rownames(obj@matDNACountsProc))
        vecidxCase <- setdiff(vecidxAllIDs, vecidxCtrl)
    }
    
    if(obj@strModel=="pointlnDNAnbRNA"){
        # First fit DNA distribution
        # Then condition RNA distribution on DNA MLE
        # These are handled as separate calls to this function.
        if(boolFitDNA) {
            if(boolVerbose) message(paste0("Fit DNA"))
            lsFits <- bplapply(seq(1,scaNGenes),function(i){
                if(boolVerbose) message(paste0("### DNA Enhancer ", i))
                fitDNA_lnDNA(
                    vecDNACounts=obj@matDNACountsProc[i,,drop=FALSE],
                    vecDNADepth=obj@vecDNADepth,
                    lsvecidxBatchDNA=lsvecidxBatchDNA,
                    MAXIT=MAXIT,
                    RELTOL=RELTOL )
            })
            names(lsFits) <- rownames(obj@matRNACountsProc)
        } else {
            if(is.null(obj@lsDNAModelFits)) {
                stop(print("ERROR: Fit DNA models first."))
            }
            if(boolVerbose) message(paste0("Fit RNA"))
            lsFits <- bplapply(vecidxCase,function(i){
                if(boolVerbose) message(paste0("### RNA Enhancer ", i))
                lsFitRNA <- fitRNA_pointDNAnbRNA(
                    matRNACounts=obj@matRNACountsProc[c(i,vecidxCtrl),,drop=FALSE],
                    matDNAEst=do.call(rbind, lapply(c(i, vecidxCtrl), function(j) 
                        obj@lsDNAModelFits[[j]]$vecFitDNAHat )),
                    scaDisp=obj@vecDispersions[i],
                    vecRNADepth=obj@vecRNADepth,
                    lsvecidxBatchRNA=lsvecidxBatchRNA,
                    boolBaselineCtrl="1" %in% vecModelFacRNACtrl,
                    lsvecidxBatchRNACtrl=lsvecidxBatchRNACtrl,
                    MAXIT=MAXIT,
                    RELTOL=RELTOL )
                # Summarise two fits
                lsFit <- list(
                    vecDNAModel=c(obj@lsDNAModelFits[[i]]$scaMuDNA, 
                                  obj@lsDNAModelFits[[i]]$scaSdDNA),
                    scaRNAModel=lsFitRNA$vecExprModel,
                    scaRNAModelCtrl=lsFitRNA$vecExprModelCtrl,
                    lsvecBatchFacDNA=obj@lsDNAModelFits[[i]]$lsvecBatchFactorsDNA,
                    lsvecBatchFacRNA=lsFitRNA$lsvecBatchFactorsRNA,
                    lsvecBatchFactorsRNACtrl=lsFitRNA$lsvecBatchFactorsRNACtrl,
                    vecFitDNA=obj@lsDNAModelFits[[i]]$vecFitDNA,
                    vecFitRNA=lsFitRNA$vecFitRNA,
                    scaDisp=lsFitRNA$scaDisp,
                    scaDF=lsFitRNA$scaDF,
                    scaLL=lsFitRNA$scaLL,
                    scaConvergence=lsFitRNA$scaConvergence)
                return(lsFit)
            })
            names(lsFits) <- rownames(obj@matRNACountsProc)[vecidxCase]
        }
    } else if(obj@strModel=="gammaDNApoisRNA"){
        lsFits <- bplapply(vecidxCase, function(i){
            if(boolVerbose) message(paste0("### Enhancer ", i))
            fitDNARNA_gammaDNApoisRNA(
                matDNACounts=obj@matDNACountsProc[c(i,vecidxCtrl),,drop=FALSE],
                matRNACounts=obj@matRNACountsProc[c(i,vecidxCtrl),,drop=FALSE],
                vecDNADepth=obj@vecDNADepth,
                vecRNADepth=obj@vecRNADepth,
                lsvecidxBatchDNA=lsvecidxBatchDNA,
                lsvecidxBatchRNA=lsvecidxBatchRNA,
                boolBaselineCtrl="1" %in% vecModelFacRNACtrl,
                lsvecidxBatchRNACtrl=lsvecidxBatchRNACtrl,
                lsDNAModelFitsCtrl=obj@lsDNAModelFitsCtrl,
                MAXIT=MAXIT,
                RELTOL=RELTOL )
        })
        names(lsFits) <- rownames(obj@matRNACountsProc)[vecidxCase]
    } else if(obj@strModel=="gammaDNApoisRNA_coordascent"){
        lsFits <- bplapply(vecidxCase, function(i){
            if(boolVerbose) message(paste0("### Enhancer ", i))
            fitDNARNA_gammaDNApoisRNA_coordascent(
                matDNACounts=obj@matDNACountsProc[c(i,vecidxCtrl),,drop=FALSE],
                matRNACounts=obj@matRNACountsProc[c(i,vecidxCtrl),,drop=FALSE],
                vecDNADepth=obj@vecDNADepth,
                vecRNADepth=obj@vecRNADepth,
                lsvecidxBatchDNA=lsvecidxBatchDNA,
                lsvecidxBatchRNA=lsvecidxBatchRNA,
                boolBaselineCtrl="1" %in% vecModelFacRNACtrl,
                lsvecidxBatchRNACtrl=lsvecidxBatchRNACtrl,
                MAXIT=MAXIT,
                RELTOL=RELTOL )
        })
        names(lsFits) <- rownames(obj@matRNACountsProc)[vecidxCase]
    } else {
        stop(paste0("obj@strModel not recognised in fitModels(): ",
                    obj@strModel))
    }
    
    return(lsFits)
}