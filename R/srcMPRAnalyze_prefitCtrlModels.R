prefitCtrlDNA_gammaDNApoisRNA_coordascent <- function(
    matDNACounts,
    matRNACounts,
    vecDNADepth,
    vecRNADepth,
    lsvecidxBatchDNA,
    lsvecidxBatchRNA,
    MAXIT=1000,
    RELTOL=10^(-8) ){
    
    # Initialise
    # DNA models
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
    
    lsModelRNA <- list(
        vecExprModel=scaSlopeRNAvsDNA,
        vecExprModelCtrl=NULL,
        lsvecBatchFactorsRNA=lsvecBatchFacRNA,
        lsvecBatchFactorsRNACtrl=NULL,
        vecRNAModelFit=rep(scaSlopeRNAvsDNA, length(vecRNADepth)),
        vecRNAModelFitCtrl=NULL )
    
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
        lsModelsDNA <- bplapply(seq(1, dim(matDNACounts)[1]), function(i){
            vecModel <- lsModelRNA$vecRNAModelFit
            
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
            vecRNAModelFitCtrl=rep(1, length(lsModelRNA$vecRNAModelFit)))
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
                boolBaselineCtrl=FALSE,
                lsvecidxBatchRNACtrl=NULL,
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
        
        vecFitRNA <- lsModelsDNA[[1]]$vecFitDNAHat*vecRNAModelFit*vecRNADepth
        lsModelRNA <- list(
            vecExprModel=scaSlopeRNAvsDNA,
            vecExprModelCtrl=NULL,
            lsvecBatchFactorsRNA=lsvecBatchFactorsRNA,
            lsvecBatchFactorsRNACtrl=NULL,
            vecFitRNA=vecFitRNA,
            vecRNAModelFit=vecRNAModelFit,
            vecRNAModelFitCtrl=NULL,
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
            vecRNAModelFitCtrl=rep(1, length(lsModelRNA$vecRNAModelFit)))
        print(paste0(scaIter, ". iteration RNA done with LL=", round(scaLLNew,5)))
    }
    print("Coordinate ascent done.")
    
    scaDF <- sum(sapply(lsModelsDNA, function(f) length(f$par) )) + length(lsFitRNA$par)
    
    return(list(
        lsModelsDNA=lsModelsDNA,
        lsModelRNA=lsModelRNA
    ))
}

prefitCtrlDNAModels <- function(
    obj, vecModelFacRNA, vecModelFacDNA,
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
    
    # this prefitting only makes sense for gammaDNApoisRNA_coordascent
    lsFits <- prefitCtrlDNA_gammaDNApoisRNA_coordascent(
        matDNACounts=obj@matDNACountsProc[obj@vecCtrlIDs,,drop=FALSE],
        matRNACounts=obj@matRNACountsProc[obj@vecCtrlIDs,,drop=FALSE],
        vecDNADepth=obj@vecDNADepth,
        vecRNADepth=obj@vecRNADepth,
        lsvecidxBatchDNA=lsvecidxBatchDNA,
        lsvecidxBatchRNA=lsvecidxBatchRNA,
        MAXIT=MAXIT,
        RELTOL=RELTOL )
    
    return(lsFits$lsModelsDNA)
}