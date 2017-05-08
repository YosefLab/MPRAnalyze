extractModel <- function(
      strModelFull=NULL,
      strModelRed=NULL,
      vecModelFacRNAFull=NULL,
      vecModelFacRNAFullCtrl=NULL,
      vecModelFacRNARed=NULL,
      vecModelFacRNARedCtrl=NULL){
    
    if(is.null(strModelFull) & is.null(vecModelFacRNAFull)){
        stop(paste0("Supply either strModelFull or vecModelFacRNAFull."))
    }
    if(is.null(strModelRed) & is.null(vecModelFacRNARed)){
        stop(paste0("Supply either strModelRed or vecModelFacRNARed."))
    }
    if(!is.null(strModelFull) & !is.null(vecModelFacRNAFull)){
        warning(paste0("Using strModelFull and not vecModelFacRNAFull."))
    }
    if(!is.null(strModelRed) & !is.null(vecModelFacRNARed)){
        warning(paste0("Using strModelRed and not vecModelFacRNARed."))
    }
    
    lsModelFac <- list()
    if(is.null(strModelFull)){
        lsModelFac$vecModelFacRNAFull <- vecModelFacRNAFull
        lsModelFac$vecModelFacRNAFullCtrl <- vecModelFacRNAFullCtrl
    } else {
        stop("building")
        vecModelFacRNAFull <- c("1")
        vecModelFacRNAFullCtrl <- c()
        
        vecElModelFull <- unlist(strsplit(strModelFull, split = "[+]"))
        vecboolElUsed <- rep(FALSE, length(vecElModelFull))
        vecIdxCtrlInteract <- grep(pattern="CtrlSequence:", x=vecElModelFull)

        if(any(vecElModelFull=="1")) {
            vecboolElUsed[which(vecElModelFull=="1")] <- TRUE
        }
        if(any(vecElModelFull=="CtrlSequence")) {
            vecModelFacRNAFullCtrl <- c(vecModelFacRNAFullCtrl, "1")
            vecboolElUsed[which(vecElModelFull=="CtrlSequence")] <- TRUE
        }
        
    }
    if(is.null(strModelRed)){
        lsModelFac$vecModelFacRNARed <- vecModelFacRNARed
        lsModelFac$vecModelFacRNARedCtrl <- vecModelFacRNARedCtrl
    } else {
        stop("building")
        vecModelFacRNARed <- c("1")
        vecModelFacRNARedCtrl <- c()
        
        vecElModelFull <- unlist(strsplit(strModelFull, split = "[+]"))
        vecboolElUsed <- rep(FALSE, length(vecElModelFull))
        vecIdxCtrlInteract <- grep(pattern="CtrlSequence:", x=vecElModelFull)

        if(any(vecElModelFull=="1")) {
            vecboolElUsed[which(vecElModelFull=="1")] <- TRUE
        }
        if(any(vecElModelFull=="CtrlSequence")) {
            vecModelFacRNAFullCtrl <- c(vecModelFacRNAFullCtrl, "1")
            vecboolElUsed[which(vecElModelFull=="CtrlSequence")] <- TRUE
        }
    }
    
    return(lsModelFac)
}