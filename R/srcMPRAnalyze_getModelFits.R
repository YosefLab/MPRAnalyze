### Get model fits

#' Return model fits
#'  
#' @param obj (MPRAnalyzeObject) Object with fits.
#' 
#' @return (matrix enhancers x samples) Matrix with DNA fit values.
#' 
#' @author David Sebastian Fischer
#' 
#' @export
getDNAFitFull <- function( obj ){
   
    matModelFitsDNA <- do.call(rbind, lapply(objMPRA@lsModelFitsFull, function(fit) {
        fit$vecFitDNA
    }))
    return(matModelFitsDNA)
}
getDNAFitRed <- function( obj ){
   
    matModelFitsDNA <- do.call(rbind, lapply(objMPRA@lsModelFitsRed, function(fit) {
        fit$vecFitDNA
    }))
    return(matModelFitsDNA)
}
getRNAFitFull <- function( obj ){
   
    matModelFitsRNA <- do.call(rbind, lapply(objMPRA@lsModelFitsFull, function(fit) {
        fit$vecFitRNA
    }))
    return(matModelFitsRNA)
}
getRNAFitRed <- function( obj ){
   
    matModelFitsRNA <- do.call(rbind, lapply(objMPRA@lsModelFitsRed, function(fit) {
        fit$vecFitRNA
    }))
    return(matModelFitsRNA)
}

getDNAModel <- function( obj ){
   
    matModelDNA <- do.call(rbind, lapply(objMPRA@lsModelFitsFull, function(fit) {
        fit$vecExprModel
    }))
    return(matModelDNA)
}
getExprModelFull <- function( obj ){
   
    matModelExpr <- do.call(rbind, lapply(objMPRA@lsModelFitsFull, function(fit) {
        fit$vecExprModel
    }))
    matModelExprCtrl <- do.call(rbind, lapply(objMPRA@lsModelFitsFull, function(fit) {
        fit$vecExprModelCtrl
    }))
    return(list(
        baseline=matModelExpr,
        control=matModelExprCtrl
    ))
}
getExprModelRed <- function( obj ){
   
    matModelExpr <- do.call(rbind, lapply(objMPRA@lsModelFitsRed, function(fit) {
        fit$vecExprModel
    }))
    matModelExprCtrl <- do.call(rbind, lapply(objMPRA@lsModelFitsRed, function(fit) {
        fit$vecExprModelCtrl
    }))
    return(list(
        baseline=matModelExpr,
        control=matModelExprCtrl
    ))
}