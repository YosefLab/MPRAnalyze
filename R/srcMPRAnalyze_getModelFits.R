### Get model fits

#' Return models or model fits.
#'  
#' @param obj (MPRAnalyzeObject) Object with fits.
#' 
#' @return (matrix enhancers x samples) Matrix with models or model fits.
#' 
#' @name summarizeModel
#' @rdname summarizeModel
#' @aliases 
#' getDNAFitFull
#' getDNAFitRed 
#' getRNAFitFull
#' getRNAFitRed 
#' getDNAModel
#' getExprModelFull 
#' getExprModelRed
#' 
#' @author David Sebastian Fischer
NULL

#' @rdname summarizeModel
#' @export
getDNAFitFull <- function( obj ){
   
    matModelFitsDNA <- do.call(rbind, lapply(objMPRA@lsModelFitsFull, function(fit) {
        fit$vecFitDNA
    }))
    return(matModelFitsDNA)
}

#' @rdname summarizeModel
#' @export
getDNAFitRed <- function( obj ){
   
    matModelFitsDNA <- do.call(rbind, lapply(objMPRA@lsModelFitsRed, function(fit) {
        fit$vecFitDNA
    }))
    return(matModelFitsDNA)
}

#' @rdname summarizeModel
#' @export
getRNAFitFull <- function( obj ){
   
    matModelFitsRNA <- do.call(rbind, lapply(objMPRA@lsModelFitsFull, function(fit) {
        fit$vecFitRNA
    }))
    return(matModelFitsRNA)
}

#' @rdname summarizeModel
#' @export
getRNAFitRed <- function( obj ){
   
    matModelFitsRNA <- do.call(rbind, lapply(objMPRA@lsModelFitsRed, function(fit) {
        fit$vecFitRNA
    }))
    return(matModelFitsRNA)
}

#' @rdname summarizeModel
#' @export
getDNAModel <- function( obj ){
   
    matModelDNA <- do.call(rbind, lapply(objMPRA@lsModelFitsFull, function(fit) {
        fit$vecExprModel
    }))
    return(matModelDNA)
}

#' @rdname summarizeModel
#' @export
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

#' @rdname summarizeModel
#' @export
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