
#' evaluate the performance of the supported models on the given data
#'
#' @param obj the MpraObject
evalModels <- function(obj) {
    ##TODO: use controls\scrambles if available, otherwise sample ~200 random
    ##rows to evaluate
}

#' plot various diagnostic plots to inform the user on model performance and
#' inform manual model selection
#'
#' @export
#'
#' @param obj the MpraObject
plotModelPerformance <- function(obj) {
    modelDiag <- evalModels()

    ##TODO
}


autoChooseModel <- function(obj) {

    modelDiag <- evalModels()

    ##TODO
}
