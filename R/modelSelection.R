
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

#' Set the distributional model used
#' @param obj the MPRAnalyze object
#' @param model the charater identifier of the model to be used. Currently 
#' supported models: "ln.nb", "gamma.pois".
#' @return the MPRAnalyze with the model set for the given value
#' @export
setModel <- function(obj, model) {
    PSEUDOCOUNT <- 10^(-5)
    if(model %in% c("ln.nb", "gamma.pois", "ln.ln")) {
        obj@model <- model
    } else {
        stop("model '", model, "' not supported.")
    }
    if(model == "ln.nb") {
        obj@dnaCounts[obj@dnaCounts==0] <- PSEUDOCOUNT
    }
    if(model == "ln.ln") {
        obj@dnaCounts[obj@dnaCounts==0] <- PSEUDOCOUNT
        obj@rnaCounts[obj@rnaCounts==0] <- PSEUDOCOUNT
    }
    return(obj)
}

autoChooseModel <- function(obj) {
    obj@model <- "gamma.pois"
    ##modelDiag <- evalModels()
    ##TODO: select the best model somehow
    return(obj)
}
