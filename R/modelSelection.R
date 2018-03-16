
#' Set the distributional model used
#' @param obj the MPRAnalyze object
#' @param model the charater identifier of the model to be used. Currently 
#' supported models: "ln.nb", "gamma.pois".
#' @return the MPRAnalyze with the model set for the given value
#' @export
setModel <- function(obj, model) {
    
    if(model %in% c("ln.nb", "gamma.pois", "ln.ln")) {
        obj@model <- model
    } else {
        stop("model '", model, "' not supported.")
    }
    
    # check that data is integer if count distribution is used
    # round if this is not the case and issue a warning
    if(obj@model %in% c("ln.nb", "gamma.pois")){
        if(!all(obj@rnaCounts%%1 == 0)) {
            warning("Negative binomial count distribution of rna model in ",
                    obj@model, " requires integer RNA observations. ",
                    "Rounded RNA count matrix to integers.")
            obj@rnaCounts <- round(obj@rnaCounts)
        }
    }
    return(obj)
}

autoChooseModel <- function(obj) {
    obj@model <- "gamma.pois"
    ##modelDiag <- evalModels()
    ##TODO: select the best model somehow
    return(obj)
}
