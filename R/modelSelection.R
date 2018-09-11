
#' Set the distributional model used. Default is gamma.pois, and is recommended.
#' Other supoprted models are ln.nb in which the DNA follows a log-normal
#' distribution and the RNA follows a negative binomial, and ln.ln in which
#' both follow log-normal distributions.
#' To use alternative distributional models, use this function before fitting
#' the model.
#' @param obj the MPRAnalyze object
#' @param model the charater identifier of the model to be used. Currently 
#' supported models: "ln.nb", "gamma.pois", "ln.ln"
#' @return the MPRAnalyze with the model set for the given value
#' @export
#' 
#' @examples
#' data <- simulateMPRA(tr = rep(2,10), da=NULL, nbatch=2, nbc=20)
#' obj <- MpraObject(dnaCounts = data$obs.dna, 
#'                   rnaCounts = data$obs.rna, 
#'                   colAnnot = data$annot)
#' obj <- estimateDepthFactors(obj, lib.factor = "batch", which.lib = "both")
#' obj <- setModel(obj, "ln.ln")
#' obj <- analyzeQuantification(obj, dnaDesign = ~ batch + barcode, 
#'                               rnaDesign = ~1)
setModel <- function(obj, model) {
    
    if(model %in% c("ln.nb", "gamma.pois", "ln.ln")) {
        obj@model <- model
    } else {
        stop("model '", model, "' not supported.")
    }
    
    # check that data is integer if count distribution is used
    # round if this is not the case and issue a warning
    if(model(obj) %in% c("ln.nb", "gamma.pois")){
        if(!all(rnaCounts(obj) %% 1 == 0)) {
            warning("Negative binomial count distribution of rna model in ",
                    model(obj), " requires integer RNA observations. ",
                    "Rounded RNA count matrix to integers.")
            obj@rnaCounts <- round(rnaCounts(obj))
        }
    }
    return(obj)
}

autoChooseModel <- function(obj) {
    obj@model <- "gamma.pois"
    return(obj)
}
