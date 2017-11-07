setClass("MPRAobject", slots = c(
    dnaCounts = "matrix",
    rnaCounts = "matrix",
    colAnnot = "data.frame",
    controls = "integer",
    design = "formula"
))

MPRAobject <- function(dnaCounts, rnaCounts, colAnnot, controls, design) {
    return(new("MPRAobject", dnaCounts=dnaCounts, rnaCounts=rnaCounts,
               colAnnot=colAnnot, controls=controls, design = design))
}

##' one condition scenario:
##' - normalize for library size
##' - choose model (?)
##' - fit model (single model!)
##' - extract model params
##' - visualizations
##'
##' Two conditions scenario:
##' - normalize for library size
##' - choose model
##' - fit model (full + reduced!)
##' - differential activity analysis
##' - Visualizations

getDesign <- function(obj, formula) {
    design <- sparse.model.matrix(formula, data = obj@colAnnot)
    # model.terms <- attr(terms(formula), "term.labels")
    # design <- do.call(cbind, lapply(model.terms, function(term) {
    #     sparse.model.matrix(object = as.formula(c("~ 0 +", term)),
    #                         data = obj@colAnnot)
    # }))
    return(design)
}

fitSingleModel <- function(obj) {}

fitComparativeModel <- function(obk, condition) {}

