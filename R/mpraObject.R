setClass("MpraObject", slots = c(
    dnaCounts = "matrix",
    rnaCounts = "matrix",
    dnaDesign = "Matrix", ##sparse
    rnaDesign = "Matrix", ##sparse
    controls = "integer", #idx of negative controls

    dneDepth = "numeric",
    rnaDepth = "numeric",

    model = "character",
    dnaDesign = "formula",
    dnaDesign.mat = "matrix",
    rnaDesign = "formula",
    rnaDesign.mat = "matrix"

))