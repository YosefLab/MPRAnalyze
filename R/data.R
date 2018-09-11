#' Sample MPRA data
#' 
#' A subset of MPRA data from Inoue et al., comparing enhancer activity of
#' episomal constructs vs. chromosomally integrated constructs (integration
#' was performed with a lentivirus). Data included negative control enhancers,
#' multiple batches and barcodes, a subsample of which are included in this
#' sample data for runtime purposes.
#'
#'@docType data
#'
#'@usage data(ChrEpi)
#'
#'@format
#'\describe{
#'\item{ce.colAnnot}{
#'Column annotations for each column (sample) in the data matrices
#'\describe{
#'\item{batch}{batch identifier, factor}
#'\item{condition}{condition identifier, factor. WT corresponds to chromosomal
#'and MT corresponds to episomal}
#'\item{barcode}{barcode identifier, factor}
#'}
#'}
#'\item{ce.dnaCounts}{DNA observations}
#'\item{ce.rnaCounts}{DNA observations}
#'\item{ce.control}{indices of control enhancers}
#'}
#'
#'@source \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5204343/}
#'
#'@name ChrEpi
NULL

#'@name ChrEpi
"ce.colAnnot"

#'@name ChrEpi
"ce.dnaCounts"

#'@name ChrEpi
"ce.rnaCounts"

#'@rdname ChrEpi
"ce.control"
