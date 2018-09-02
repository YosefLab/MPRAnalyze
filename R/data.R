#'Column annotations for the chromosomal-episomal dataset
#'@format
#'\describe{
#'\item{batch}{batch identifier, factor}
#'\item{condition}{condition identifier, factor. WT corresponds to chromosomal
#'and MT corresponds to episomal}
#'\item{barcode}{barcode identifier, factor}
#'}
#'@source \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5204343/}
"colAnnot"

#' DNA observations for the chromosomal-episomal dataset
#' @format integer matrix
#' @source \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5204343/}
"dnaCounts"

#' DNA observations for the chromosomal-episomal dataset
#' @format integer matrix
#' @source \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5204343/}
"rnaCounts"

#' indices of control enhancers in the chromosomal-episomal dataset
#' @format logical vector, TRUE means the corresponding row in the observation
#' matrices are control features
#' @source #' DNA observations for the chromosomal-episomal dataset
#' @format integer matrix
#' @source \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5204343/}
"control"
