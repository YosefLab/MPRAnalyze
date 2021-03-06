% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/accessFunctions.R
\name{getFits_RNA}
\alias{getFits_RNA}
\title{Get RNA model-based estimates from an MpraObject (the expected values based 
on the model). These can be compared with the observed counts to assess 
goodness of fit.}
\usage{
getFits_RNA(obj, enhancers = NULL, depth = TRUE, full = TRUE, rnascale = TRUE)
}
\arguments{
\item{obj}{MpraObject to extract from}

\item{enhancers}{which enhancers to get the fits for. Can be character 
vectors with enhancer names, logical or numeric enhancer indices, or NULL if 
all enhancers are to be extracted (default)}

\item{depth}{include depth correction in the model fitting (default TRUE)}

\item{full}{if LRT modeling was used, TRUE (default) would return the fits
of the full model, FALSEwould return the reduced model fits.}

\item{rnascale}{if controls were used to correct the fitting (in comparative
analyses), use these factors to re-adjust the estimates back.}
}
\value{
RNA fits (numeric, enhancers x samples)
}
\description{
Get RNA model-based estimates from an MpraObject (the expected values based 
on the model). These can be compared with the observed counts to assess 
goodness of fit.
}
\examples{
data <- simulateMPRA(tr = rep(2,5), nbatch=2, nbc=15)
obj <- MpraObject(dnaCounts = data$obs.dna, 
                  rnaCounts = data$obs.rna, 
                  colAnnot = data$annot)
obj <- estimateDepthFactors(obj, lib.factor = "batch", which.lib = "both")
obj <- analyzeQuantification(obj, dnaDesign = ~ batch + barcode, 
                              rnaDesign = ~1)
rna.fits <- getFits_RNA(obj)
}
