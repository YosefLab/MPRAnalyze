#' Simulate MPRA data
#' 
#' @param n.case number of case enhancers
#' @param n.ctrl number of control/scrambled enhancers
#' @param n.cond number of conditions
#' @param n.bc number of barcodes observed per condition
#' @param mu.dna mean dna level
#' @param sd.dna.cond stddev to draw dna fold change 
#' by condition from (centred at 1)
#' @param sd.dna.bc stddev to draw dna fold change 
#' by barcode from (centred at 1)
#' @param mu.rna mean to draw dna fold change from
#' @param sd.rna stddev to draw rna slope from (centred at mu.rna)
#' @param sd.rna.cond stddev to draw rna slope fold change 
#' by condition from (centred at 1)
#' 
#' @return (list)
#' \itemize{
#' \item dna dna count matrix
#' \item rna rna count matrix
#' \item colData column annotation of samples (meta data)
#' \item rowAnnot row annotation of enhancers (meta data)
#' }
#' 
#' @export
simulate.mpra <- function(n.case=100, n.ctrl=20,
                          n.cond=2, n.bc=25,
                          mu.dna=100, sd.dna.cond=2, sd.dna.bc=0.1,
                          mu.rna=1, sd.rna=0.2, sd.rna.cond=1) {
    
    n.samples <- n.cond*n.bc
    n.enhancers <- n.case+n.ctrl
    # sample allocation to condition
    idx.cond <- sapply(seq(1, n.cond), function(x) rep(x, n.bc) )
    
    ## draw dna counts
    fc.bc <- rep(rnorm(n = n.bc, mean = 1, sd = sd.dna.bc), n.cond)
    fc.bc[fc.bc] <- 0.2 # Threshold
    dcounts.case <- do.call(rbind, lapply(
        sample(x=mu.dna, size=n.case), 
        function(mu) {
            round(rlnorm(n=n.samples, meanlog = log(mu), sdlog = sd.dna.cond)*fc.bc)
        }))
    rownames(dcounts.case) <- paste0("case_", seq_len(n.case))
    dcounts.ctrl <- do.call(rbind, lapply(
        sample(x=mu.dna, size=n.ctrl), function(mu) {
            round(rlnorm(n=n.samples, meanlog = log(mu), sdlog = sd.dna.cond)*fc.bc)
        }))
    rownames(dcounts.ctrl) <- paste0("ctrl_", seq_len(n.ctrl))
    dcounts <- rbind(dcounts.case, dcounts.ctrl)
    colnames(dcounts) <- paste0("S", seq_len(n.samples))
    dcounts[dcounts <= 1] <- 1
    
    ## draw rna counts
    fc.rna <- rnorm(n=n.enhancers, mean=mu.rna, sd=sd.rna)
    fc.rna[fc.rna < 0.1] <- 0.1 # Threshold
    fc.rna.cond <- matrix(1, nrow=n.enhancers, ncol=n.cond)
    fc.rna.cond[,-1] <- rnorm(n=n.enhancers*(n.cond-1), 
                              mean=mu.rna, sd=sd.rna)
    fc.rna.cond[fc.rna.cond < 0.1] <- 0.1 # Threshold
    rcounts <- do.call(rbind, lapply(seq_len(n.enhancers), function(i) {
        sapply(dcounts[i,]*fc.rna.cond[i,idx.cond], function(x) 
            round(rpois(n=1, lambda = x*fc.rna[i])) )
    }))
    rownames(dcounts) <- rownames(dcounts)
    colnames(rcounts) <- colnames(rcounts)
    
    ## create annotation
    colAnnot <- data.frame(
        sample=colnames(dcounts),
        cond=paste0("cond_", seq(1, n.cond))[idx.cond],
        barcode=rep(paste0("BC", seq(1, n.bc)), n.cond),
        stringsAsFactors = FALSE
    )
    colAnnot$experiment <- colAnnot$cond # for depth factors
    rowAnnot <- data.frame(
        type=c(rep("case", n.case),
               rep("ctrl", n.ctrl)),
        stringsAsFactors = FALSE
    )
    
    return(list( dna=dcounts,
                 rna=rcounts,
                 colAnnot=colAnnot,
                 rowAnnot=rowAnnot))
}