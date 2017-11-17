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
simulateMPRA <- function(n.case=100, n.ctrl=20,
                          n.cond=2, n.bc=25, n.reps=2,
                          mu.dna=100, sd.dna.cond=2, sd.dna.bc=0.1,
                          mu.rna=1, sd.rna=0.2, sd.rna.cond=1) {
    
    n.samples <- n.cond*n.bc*n.reps
    n.enhancers <- n.case+n.ctrl
    # sample allocation to condition
    idx.cond <- rep(sapply(seq(1, n.cond), function(x) rep(x, n.bc) ), n.reps)
    idx.rep <- sapply(seq(1, n.reps), function(x) rep(x, n.bc*n.cond) )
    
    ## draw dna counts
    fc.bc <- rep(rep(rnorm(n = n.bc, mean = 1, sd = sd.dna.bc), n.cond), n.reps)
    fc.bc[fc.bc] <- 0.2 # Threshold
    dcounts.case <- do.call(rbind, lapply(
        sample(x=mu.dna, size=n.case), 
        function(mu) {
            round(rlnorm(n=n.samples, meanlog = log(mu), sdlog = sd.dna.cond)*fc.bc)
        }))
    rownames(dcounts.case) <- paste0("case_", seq_len(n.case))
    if(n.ctrl > 0){
        dcounts.ctrl <- do.call(rbind, lapply(
            sample(x=mu.dna, size=n.ctrl), function(mu) {
                round(rlnorm(n=n.samples, meanlog = log(mu), sdlog = sd.dna.cond)*fc.bc)
            }))
        rownames(dcounts.ctrl) <- paste0("ctrl_", seq_len(n.ctrl))
    } else {
        dcounts.ctrl <- NULL
    }
    dcounts <- rbind(dcounts.case, dcounts.ctrl)
    colnames(dcounts) <- paste0("S", seq_len(n.samples))
    dcounts[dcounts <= 1] <- 1
    
    ## draw rna counts
    fc.rna <- rnorm(n=n.case, mean=mu.rna, sd=sd.rna)
    fc.rna[fc.rna < 0.1] <- 0.1 # Threshold
    fc.rna.cond <- matrix(1, nrow=n.case, ncol=n.cond)
    fc.rna.cond[,-1] <- rnorm(n=n.case*(n.cond-1), 
                              mean=mu.rna, sd=sd.rna)
    fc.rna.cond[fc.rna.cond < 0.1] <- 0.1 # Threshold
    rcounts <- do.call(rbind, lapply(seq_len(n.case), function(i) {
        sapply(dcounts[i,]*fc.rna.cond[i,idx.cond], function(x) 
            round(rpois(n=1, lambda = x*fc.rna[i])) )
    }))
    if(n.ctrl > 0){
        fc.rna.ctrl <- rnorm(n=1, mean=mu.rna, sd=sd.rna)
        fc.rna.ctrl[fc.rna.ctrl < 0.1] <- 0.1 # Threshold
        rcounts.ctrl <- do.call(rbind, lapply(n.case+seq_len(n.ctrl), function(i) {
            sapply(dcounts[i,]*fc.rna.ctrl, function(x) 
                round(rpois(n=1, lambda=x )) )
        }))
        rcounts <- rbind(rcounts, rcounts.ctrl)
    }
    rownames(rcounts) <- rownames(dcounts)
    colnames(rcounts) <- colnames(dcounts)
    
    ## create annotation
    colAnnot <- data.frame(
        sample=colnames(dcounts),
        cond=paste0("cond_", seq(1, n.cond))[idx.cond],
        rep=paste0("rep_", seq(1, n.reps))[idx.rep],
        barcode=rep(paste0("BC", seq(1, n.bc)), n.cond*n.reps),
        stringsAsFactors = FALSE
    )
    colAnnot$experiment <- colAnnot$cond # for depth factors
    rowAnnot <- data.frame(
        type=c(rep("case", n.case),
               rep("ctrl", n.ctrl)),
        cond2_effect=1,
        stringsAsFactors = FALSE
    )
    if(n.cond > 1) {
        rowAnnot[rowAnnot$type == "case",]$cond2_effect <- fc.rna.cond[,2]
    }
    
    return(list( dna=dcounts,
                 rna=rcounts,
                 colAnnot=colAnnot,
                 rowAnnot=rowAnnot))
}