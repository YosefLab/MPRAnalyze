#' Simulate MPRA data
#' 
#' @import stats
#' @import methods
#' 
#' @param model noise model out of {"gamma.pois", "ln.nb", "ln.ln"}
#' @param n.case number of case enhancers
#' @param n.ctrl number of control/scrambled enhancers
#' @param n.cond number of conditions
#' @param n.bc number of barcodes observed per condition
#' @param n.reps number of repeats to generate
#' @param mu.dna mean dna level
#' @param sd.dna sd to draw DNA from
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
simulateMPRA <- function(
    model="gamma.pois", n.case=10, n.ctrl=10,
    n.cond=2, n.bc=4, n.reps=2, frac.de=0.5,
    mu.dna=100, mu.sd.dna=0.3, sd.dna=10, sd.dna.cond=0.2, sd.dna.bc=0.2,
    mu.rna=1, mu.sd.rna=0.1, sd.rna=0.2, sd.rna.cond=0.1, sd.rna.cond_case=0.1) {
    
    if(!model %in% c("gamma.pois", "ln.nb", "ln.ln")) {
        stop("model not recognized")
    }
    if(model=="gamma.pois" & !is.null(sd.rna)){
        warning("sd.rna is ignored in gamma.pois model")
    }
    if(model=="ln.nb"){
        warning("sd.rna is interpreted as size parameter of negative binomial")
    }
    if(n.cond==1 & (sd.rna.cond!=0 | sd.dna.cond!=0)){
        warning("non-zero sd.dna.cond or sd.rna.cond supplied",
                " with only one condition, setting these to 0")
        sd.dna.cond <- 0
        sd.rna.cond <- 0
    }
    n.samples <- n.cond*n.bc*n.reps
    # sample allocation to condition
    idx.cond <- rep(sapply(seq(1, n.cond), function(x) rep(x, n.bc) ), n.reps)
    idx.rep <- sapply(seq(1, n.reps), function(x) rep(x, n.bc*n.cond) )
    idx.bc <- rep(seq_len(n.bc), n.cond*n.reps)
    
    ## draw dna counts
    par.dna.mu.case <- rlnorm(
        n=n.case, meanlog=log(mu.dna / sqrt(1+mu.sd.dna^2/mu.dna^2)), 
        sdlog=sqrt(log(1+mu.sd.dna^2/mu.dna)))
    par.dna.bc.case <- do.call(rbind, lapply(seq_len(n.case), function(i) {
        exp(rnorm(n = n.bc, mean = 0, sd = sd.dna.bc))[idx.bc]
    }))
    par.dna.cond.case <- do.call(rbind, lapply(seq_len(n.case), function(i) {
        exp(rnorm(n = n.cond, mean = 0, sd = sd.dna.bc))[idx.cond]
    }))
    dcounts.case.true <- do.call(rbind, lapply(seq_len(n.case), function(i) {
        par.dna.mu.case[i]*par.dna.bc.case[i,]*par.dna.cond.case[i,]
    }))
    dcounts.case.obs <- do.call(rbind, lapply(seq_len(n.case), function(i) {
        if(model %in% c("ln.nb", "ln.ln")) {
            rlnorm(n=n.samples, meanlog=
                       log(dcounts.case.true[i,] /
                               sqrt(1+sd.dna^2/dcounts.case.true[i,]^2)),
                   sdlog=sqrt(log(1+sd.dna^2/dcounts.case.true[i,]^2)))
        } else if(model %in% c("gamma.pois")) {
            rgamma(n=n.samples, shape=(dcounts.case.true[i,])^2/sd.dna.cond^2, 
                   scale=sd.dna.cond^2/dcounts.case.true[i,])
        }
    }))
    rownames(dcounts.case.true) <- paste0("case_", seq_len(n.case))
    rownames(dcounts.case.obs) <- paste0("case_", seq_len(n.case))
    if(n.ctrl > 0){
        par.dna.mu.ctrl <- rlnorm(
            n=n.ctrl, meanlog=log(mu.dna / sqrt(1+mu.sd.dna^2/mu.dna^2)), 
            sdlog=sqrt(log(1+mu.sd.dna^2/mu.dna)))
        par.dna.bc.ctrl <- do.call(rbind, lapply(seq_len(n.ctrl), function(i) {
            exp(rnorm(n = n.bc, mean = 0, sd = sd.dna.bc))[idx.bc]
        }))
        par.dna.cond.ctrl <- do.call(rbind, lapply(seq_len(n.ctrl), function(i) {
            exp(rnorm(n = n.cond, mean = 0, sd = sd.dna.bc))[idx.cond]
        }))
        dcounts.ctrl.true <- do.call(rbind, lapply(seq_len(n.ctrl), function(i) {
            par.dna.mu.ctrl[i]*par.dna.bc.ctrl[i,]*par.dna.cond.ctrl[i,]
        }))
        dcounts.ctrl.obs <- do.call(rbind, lapply(seq_len(n.ctrl), function(i) {
            if(model %in% c("ln.nb", "ln.ln")) {
                rlnorm(n=n.samples, meanlog=
                           log(dcounts.ctrl.true[i,] /
                                   sqrt(1+sd.dna^2/dcounts.ctrl.true[i,]^2)),
                       sdlog=sqrt(log(1+sd.dna^2/dcounts.ctrl.true[i,]^2)))
            } else if(model %in% c("gamma.pois")) {
                rgamma(n=n.samples, shape=(dcounts.ctrl.true[i,])^2/sd.dna.cond^2, 
                       scale=sd.dna.cond^2/dcounts.ctrl.true[i,])
            }
        }))
        rownames(dcounts.ctrl.obs) <- paste0("ctrl_", seq_len(n.ctrl))
        rownames(dcounts.ctrl.true) <- paste0("ctrl_", seq_len(n.ctrl))
    } else {
        par.dna.mu.ctrl <- NULL
        par.dna.cond.ctrl <- NULL
        par.dna.bc.ctrl <- NULL
        dcounts.ctrl.true <- NULL
        dcounts.ctrl.obs <- NULL
    }
    par.dna.mu <- c(par.dna.mu.case, par.dna.mu.ctrl)
    par.dna.cond <- rbind(par.dna.cond.case, par.dna.cond.ctrl)
    par.dna.bc <- rbind(par.dna.bc.case, par.dna.bc.ctrl)
    dcounts.true <- rbind(dcounts.case.true, dcounts.ctrl.true)
    dcounts.obs <- rbind(dcounts.case.obs, dcounts.ctrl.obs)
    colnames(dcounts.true) <- paste0("S", seq_len(n.samples))
    colnames(dcounts.obs) <- paste0("S", seq_len(n.samples))
    #dcounts.obs[dcounts.obs <= 1] <- 1
    
    ## draw rna counts
    if(n.ctrl > 0){
        # only draw one underlying mean and one vector of condition effects
        # as the model assumes that the controls follow the same
        # rna generating process
        par.rna.mu.ctrl <- mu.rna*exp(rnorm(n=1, mean=0, sd=mu.sd.rna))
        par.rna.cond.ctrl <- exp(rnorm(n=n.cond, mean=0, sd=sd.rna.cond))[idx.cond]
        rcounts.ctrl.true <- do.call(rbind, lapply(seq_len(n.ctrl), function(i) {
            dcounts.ctrl.true[i,]*par.rna.mu.ctrl*par.rna.cond.ctrl
        }))
        rcounts.ctrl.obs <- do.call(rbind, lapply(seq_len(n.ctrl), function(i) {
            if(model %in% c("ln.ln")) {
                rlnorm(n=n.samples, meanlog=
                           log(rcounts.ctrl.true[i,] /
                                   sqrt(1+sd.rna^2/rcounts.ctrl.true[i,]^2)),
                       sdlog=sqrt(log(1+sd.rna^2/rcounts.ctrl.true[i,]^2)))
            } else if(model %in% c("ln.nb")) {
                rnbinom(n=n.samples, mu = rcounts.ctrl.true[i,], 
                        size=sd.rna )
            } else if(model %in% c("gamma.pois")) {
                rnbinom(n=n.samples, mu = rcounts.ctrl.true[i,], 
                        size=(dcounts.ctrl.true[i,])^2/sd.dna.cond^2)
            }
        }))
        
        par.rna.de.case <- rbinom(n=n.case, size=1, prob=frac.de)==1
        par.rna.mu.case <- mu.rna*exp(rnorm(n=n.case, mean=0, sd=mu.sd.rna))
        par.rna.mu.case[par.rna.de.case] <- par.rna.mu.ctrl
        par.rna.cond.case <- do.call(rbind, lapply(seq_len(n.case), function(i) {
            if(par.rna.de.case[i] & n.cond > 1){
                par.rna.cond.ctrl * 
                    exp(rnorm(n=n.cond, mean=0, sd=sd.rna.cond_case))[idx.cond]
            } else {
                par.rna.cond.ctrl
            }
        }))
        rcounts.case.true <- do.call(rbind, lapply(seq_len(n.case), function(i) {
            dcounts.case.true[i,]*par.rna.mu.case[i]*par.rna.cond.case[i,]
        }))
        rcounts.case.obs <- do.call(rbind, lapply(seq_len(n.case), function(i) {
            if(model %in% c("ln.ln")) {
                rlnorm(n=n.samples, meanlog=
                           log(rcounts.case.true[i,] /
                                   sqrt(1+sd.rna^2/rcounts.case.true[i,]^2)),
                       sdlog=sqrt(log(1+sd.rna^2/rcounts.case.true[i,]^2)))
            } else if(model %in% c("ln.nb")) {
                rnbinom(n=n.samples, mu = rcounts.case.true[i,], 
                        size=sd.rna )
            } else if(model %in% c("gamma.pois")) {
                rnbinom(n=n.samples, mu = rcounts.case.true[i,], 
                        size=(dcounts.case.true[i,])^2/sd.dna.cond^2)
            }
        }))
    } else {
        par.dna.mu.ctrl <- NULL
        par.dna.cond.ctrl <- NULL
        par.dna.bc.ctrl <- NULL
        dcounts.ctrl.true <- NULL
        dcounts.ctrl.obs <- NULL
        
        par.rna.mu.case <- mu.rna*exp(rnorm(n=n.case, mean=0, sd=mu.sd.rna))
        par.rna.de.case <- rbinom(n=n.case, size=1, prob=frac.de)==1
        par.rna.cond.case <- do.call(rbind, lapply(seq_len(n.case), function(i) {
            if(par.rna.de.case[i] & n.cond > 1){
                exp(rnorm(n=n.cond, mean=0, sd=sd.rna.cond))[idx.cond]
            } else {
                rep(1, length(idx.cond))
            }
        }))
        rcounts.case.true <- do.call(rbind, lapply(seq_len(n.case), function(i) {
            dcounts.case.true[i,]*par.rna.mu.case[i]*par.rna.cond.case[i,]
        }))
        rcounts.case.obs <- do.call(rbind, lapply(seq_len(n.case), function(i) {
            if(model %in% c("ln.ln")) {
                rlnorm(n=n.samples, meanlog=
                           log(rcounts.case.true[i,] /
                                   sqrt(1+sd.rna^2/rcounts.case.true[i,]^2)),
                       sdlog=sqrt(log(1+sd.rna^2/rcounts.case.true[i,]^2)))
            } else if(model %in% c("ln.nb")) {
                rnbinom(n=n.samples, mu = rcounts.case.true[i,], 
                        size=sd.rna )
            } else if(model %in% c("gamma.pois")) {
                rnbinom(n=n.samples, mu = rcounts.case.true[i,], 
                        size=(dcounts.case.true[i,])^2/sd.dna.cond^2)
            }
        }))
    }
    par.rna.mu <- c(par.rna.mu.case, par.rna.mu.ctrl)
    par.rna.cond <- rbind(par.rna.cond.case, par.rna.cond.ctrl)
    rcounts.true <- rbind(rcounts.case.true, rcounts.ctrl.true)
    rcounts.obs <- rbind(rcounts.case.obs, rcounts.ctrl.obs)
    
    rownames(rcounts.true) <- rownames(dcounts.true)
    colnames(rcounts.true) <- colnames(dcounts.true)
    rownames(rcounts.obs) <- rownames(dcounts.obs)
    colnames(rcounts.obs) <- colnames(dcounts.obs)
    
    ## create annotation
    colAnnot <- data.frame(
        sample=colnames(dcounts.obs),
        cond=paste0("cond_", seq(1, n.cond))[idx.cond],
        rep=paste0("rep_", seq(1, n.reps))[idx.rep],
        barcode=paste0("BC", seq(1, n.bc))[idx.bc],
        stringsAsFactors = FALSE
    )
    colAnnot$experiment <- colAnnot$cond # for depth factors
    rowAnnot <- data.frame(
        type=c(rep("case", n.case),
               rep("ctrl", n.ctrl)),
        stringsAsFactors = FALSE
    )
    
    return(list(dcounts.obs=dcounts.obs,
                dcounts.true=dcounts.true,
                par.dna=list(par.dna.mu=par.dna.mu,
                             par.dna.cond=par.dna.cond,
                             par.dna.bc=par.dna.bc),
                rcounts.obs=rcounts.obs,
                rcounts.true=rcounts.true,
                par.rna=list(par.rna.mu=par.rna.mu,
                             par.rna.cond=par.rna.cond,
                             par.rna.de.case=par.rna.de.case),
                colAnnot=colAnnot,
                rowAnnot=rowAnnot))
}
