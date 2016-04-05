#'Run alternative methods for binomial simulations
#'@description Run alternative methods BSmooth,
#'MethylKit and some t-test based methods on binomial data
#'@param file.prefix File prefixe
#'@param bs.factor Separation between sites used for BSmooth
#'@param tmpdir Directory to use for temporary files.
#'@return Nothing
#'@export
alternatives_binomial <- function(file.prefix, bs.factor=5, tmpdir=".",
                                  sites=1:p, bs.ns=70){
  data.file <- paste0("data/", file.prefix, "_data.RData")
  out.file <- paste0("alt/", file.prefix, "_altpvals.RData")
  #cat(out.file, "\n")
  R <- getobj(data.file)

  p <- dim(R$Y)[1]
  k <- dim(R$Y)[2]
  K <- length(R$sample.size)

  y <- reads <- yhat <- matrix(nrow=p, ncol=K)
  strt <- c(1)
  stp <- c()
  for(i in 1:K){
    stp <- c(stp, strt[i] + R$sample.size[i] -1)
    y[, i] <- rowSums(R$Y[, strt[i]:stp[i]])
    reads[, i] <- rowSums(R$READS[, strt[i]:stp[i]])
    yhat[,i] <- (y[,i] + 0.5)/(reads[,i] + 1)
    if(i < K) strt <- c(strt, strt[i] + R$sample.size[i])
  }

  #MethylKit
  mkit.dataframes <- list()
  m.tmp.files <- c()
  for(i in 1:K){
    m.tmp.files[i] <- tempfile(pattern="mkit", tmpdir=tmpdir)
    mkit.dataframes[[i]] <- data.frame(chrBase=paste("1", 1:p, sep=""),
                                       chr=rep(1, p), base=1:p, strand=rep("R", p),
                                       coverage=reads[,i],
                                       freqC=round(y[,i]/reads[,i]*100, digits=2),
                                       freqT=round(100-(y[,i]/reads[,i]*100), digits=2))
    write.table(mkit.dataframes[[i]], file=m.tmp.files[i], quote=FALSE, row.names=FALSE)
  }
  myobj <- read(list(m.tmp.files[1], m.tmp.files[2]),
                sample.id=list("H", "D"), treatment=c(0, 1), assembly="")
  merge_myobj <- unite(myobj, destrand=TRUE)
  agg.dmr <- calculateDiffMeth(merge_myobj)
  unlink(m.tmp.files)

  ##Run methylKit on individual data
  mfiles <- as.list(replicate(k, tempfile(pattern="mkit", tmpdir=tmpdir)))
  for(i in 1:k){
    mkit <- data.frame(chrBase=paste("1", 1:p, sep=""),
                       chr=rep(1, p),
                       base=1:p,
                       strand=rep("R", p),
                       coverage=R$READS[,i],
                       freqC=round(R$Y[,i]/R$READS[,i]*100, digits=2),
                       freqT=round(100-(R$Y[,i]/R$READS[,i]*100), digits=2))
    write.table(mkit, file=mfiles[[i]], quote=FALSE, row.names=FALSE)
  }
  myobj <- read(mfiles, sample.id=as.list(c(paste0("H", 1:R$sample.size[1]), paste0("D", 1:R$sample.size[2]))),
                treatment=c(rep(0, R$sample.size[1]), rep( 1, R$sample.size[2])), assembly="")
  merge_myobj <- unite(myobj, destrand=TRUE)
  ind.dmr <- calculateDiffMeth(merge_myobj)
  for(i in 1:k) unlink(mfiles[[i]])

  mkit.pvals <- data.frame(agg.dmr$pvalue, agg.dmr$qvalue, ind.dmr$pvalue, ind.dmr$qvalue)

  #MethylKit aggregated and individual, pvals and qvals
  names(mkit.pvals) <- c("mk.agg.pval", "mk.agg.qval", "mk.ind.pval", "mk.ind.qval")
  #4

  #BSSEQ
  B.raw <- BSseq( chr=rep(1, p), pos=round(sites*bs.factor, digits=0),
                  M=R$Y, Cov=R$READS,
                  sampleNames=c(paste0("H", 1:R$sample.size[1]), paste0("D", 1:R$sample.size[2])))
  B.smooth <- BSmooth(B.raw, ns=bs.ns)
  B.tstat <- BSmooth.tstat(B.smooth, group1=paste0("H", 1:R$sample.size[1]) ,
                           group2=paste0("D", 1:R$sample.size[2]))
  bsseq.tstat <- B.tstat@stats[,"tstat"]
  bsseq.pvals <- dt(bsseq.tstat, df=k-1)
  bsseq.slim.pi0 <- SLIMfunc(bsseq.pvals)
  bsseq.slim <- QValuesfun(bsseq.pvals, bsseq.slim.pi0$pi0_Est)
  bsseq.stats <- data.frame(bsseq.tstat, bsseq.pvals, bsseq.slim)
  names(bsseq.stats) <- c("bss.tstat", "bss.pval", "bss.qval")
  #3


  #Spline + T-Tests
  spline.Y <- apply(R$Y/R$READS, MARGIN=2, FUN=spline.func, sites=sites)
  #Weighted spline t-tests using resampling (stat only, no pvalue)
  spline.wt.ttests <- fit_withttest_binom(yhat, reads, sites, B=100, type="spline")
  #Regular t-tests
  spline.pvals <- apply(spline.Y, MARGIN=1, FUN=ttest.func, sample.size=R$sample.size)

  #Locfit + T-test
  #locfit.binom.func uses the same parameters as bsmooth
  locfit.binom.Y <-  apply(rbind(R$Y, R$READS), MARGIN=2,
                           FUN=locfit.binom.func, sites=round(sites*bs.factor, digits=0), ns=bs.ns)
  locfit.wt.ttests <- fit_withttest_binom(yhat, reads, sites, B=100, type="locfit", bs.ns=bs.ns)
  locfit.pvals <- apply(locfit.binom.Y, MARGIN=1, FUN=ttest.func, sample.size=R$sample.size)


  #Adjust all pvalues using SLIM, BH

  spline.slim.pi0 <- SLIMfunc(spline.pvals)
  spline.slim <- QValuesfun(spline.pvals, spline.slim.pi0$pi0_Est)

  locfit.slim.pi0 <- SLIMfunc(locfit.pvals)
  locfit.slim <- QValuesfun(locfit.pvals, locfit.slim.pi0$pi0_Est)


  my.ttests <- data.frame(spline.wt.ttests, spline.pvals, spline.slim,
                          locfit.wt.ttests, locfit.pvals, locfit.slim)
  #6
  all.stats <- cbind(my.ttests, mkit.pvals, bsseq.stats)
  save(all.stats , file=out.file)
}


#Utility functions for alternative methods

#Weights t-test to reflect variability of fit resampling binomial
fit_withttest_binom <- function(phat, reads, sites, B=100, type=c("spline", "locfit"), bs.ns=70){
  type <- match.arg(type)
  p <- dim(reads)[1]
  avg <- matrix(0, p , 2)
  var <- matrix(0, p , 2)
  for(j in 1:2){
    fits <- matrix(0, p, B)
    i <- 1
    while(i <=B){
      myy <- rbinom(n=p, size=reads[,j], prob=phat[,j])
      if(type=="spline"){
        fits[,i] <- spline.func(myy/reads[,j], sites)
      }else if(type=="locfit"){
        fits[,i] <- locfit.binom.func(c(myy,reads[,j]), round(sites*1e6, digits=0), ns=bs.ns)
      }
      i <- i+1
    }
    avg[,j] <- rowMeans(fits)
    var[,j]<- (1/(B-1))*rowSums(apply(fits, MARGIN=2, FUN=function(x, avg){
      return((x-avg)^2)
    }, avg=avg[,j]))
  }
  ttests <- unlist(lapply(1:p, FUN=function(idx, avg, var){
    t <- (avg[idx,1]-avg[idx,2])/sqrt(var[idx,1] + var[idx,2])
    return(t)}, avg=avg, var=var))
  return(ttests)

}

ttest.func <- function(x, sample.size){
  labs <- c(rep(0, sample.size[1]), rep(1,sample.size[2]))
  s1 <- sum(!is.na(x[labs==0]))
  s2 <- sum(!is.na(x[labs==1]))
  if(s1 < 3 | s2 < 3) return(NA)
  t <- t.test(x~ labs)
  return(t$p.value)
}

spline.func <- function(x, sites){
  ix <- which(!is.na(x))
  f <- smooth.spline(x=sites[ix], y=x[ix], cv=TRUE)
  y <- predict(f, x=sites)$y
  return(y)
}

locfit.binom.func <- function(x, sites, ns){
  p <- length(sites)
  y <- x[1:p]; reads <- x[(p+1):(2*p)]

  sdata <- data.frame(pos = sites, M = pmin(pmax(y, 0.01), reads - 0.01), Cov = reads)
  fit <- locfit(M~lp(pos, nn=(ns/p), h=1000),
                weights=reads, family="binomial", data=sdata, maxk=10000)
  pp <- preplot(fit, where="data")
  return(expit(pp$fit))

}
