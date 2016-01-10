#'Run alternative methods for normal simulations
#'@description Run alternative t-test based methods on normal data
#'@param file.prefix File prefixe
#'@return Nothing
#'@export
alternatives_normal <- function(file.prefix){
  data.file <- paste("data/", file.prefix, "_data.RData", sep="")
  out.file <- paste("alt/", file.prefix, "_altpvals.RData", sep="")
  #cat(out.file, "\n")
  R <- jadeTF:::getobj(data.file)

  p <- dim(R$Y)[1]
  K <- length(R$sample.size)
  stopifnot(dim(R$Y)[2] == sum(R$sample.size))
  stopifnot(K==2)
  sites <- 1:p

  strt <- c(1)
  stp <- c()
  for(i in 1:K){
    stp <- c(stp, strt[i] + R$sample.size[i] -1)
    if(i < K) strt <- c(strt, strt[i] + R$sample.size[i])
  }


  #No smoothing + T-test
  ttests <- apply(R$Y, MARGIN=1, FUN=ttest.func, sample.size=R$sample.size)

  #Spline + T-Tests
  spline.Y <- matrix(nrow=p, ncol=sum(R$sample.size))
  for(i in 1:K){
    spline.Y[, strt[i]:stp[i]] <- apply(R$Y[, strt[i]:stp[i]], MARGIN=2,
                                        FUN=spline.func, sites=sites)
  }
  spline.ttests <- apply(spline.Y, MARGIN=1, FUN=ttest.func, sample.size=R$sample.size)

  spline.nv.ttests <- fit_withttest(R$Y, R$sample.size, sites, B=100, type="spline")

  #Locfit + T-test
  locfit.Y <- matrix(nrow=p, ncol=sum(R$sample.size))
  for(i in 1:K){
    locfit.Y[, strt[i]:stp[i]] <- apply(R$Y[, strt[i]:stp[i]], MARGIN=2, FUN=locfit.func, sites=sites)
  }
  locfit.ttests <- apply(locfit.Y, MARGIN=1, FUN=ttest.func, sample.size=R$sample.size)


  locfit.nv.ttests <- fit_withttest(R$Y, R$sample.size, sites, B=100, type="locfit")


  #Adjust all pvalues using SLIM, BH
  tt.slim.pi0 <- SLIMfunc(ttests)
  tt.slim <- QValuesfun(ttests, tt.slim.pi0$pi0_Est)
  tt.bh <- p.adjust(ttests, method="BH")

  spline.tt.slim.pi0 <- SLIMfunc(spline.ttests)
  spline.tt.slim <- QValuesfun(spline.ttests, spline.tt.slim.pi0$pi0_Est)
  spline.tt.bh <- p.adjust(spline.ttests, method="BH")

  locfit.tt.slim.pi0 <- SLIMfunc(locfit.ttests)
  locfit.tt.slim <- QValuesfun(locfit.ttests, locfit.tt.slim.pi0$pi0_Est)
  locfit.tt.bh <- p.adjust(locfit.ttests, method="BH")

  #Pvals
  stats <- data.frame(ttests, tt.slim, tt.bh, spline.ttests, spline.tt.slim, spline.tt.bh, spline.nv.ttests, locfit.ttests, locfit.tt.slim, locfit.tt.bh, locfit.nv.ttests)
  save(stats, file=out.file)
}


locfit.gcv.func <- function(alpha, x, sites){
  z <- gcv(x~sites, alpha=alpha)
  return(z["gcv"])
}

locfit.func <- function(x, sites){
  alpha.opt <- optimize(f=locfit.gcv.func, interval=c(0, 0.1), x=x, sites=sites, maximum=FALSE)
  alpha <- alpha.opt$minimum
  fit <- locfit(x~sites, alpha=alpha)
  pp <- preplot(fit, where="data")
  return(pp$fit)
}

fit_withttest <- function(Y, sample.size, sites, B=100, type=c("spline", "locfit")){
  type <- match.arg(type)
  p <- dim(Y)[1]
  avg <- matrix(0, p , 2)
  var <- matrix(0, p , 2)
  strt <- 1
  for(j in 1:2){
    stp <- strt + sample.size[j] -1
    y <- rowMeans(Y[, strt:stp])
    sd_hat <- apply(Y[, strt:stp], MARGIN=1, FUN=sd)
    fits <- matrix(0, p, B)
    i <- 1
    while(i <=B){
      myy <- rnorm(n=p, mean=y, sd=sd_hat)
      if(type=="spline"){
        fits[,i] <- spline.func(myy, sites)
      }else if(type=="locfit"){
        fits[,i] <- locfit.func(myy, sites)
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
