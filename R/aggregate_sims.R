
#'@title Aggregate simulations
#'@description Aggregate results of simulations
#'@param file.prefix File prefix
#'@param profiles Mean values
#'@param save.file Optional output file
#'@param which.rep Which simulations (length B)
#'@param tol Threshold for declaring sites separated
#'@param level Level at which to calcluate average true and false positive rate
#'@param run.cv For JADE, calculate cv error if it hasn't been run already.
#'@return A list with elements
#'#' \describe{
#'  \item{\code{names}}{Name of each method in order}
#'  \item{\code{all.stats}}{B x n x F array of statistic values for alternative methods.}
#'   \item{\code{site.labels}}{Which sites are separated}
#'   \item{\code{tol, level}}{Parameters passed in}
#'   \item{\code{avg.tpr, avg.fpr}}{Average true and false positive rates at given level and for JADE (last column)}
#' }
#'@export
aggregate_sims <- function(file.prefix, profiles, which.alt, stat.thresh,
                           save.file=NULL,
                            which.reps=1:100, tol=5e-3,
                            run.cv=TRUE, clean.jade=FALSE, use.cv=TRUE,
                            merge.margin=0, min.length=1, min.acc=0, max.prop=0.5){
  stopifnot(length(which.alt) == length(stat.thresh))
  n.alt = length(which.alt)
  n.sims <- length(which.reps)
  p <- dim(profiles)[1]
  K <- dim(profiles)[2]
  site.labels <- abs(profiles[,1]-profiles[,2]) > 1e-6

  all.sep <- list()
  all.stats <- array(dim=c(n.sims, p, n.alt))
  all.tpr <- all.fpr <- matrix(nrow=n.sims, ncol=n.alt+1)
  all.region.tpr <- all.region.fct <- matrix(nrow=n.sims, ncol=n.alt+1)


  fdr.cols <- which(which.alt %in% c("mk.agg.qval", "mk.ind.qval",
                                     "bss.qval" ,"spline.slim", "locfit.slim","tt.slim"))
  p.cols <- which(which.alt %in% c("mk.agg.pval", "mk.ind.pval", "bss.pval",
                                   "spline.pvals", "locfit.pvals", "tt.pvals"))

  for(j in 1:length(which.reps)){
    rep <- which.reps[j]
    jade.cv.file <- paste0("cv/", file.prefix , "_n", rep , "_cv.RData")
    orig.path.file <- paste0("path/", file.prefix, "_n", rep, "_path.0.RData")
    path.file.list <- paste0("path/", file.prefix, "_n", rep, "_path.", 1:5, ".RData")
    data.file <- paste0("data/", file.prefix, "_n", rep, "_data.RData")
    alt.file <- paste0("alt/", file.prefix, "_n", rep, "_altpvals.RData")

    #Pvals for alternatives
    stats <- getobj(alt.file)
    ix <- match(which.alt, names(stats))
    all.stats[j, , ] <- as.matrix(stats)[, ix]

    #FPR and TPR

    #stat.cols <- which(pnames %in% c("bss.tstat", "spline.wt.ttests", " locfit.wt.ttests"))

    for(i in 1:length(which.alt)){
      if(is.na(stat.thresh[i])) next
      my.labs <- rep(0, p)
      if((i %in% fdr.cols | i %in% p.cols) & any(stats[,ix[i]] <=0)){
          cat("Failed: ", which.alt[i], "\n")
      }else{
        my.labs[abs(stats[, ix[i]]) < stat.thresh[i]] <- 1
      }
      all.tpr[j, i] <- tpr.func(my.labs, site.labels)
      all.fpr[j, i] <- fpr.func(my.labs, site.labels)
      rr <- rates_by_region(my.labs, site.labels)
      all.region.tpr[j, i] <- rr["tpr"]
      all.region.fct[j, i] <- rr["fp.ct"]
    }
    i <- i  + 1
    #JADE
    path <- getobj(orig.path.file)
    u <- order(path$gammas[-1])
    sep <- matrix(unlist(lapply(path$JADE_fits[-1], FUN=function(f, tol){
      z <- get_sep(f$beta, tol=tol)
      return(z[[1]][[1]])
    }, tol=tol)), nrow=p)
    sep <- sep[,u]
    sep <- cbind( unlist(get_sep(path$JADE_fits[[1]]$fits, tol)), sep)
    if(clean.jade){
      sep.clean <- t( apply(sep, MARGIN=1, FUN=cummin))
      sep <- sep.clean
    }
    all.sep[[j]] <- sep

    #CV
    if(use.cv){
      if(run.cv){
        cv.obj <- cv_err_wts(orig.path.file, path.file.list,
                         use.converged.only=TRUE, control.l1=TRUE)
        save(cv.obj, file=jade.cv.file)
      }else{
        cv.obj <- getobj(jade.cv.file)
      }
      my.labs <- sep[, cv.obj$cv.1se.l1]
      all.tpr[j, i] <- tpr.func(my.labs, site.labels)
      all.fpr[j, i] <- fpr.func(my.labs, site.labels)
      rr <- rates_by_region(my.labs, site.labels)
      all.region.tpr[j, i] <- rr["tpr"]
      all.region.fct[j, i] <- rr["fp.ct"]
    }
  }


  avg.tpr <- colMeans(all.tpr, na.rm=TRUE)
  avg.fpr <- colMeans(all.fpr, na.rm=TRUE)
  avg.region.tpr <- colMeans(all.region.tpr, na.rm=TRUE)
  avg.region.fct <- colMeans(all.region.fct, na.rm=TRUE)
  names <- c(which.alt, "jade")
  R <- list("avg.tpr"=avg.tpr, "avg.fpr"=avg.fpr,
            "avg.region.tpr"=avg.region.tpr, "avg.region.fct"=avg.region.fct,
            "names"=names, "all.sep"=all.sep, "tol"=tol, "stat.thresh"=stat.thresh,
            "all.stats"=all.stats, "site.labels"=site.labels)
  if(!is.null(save.file)) save(R, file=save.file)
  return(R)
}
