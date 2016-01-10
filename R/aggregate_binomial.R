
#'@title Aggregate binomial simulations
#'@description Aggregate results of binomial simulations
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
#'  \item{\code{all.stats}}{B x n x 13 array of statistic values for alternative methods.}
#'   \item{\code{site.labels}}{Which sites are separated}
#'   \item{\code{tol, level}}{Parameters passed in}
#'   \item{\code{avg.tpr, avg.fpr}}{Average true and false positive rates at given level and for JADE (14th column)}
#' }
#'@export
aggregate_binomial <- function(file.prefix, profiles, save.file=NULL,
                                 which.rep=1:60, tol=5e-3, level=0.1,
                                 run.cv=TRUE){
  n.sims <- length(which.rep)
  p <- dim(profiles)[1]
  K <- dim(profiles)[2]
  site.labels <- abs(profiles[,1]-profiles[,2]) > tol
  all.stats <- array(dim=c(n.sims, p, 13))
  all.sep <- list()
  all.tpr <- matrix(nrow=n.sims, ncol=14)
  all.fpr <- matrix(nrow=n.sims, ncol=14)

  j <- 1
  for(rep in which.reps){
    cat(file.prefix_n, "\n")
    jade.cv.file <- paste0("cv/", file.prefix , "_n", rep , "_cv.RData")
    orig.path.file <- paste0("path/", file.prefix, "_n", rep, "_path.0.RData")
    path.file.list <- paste0("path/", file.prefix, "_n", rep, "_path.", 1:5, ".RData")
    data.file <- paste0("data/", file.prefix, "_n", rep, "_data.RData")
    alt.file <- paste0("alt/", file.prefix, "_n", rep, "_altpvals.RData")

    path <- getobj(orig.path.file)
    #stopifnot(any(abs(path$gammas - gammas) > 1e-5))

    #Separation

    sep <- matrix(unlist(lapply(path$JADE_fits[-1], FUN=function(f, tol){
      z <- get_sep(f$beta, tol=tol)
      return(z[[1]][[1]])
    }, tol=tol)), nrow=p)
    sep <- cbind( get_sep(path$JADE_fits[[1]]$fits, tol), sep)
    all.sep[[j]] <- sep

    #Cross validation for JADE
    if(run.cv){
      cv.obj <- cv_err_wts(orig.path.file, path.file.list,
                           use.converged.only=TRUE, control.l1=TRUE)
      save(cv.obj, file=jade.cv.file)
      all.tpr[j, 14] <- tpr.func(sep[, cv.obj$cv.1se.l1], site.labels)
      all.fpr[j, 14] <- fpr.func(sep[, cv.obj$cv.1se.l1], site.labels)
    }

    #Pvals for alternatives
    stats <- getobj(alt.file)
    all.stats[j, , ] <- as.matrix(stats)

    if(j==1){
      pnames <- names(stats)
    }

    #FPR and TPR
    fdr.cols <- which(pnames %in% c("mk.agg.qval", "mk.ind.qva","bss.qval" , "spline.ind.tt.slim", "locfit.ind.tt.slim"))
    p.cols <- which(pnames %in% c("mk.agg.pval", "mk.ind.pval", "bss.pval", "spline.ind.ttests", "locfit.ind.ttests"))
    stat.cols <- which(pnames %in% c("bss.tstat", "spline.ttests", " locfit.ttests"))
    for(i in 1:13){
      my.labs <- rep(0, p)
      if(i %in% fdr.cols | i %in% p.cols){
        if(any(stats[,i] <=0)){
          cat("Failed: ", pnames[i], "\n")
        }else{
          my.labs[stats[, i] < level] <- 1
        }
      }
      all.tpr[j, i] <- tpr.func(my.labs, site.labels)
      all.fpr[j, i] <- fpr.func(my.labs, site.labels)
    }
    j=j+1
  }
  avg.tpr <- colMeans(all.tpr, na.rm=TRUE)
  avg.fpr <- colMeans(all.fpr, na.rm=TRUE)
  names <- c(pnames, "jade")
  R <- list("avg.tpr"=avg.tpr, "avg.fpr"=avg.fpr,
            "names"=names, "all.sep"=all.sep, "tol"=tol, "level"=level,
            "all.stats"=all.stats, "site.labels"=site.labels)
  if(!is.null(save.file)) save(R, file=save.file)
  return(R)
}
