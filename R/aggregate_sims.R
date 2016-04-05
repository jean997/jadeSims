
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
#'   \item{\code{avg.tpr, avg.fpr}}{Average true and false positive rates at given level and for JADE (14th column)}
#' }
#'@export
aggregate_sims <- function(file.prefix, profiles, save.file=NULL,
                                 which.reps=1:60, tol=5e-3, level=0.1,
                                 run.cv=TRUE, clean.jade=FALSE, use.cv=TRUE){
  n.sims <- length(which.reps)
  p <- dim(profiles)[1]
  K <- dim(profiles)[2]
  site.labels <- abs(profiles[,1]-profiles[,2]) > tol

  all.sep <- list()

  j <- 1
  for(rep in which.reps){
    jade.cv.file <- paste0("cv/", file.prefix , "_n", rep , "_cv.RData")
    orig.path.file <- paste0("path/", file.prefix, "_n", rep, "_path.0.RData")
    path.file.list <- paste0("path/", file.prefix, "_n", rep, "_path.", 1:5, ".RData")
    data.file <- paste0("data/", file.prefix, "_n", rep, "_data.RData")
    alt.file <- paste0("alt/", file.prefix, "_n", rep, "_altpvals.RData")

    #Pvals for alternatives
    stats <- getobj(alt.file)
    if(j==1){
      pnames <- names(stats)
      all.stats <- array(dim=c(n.sims, p, length(pnames)))
      all.tpr <- all.fpr <- matrix(nrow=n.sims, ncol=length(pnames)+1)
    }
    all.stats[j, , ] <- as.matrix(stats)

    #FPR and TPR
    fdr.cols <- which(pnames %in% c("mk.agg.qval", "mk.ind.qval",
                                    "bss.qval" ,
                                    "spline.slim", "locfit.slim",
                                    "tt.slim"))
    p.cols <- which(pnames %in% c("mk.agg.pval", "mk.ind.pval", "bss.pval",
                                  "spline.pvals", "locfit.pvals", "tt.pvals"))
    stat.cols <- which(pnames %in% c("bss.tstat", "spline.wt.ttests", " locfit.wt.ttests"))

    for(i in 1:length(pnames)){
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
      all.tpr[j, length(pnames)+1] <- tpr.func(sep[, cv.obj$cv.1se.l1], site.labels)
      all.fpr[j, length(pnames)+1] <- fpr.func(sep[, cv.obj$cv.1se.l1], site.labels)
    }
    j <- j+1
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
