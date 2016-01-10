perm_results_normal <- function(file.prefix.list, profiles, tol,
                                levels=c(0.05, 0.1),
                                save.file=NULL, run.cv=FALSE){
  site.labels <- abs(profiles[,1]-profiles[,2]) > tol
  all.tpr <- matrix(nrow=length(file.prefix.list), ncol=length(levels))
  all.fpr <- matrix(nrow=length(file.prefix.list), ncol=length(levels))
  for(i in 1:length(file.prefix.list)){
    f <- file.prefix.list[i]
    res <- getobj(paste0("perm_res/", f, "_perm.res.RData"))
    path <- getobj(paste0("path/", f, "_path.0.RData"))
    p <- dim(path$JADE_fits[[1]]$y)[1]
    path_sep <- get_sep_total(path, tol=tol)

    sep <- matrix(unlist(lapply(path$JADE_fits[-1], FUN=function(f, tol){
      z <- get_sep(f$beta, tol=tol)
      return(z[[1]][[1]])
    }, tol=tol)), nrow=p)
    sep <- cbind( get_sep(path$JADE_fits[[1]]$fits, tol)[[1]][[1]], sep)


    res$gammas <- c(0, res$gammas)
    idx <- match(round(res$gammas, digits=8), round(path$gammas, digits=8))
    Es <- rowMeans(res$sep.total)
    fdr <- Es/path_sep[idx]
    fdr[Es==0 & path_sep[idx]==0] <- 0

    for(l in 1:length(levels)){
      alph <- levels[l]
      if(any(fdr < alph)){
        gamma_fdr <- min(res$gammas[fdr < alph & res$gammas > 0])
        g_idx <- match(round(gamma_fdr, digits=8), round(path$gammas, digits=8))
        all.tpr[i, l] <- tpr.func(sep[, g_idx], site.labels)
        all.fpr[i, l] <- fpr.func(sep[, g_idx], site.labels)
      }else{
        all.tpr[i, l] <- 0; all.fpr[i, l] <- 0
      }
    }
  }
  all.tpr <- data.frame(all.tpr)
  names(all.tpr) <- levels
  all.fpr <- data.frame(all.fpr)
  names(all.fpr) <- levels
  return(list("tpr"=all.tpr, "fpr"=all.fpr))
}
