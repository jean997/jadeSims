#'@title Plot ROC curves for simulations
#'@description Plot ROC curves for JADE and any other methods
#'@param agg.obj An object produced by aggregate_sims
#'@param which.stats Names of statistics to plot
#'@param cols Colors
#'@param direction Which direction should the error bars go
#'@param main Plot title
#'@param make.legend Make a legend?
#'@return Nothing
#'@export
plot_roc_curves <- function(agg.obj,
                            which.stats, cols, ltys, lwd=1.5,
                            points.pch=16,
                            direction="vertical", n.bars=11, point.cex=2.5,
                            main="", make.legend=FALSE){
  if(direction=="vertical") err="y"
  else err="x"
  if(length(points.pch)==1) points.pch=rep(points.pch, length(which.stats))

  n.sims <- dim(agg.obj$all.stats)[1]
  cat(n.sims, "\n")
  p <- dim(agg.obj$all.stats)[2]
  labels <- matrix(rep(agg.obj$site.labels, n.sims), byrow=FALSE, nrow=p)
  whichCI <- seq(1, 200, length.out=n.bars)
  ct <- 1
  #plot JADE
  stopifnot("jade" %in% which.stats)
  c <- cols[which(which.stats=="jade")]
  l <- ltys[which(which.stats=="jade")]
  p <- points.pch[which(which.stats=="jade")]
  tpr.list <- list()
  fpr.list <- list()
  for(i in 1:n.sims){
    tf <- get_tpr_fpr(sep=agg.obj$all.sep[[i]], labels=agg.obj$site.labels)
    fpr.list[[i]] <- tf$fpr
    tpr.list[[i]] <- tf$tpr
  }
  avg_jade <- jadeSims:::avg_by_interp(tpr.list=tpr.list, fpr.list=fpr.list,
                                       direction=direction, npoints=200)
  plotCI(x=avg_jade$fpr[whichCI], y=avg_jade$tpr[whichCI],
         uiw=avg_jade$s.e[whichCI], err=err, ylim=c(0, 1),
         xlim=c(0, 1), pch=0, main=main, xlab="FPR", ylab="TPR",
         cex.lab=1.5, col=c, cex.main=2.5)
  lines(x=avg_jade$fpr, y=avg_jade$tpr, lwd=lwd, col=c, lty=l)
  j <- which(agg.obj$names=="jade")
  points(agg.obj$avg.fpr[j], agg.obj$avg.tpr[j], pch=p, col=c, cex=point.cex)
  ct <- ct+1


  stat.cols <- which(agg.obj$names %in% c("bss.tstat", "spline.wt.ttests", " locfit.wt.ttests",
                                          "cfdr-pois", "cfdr-huber", "cfdr-tt"))
  #Plot the other methods
  for(j in 1:(length(agg.obj$names)-1)){
    if(!(agg.obj$names[j] %in% which.stats)) next
    cat(agg.obj$names[j], "..")
    c <- cols[which(which.stats==agg.obj$names[j])]
    l <- ltys[which(which.stats==agg.obj$names[j])]
    p <- points.pch[which(which.stats==agg.obj$names[j])]
    if( j %in% stat.cols){
      my.pred <- prediction(prediction=abs(t(agg.obj$all.stats[, , j])), labels)
    }else{
      miss.idx <- rowSums(agg.obj$all.stats[, , j] <=0)
      if(any(miss.idx > 0)) cat("Missing: ", which(miss.idx > 0), "\n")
      my.pred <- prediction(prediction=-log10(t(agg.obj$all.stats[which(miss.idx ==0), , j])),
                            labels[, which(miss.idx ==0)])
    }
    my.perf <- performance(my.pred, 'tpr', 'fpr')
    my.interp <- avg_by_interp(tpr.list=my.perf@y.values, fpr.list=my.perf@x.values,
                               direction=direction, npoints=200)
    plotCI(x=my.interp$fpr[whichCI], y=my.interp$tpr[whichCI],
           uiw=my.interp$s.e[whichCI], err=err, col=c, add=TRUE, pch=0)
    lines(x=my.interp$fpr, y=my.interp$tpr, col=c, lwd=lwd, lty=l)
    points(agg.obj$avg.fpr[j], agg.obj$avg.tpr[j], pch=p, col=c, cex=point.cex)
    ct <- ct+ 1
  }
  if(make.legend){
    legend("bottomright", legend=which.stats, col=cols,
                         lty=rep(1, length(which.stats)))
  }
}
