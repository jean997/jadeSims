#'@title Plot ROC curves for simulations
#'@description Plot ROC curves for JADE and any other methods
#'@param rate.list An object produced by avg_by_ct
#'@param cols Colors
#'@param ltys Line types
#'@param main Plot title
#'@param lwd Line Width
#'@return Nothing
#'@export
plot_rates2 <- function(rate.list, cols, ltys, main="",
                        lwd=1.5, make.legend=FALSE, nbars=5,
                        cex.axis=1, cex.lab=1.5, cex.main=2.5){
  N <- length(rate.list)-1

  whichCI <- seq(1, length(rate.list[[1]]$fct), length.out=nbars)

  m <- max(rate.list[[1]]$fct)
  for(i in 2:N) m <- max(m, rate.list[[i]]$fct)
  plotCI(x=rate.list[[1]]$fct[whichCI], y=rate.list[[1]]$tpr[whichCI],
         uiw=rate.list[[1]]$s.e[whichCI], err="y", xlim=c(0, m), ylim=c(0, 1),
         pch=0, main=main, xlab="False Positives", cex.axis=cex.axis,
         ylab="TPR", cex.lab=cex.lab, col=cols[1], cex.main=cex.main)

  lines(rate.list[[1]]$fct, rate.list[[1]]$tpr, lwd=lwd, col=cols[1], lty=ltys[1])
  for(i in 2:N){
    plotCI(x=rate.list[[i]]$fct[whichCI], y=rate.list[[i]]$tpr[whichCI],
           uiw=rate.list[[i]]$s.e[whichCI], err="y", col=cols[i], add=TRUE, pch=0)
    lines(rate.list[[i]]$fct, rate.list[[i]]$tpr, col=cols[i], lwd=lwd, lty=ltys[i])
  }
  if(make.legend) legend("bottomright", legend=rate.list$names, lty=1, col=cols)
}
