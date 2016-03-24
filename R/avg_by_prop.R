#For ROC curves
#'@title Get average tpr and fpr rates over replicates. Align rates for fits with the same proportion detected.
#'@param tpr.list A list of vectors giving true positive rates
#'@param fpr.list A list of vectors giving false positive rates
#'@param prop.list A list of vectors giving proportion discovered
#'@param npoints How many points to evaluate at
#'@return list with fpr, tpr and s.e
#'@export
avg_by_prop <- function(tpr.list, fpr.list, prop.list, npoints=200){
  B <- length(tpr.list)
  prop.out <- seq(0, 1, length.out=npoints)
  fpr.mat <-tpr.mat <-  matrix(nrow=npoints, ncol=B)
  for(i in 1:B){
    apprx.tpr <- approx(x=prop.list[[i]], y=tpr.list[[i]],
                        xout=prop.out, rule=2, ties=max)
    tpr.mat[,i] <- apprx.tpr$y
    apprx.fpr <- approx(x=prop.list[[i]],
                        y=fpr.list[[i]], xout=prop.out, rule=2, ties=max)
    fpr.mat[,i] <- apprx.fpr$y
  }
  avg.tpr <- rowMeans(tpr.mat, na.rm=TRUE)
  avg.fpr <- rowMeans(fpr.mat, na.rm=TRUE)
  tot.obs.tpr <- rowSums(!is.na(tpr.mat))
  var.tpr <- (1/(tot.obs.tpr-1))*rowSums((tpr.mat-avg.tpr)^2, na.rm=TRUE)
  tot.obs.fpr <- rowSums(!is.na(fpr.mat))
  var.fpr <- (1/(tot.obs.fpr-1))*rowSums((fpr.mat-avg.fpr)^2, na.rm=TRUE)
  return(list("fpr"=avg.fpr, "tpr"=avg.tpr, "prop"=prop.out,
              "s.e.tpr"=sqrt(var.tpr), "s.e.fpr"=sqrt(var.fpr)))
}
