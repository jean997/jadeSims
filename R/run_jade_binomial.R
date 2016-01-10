#'@import jadeTF
#'@import locfit
#'@import methylKit
#'@import bsseq
#'@import Matrix
#'@import ROCR

#'@title Run JADE for one binomial simulation
#'@description Run one binomial simulation with existing data
#'@param file.prefix File prefix
#'@param run.f0 Run at gamma=0
#'@param folds Which folds to calculate path for (0 is full data)
#'@return A data frame with columns y and reads
#'@export
run_jade_binomial <- function( file.prefix, run.f0=TRUE, folds=0:5, log.gamma.min=-5){
	fit0.file <- paste("f0/", file.prefix, "_f0.0.RData", sep="")
	data.file <- paste("data/", file.prefix, "_data.RData", sep="")

	R <- getobj(data.file)
	p <- dim(R$Y)[1]
	K <- length(R$sample.size)


	sds <- yhat <- y <- matrix(NA, p, K)
	strt <- 1
	for(i in 1:K){
		stp <- strt + R$sample.size[i] -1
		my.counts <- R$Y[, strt:stp]
		my.reads <- R$READS[, strt:stp]
		y[,i] <- rowSums(my.counts)/rowSums(my.reads)
		yhat[,i] <- (rowSums(my.counts)+0.5)/(rowSums(my.reads)+1)
		sds[,i] <- sqrt(yhat[,i]*(1-yhat[,i])/rowSums(my.reads))
		strt <- strt + R$sample.size[i]
	}
	if(run.f0){
	  f0 <- jade_admm(y=y, sds=sds, gamma=0, lambda=NULL, sample.size=c(1,1), ord=2)
	  save(f0, file=fit0.file)
  	cv_fit0(orig.fit=fit0.file, which.fold=1:5, n.fold=5,
  	        save.prefix=paste0("f0/", file.prefix, "_f0"), return.objects=FALSE)
	}
	for(j in folds){
			#path
			fit0.file <- paste0("f0/", file.prefix, "_f0.", j, ".RData")
			out.file <- paste0("path/", file.prefix, "_path.", j, ".RData")
			#Setting tol=5e-3 equal to that used in building the ROC curves ensures the most complete coverage
			path <- jade_path(fit0=fit0.file, n.fits=100, out.file=out.file,
			                  max.it=10000, tol=5e-3, adjust.rho.alpha=TRUE, log.gamma.min=log.gamma.min)
	}
}


plot_roc_curves <- function(agg.obj, site.labels,  which.stats, cols, direction="vertical", main="", legend=FALSE){
	if(direction=="vertical") err="y"
		else err="x"

	n.sims <- dim(agg.obj$all.stats)[1]
	cat(n.sims, "\n")
	p <- dim(agg.obj$all.stats)[2]
	labels <- matrix(rep(agg.obj$site.labels, n.sims), byrow=FALSE, nrow=p)
	whichCI <- seq(1, 200, length.out=11)
	ct <- 1
	c <- cols[which(which.stats=="jade")]
	tpr.list <- list(); fpr.list <- list()
  for(i in 1:n.sims){
    tf <- get_tpr_fpr(sep=agg.obj$all.sep[[i]], labels=agg.obj$site.labels)
    fpr.list[[i]] <- tf$fpr
    tpr.list[[i]] <- tf$tpr
  }
  avg_jade <- avg_by_interp(tpr.list=tpr.list, fpr.list=fpr.list, direction=direction, npoints=200)
	plotCI(x=avg_jade$fpr[whichCI], y=avg_jade$tpr[whichCI], uiw=avg_jade$s.e[whichCI], err=err, ylim=c(0, 1), xlim=c(0, 1), pch=0, main=main, xlab="FPR", ylab="TPR", cex.lab=1.5, col=c, cex.main=2.5)
  lines(x=avg_jade$fpr, y=avg_jade$tpr, lwd=1.5, col=c)
	j <- which(agg.obj$names=="jade")
	points(agg.obj$avg.fpr[j], agg.obj$avg.tpr[j], pch=16, col=c, cex=2.5)
	ct <- ct+1


	stat.cols <- which(agg.obj$names %in% c("bss.tstat", "spline.ttests", " locfit.ttests"))
	for(j in 1:13){
		if(!(agg.obj$names[j] %in% which.stats)) next
		cat(agg.obj$names[j], "..")
		c <- cols[which(which.stats==agg.obj$names[j])]
		if( j %in% stat.cols){
			my.pred <- prediction(prediction=abs(t(agg.obj$all.stats[, , j])), labels)
		}else{
			miss.idx <- rowSums(agg.obj$all.stats[, , j] <=0)
			if(any(miss.idx > 0)) cat("Missing: ", which(miss.idx > 0), "\n")
			my.pred <- prediction(prediction=-log10(t(agg.obj$all.stats[which(miss.idx ==0), , j])), labels[, which(miss.idx ==0)])
		}
		my.perf <- performance(my.pred, 'tpr', 'fpr')
		my.interp <- avg_by_interp(tpr.list=my.perf@y.values, fpr.list=my.perf@x.values, direction=direction, npoints=200)
		plotCI(x=my.interp$fpr[whichCI], y=my.interp$tpr[whichCI], uiw=my.interp$s.e[whichCI], err=err, col=c, add=TRUE, pch=0)
  	lines(x=my.interp$fpr, y=my.interp$tpr, col=c, lwd=1.5)
		points(agg.obj$avg.fpr[j], agg.obj$avg.tpr[j], pch=16, col=c, cex=2.5)
		ct <- ct+ 1
	}
	if(legend) legend("bottomright", legend=which.stats, col=cols, lty=rep(1, length(which.stats)))


}


