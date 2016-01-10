#'@title Run JADE for one normal simulation
#'@description Run one normal simulation with existing data
#'@param file.prefix File prefix
#'@param run.f0 Run at gamma=0
#'@param folds Which folds to calculate path for (0 is full data)
#'@return Nothing
#'@export
run_jade_normal <- function(file.prefix, run.f0=TRUE, folds=1:5, log.gamma.min=-5){
	fit0.file <- paste0("fit0/", file.prefix, "_f0.0.RData")
	data.file <- paste0("data/", file.prefix, "_data.RData")
	R <- getobj(data.file)
  K <- length(R$sample.size)
	p <- dim(R$Y)[1]
	y <- matrix(nrow=p, ncol=K)
	strt <- 1
	for(i in 1:K){
	  stp <- strt + R$sample.size[i]-1
		my.y <- R$Y[, strt:stp]
		y[,i] <- rowMeans(my.y)
		strt <- strt + R$sample.size[i]
	}
  if(run.f0){
	  f0 <- jade_admm(y=y, gamma=0, lambda=lambda, sample.size=sample.size)
	  save(f0, file=fit0.file)
	  #Fit0 for folds
	  cv_fit0(orig.fit=fit0.file, which.fold=1:5, n.folds=5,
	          save.prefix = paste0("fit0/", file.prefix, "_f0"), return.objects=FALSE)
  }
	#Fit0 for folds and fit_var and path submit
	for(j in folds){
	  fit0.file <- paste0("f0/", file.prefix, "_f0.", j, ".RData")
	  path <- jade_path(fit0=fit0.file, n.fits=100, out.file=out.file,
	                    max.it=10000, tol=5e-3, adjust.rho.alpha=TRUE, log.gamma.min=log.gamma.min)
	}
}


plot_roc_curves <- function(agg.obj, which.stats, cols, direction="vertical", main="", make.legend=FALSE){
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
	plotrix::plotCI(x=avg_jade$fpr[whichCI], y=avg_jade$tpr[whichCI], uiw=avg_jade$s.e[whichCI], err=err, ylim=c(0, 1), xlim=c(0, 1), pch=0, main=main, xlab="FPR", ylab="TPR", cex.lab=1.5, col=c, cex.main=2.5)
  lines(x=avg_jade$fpr, y=avg_jade$tpr, lwd=1.5, col=c)
	points(agg.obj$avg.fpr[which(agg.obj$names=="jade")], agg.obj$avg.tpr[which(agg.obj$names=="jade")], pch=16, col=c, cex=2.5)
	ct <- ct+1


	stat.cols <- which(agg.obj$names %in% c("spline.nv.ttests", "locfit.nv.ttests"))
	for(j in 1:9){
		if(agg.obj$names[j] %in% which.stats){
			c <- cols[which(which.stats==agg.obj$names[j])]
			if( j %in% stat.cols){
				my.pred <- prediction(prediction=abs(t(agg.obj$all.stats[, , j])), labels)
			}else{
				my.pred <- prediction(prediction=-log10(t(agg.obj$all.stats[, , j])), labels)
			}
			my.perf <- performance(my.pred, 'tpr', 'fpr')
			my.interp <- avg_by_interp(tpr.list=my.perf@y.values, fpr.list=my.perf@x.values, direction=direction, npoints=200)
			plotrix::plotCI(x=my.interp$fpr[whichCI], y=my.interp$tpr[whichCI], uiw=my.interp$s.e[whichCI], err=err, col=c, add=TRUE, pch=0)
  		lines(x=my.interp$fpr, y=my.interp$tpr, col=c, lwd=1.5)
			points(agg.obj$avg.fpr[j], agg.obj$avg.tpr[j], pch=16, col=c, cex=2.5)
			ct <- ct+ 1
		}
	}
	if(make.legend) legend("bottomright", legend=which.stats, col=cols, lty=rep(1, length(which.stats)), cex=2)


}



