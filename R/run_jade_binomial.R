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
#'@return Nothing
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


