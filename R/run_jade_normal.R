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
	  f0 <- jade_admm(y=y, gamma=0, lambda=NULL, sample.size=R$sample.size)
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

