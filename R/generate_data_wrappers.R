#'@title Wrapper to generate binomial data
#'@description Generate binomial data. Run form a directory with a "data/" subdirectory.
#'@export
generate_data_binomial <- function(seed=1111111, start.rep=1, stop.rep=60){
  set.seed(seed)
  sample.size=c(10, 10)
  data("binom_profiles", package="jadeSims")
  p <- dim(binom_profiles)[1]
  K <- dim(binom_profiles)[2]
  for(re in c(0, 0.02, 0.05, 0.07)){
    for(rep in start.rep:stop.rep){
      file.prefix <- paste0("binom_", re, "_n", rep)
      data.file <- paste0("data/", file.prefix, "_data.RData")
      full.counts <- matrix(nrow=p, ncol=0)
      full.reads <- matrix(nrow=p, ncol=0)
      for(i in 1:K){
        my.data <- replicate(n=sample.size[i], expr={
            z <-  tile_binomial_data(profile=binom_profiles[,i], mean.reads=10,
                                     reff.sd=re, mean.tile.length=30)
            cbind(z$y, z$reads)})
        my.counts <-  my.data[, 1, ]
        my.reads <-  my.data[, 2, ]
        full.counts <- cbind(full.counts, my.counts)
        full.reads <- cbind(full.reads, my.reads)
      }
      R <- list("Y"=full.counts, "READS"=full.reads, "sample.size"=sample.size)
      save(R, file=data.file)
      alternatives_binomial(file.prefix=file.prefix)
    }
  }
}

#'@title Wrapper to generate normal auto-regressive data
#'@description Generate normal AR  data. Run form a directory with a "data/" subdirectory.
#'@export
generate_data_ar <- function(seed=2222222, start.rep=1, stop.rep=60){
  set.seed(seed)
  sample.size=c(10, 10)
  data("normal_profiles", package="jadeSims")
  p <- dim(normal_profiles)[1]
  K <- dim(normal_profiles)[2]
  for(arsd in c(0.5, 1, 2)){
    for(rho in c(0, 0.2, 0.4)){
      for(rep in start.rep:stop.rep){
        file.prefix <- paste0("ar_sd", arsd, "_rho_", rho, "_n", rep)
        data.file <- paste0("data/", file.prefix, "_data.RData")
        y <- matrix(nrow=p, ncol=K)
        full.data <- matrix(nrow=p, ncol=0)
        for(i in 1:K){
          my.y <- replicate(n=sample.size[i],
                            expr=normal_data(normal_profiles[,i], arsd, rho, 0))
          full.data <- cbind(full.data, my.y)
        }
        R <- list("Y"=full.data, "sample.size"=sample.size)
        save(R, file=data.file)
        alternatives_normal(file.prefix=file.prefix)
      }
    }
  }
}

#'@title Wrapper to generate normal data with random effects
#'@description Generate normal data with random effects. Run form a directory with a "data/" subdirectory.
#'@export
generate_data_re <- function(seed=3333333, start.rep=1, stop.rep=60){
  set.seed(seed)
  sample.size=c(10, 10)
  data("normal_profiles", package="jadeSims")
  p <- dim(normal_profiles)[1]
  K <- dim(normal_profiles)[2]

  sds <- c(2.18,2.12,  2.06, 2)
  res <- c(0.5, 0.707, .866, 1)
  perc <- c(5, 10, 15, 20)
  for(l in 1:4){
    for(rep in start.rep:stop.rep){
      file.prefix <- paste0("re_", perc[l], "_n", rep)
      data.file <- paste0("data/", file.prefix, "_data.RData")
      y <- matrix(nrow=p, ncol=K)
      full.data <- matrix(nrow=p, ncol=0)
      for(i in 1:K){
        my.y <- replicate(n=sample.size[i],
                      expr=normal_data(normal_profiles[,i], sds[l], 0, res[l]))
        full.data <- cbind(full.data, my.y)
      }
      R <- list("Y"=full.data, "sample.size"=sample.size)
      save(R, file=data.file)
      alternatives_normal(file.prefix=file.prefix)
    }
  }
}



#'@title Wrapper to generate normal auto-regressive data with long profiles
#'@description Generate normal AR  data. Run form a directory with a "data/" subdirectory.
#'@export
generate_data_long <- function(sig, rho, re, start.rep, stop.rep, sample.size=c(10, 10), seed=NULL){
  if(!is.null(seed))set.seed(seed)
  p <- dim(long_profiles3)[1]
  K <- dim(long_profiles3)[2]

  for(rep in start.rep:stop.rep){
    file.prefix <- paste0("ar_long_sd_", sig, "_rho_", rho,"_re_", re,  "_n", rep)
    data.file <- paste0("data/", file.prefix, "_data.RData")
    y <- matrix(nrow=p, ncol=K)
    full.data <- matrix(nrow=p, ncol=0)
    for(i in 1:K){
      my.y <- replicate(n=sample.size[i],
                expr=normal_data(long_profiles3[,i], sig, rho, re))
      full.data <- cbind(full.data, my.y)
    }
    R <- list("Y"=full.data, "sample.size"=sample.size)
    save(R, file=data.file)
    alternatives_normal(file.prefix = file.prefix)
  }
}

#'@title Wrapper to generate normal auto-regressive data with long profiles
#'@description Generate normal AR  data. Run form a directory with a "data/" subdirectory.
#'@export
generate_data_long2 <- function(sig, rho, re, start.rep, stop.rep,
                                file.prefix, bandwidth=20,
                                sample.size=c(10, 10), seed=NULL){
  if(!is.null(seed))set.seed(seed)
  p <- dim(long_profiles3)[1]
  K <- dim(long_profiles3)[2]
  stopifnot(K==2)
  stopifnot(length(sig)==length(sample.size))
  labs <- rep(c(0, 1), sample.size)
  for(rep in start.rep:stop.rep){
    data.file <- paste0("data/", file.prefix, "_n", rep, "_data.RData")
    y <- matrix(nrow=p, ncol=K)
    full.data <- matrix(nrow=p, ncol=0)
    for(i in 1:K){
      my.y <- sapply(sig[[i]], FUN=function(s){
          normal_data(long_profiles3[,i], s, rho, re)})
      full.data <- cbind(full.data, my.y)
    }

    #dat.sm <- apply(full.data, MARGIN=2, FUN=function(y){
    #  ksmooth(1:500, y, x.points = 1:500, bandwidth = bandwidth)$y
    #})
    #sm.tests <- cfdrSims:::t_stats(dat.sm, labs, s0=0)
    #tests <- cfdrSims:::t_stats(full.data, labs=labs, s0=0)
    #tests.sm <- ksmooth(1:500, tests, bandwidth=bandwidth, x.points=1:500)$y

    R <- list("Y"=full.data, "sample.size"=sample.size)
    save(R, file=data.file)
    alternatives_normal(file.prefix = file.prefix)
  }
}




#'@title Wrapper to generate normal data with short profiles
#'@description Generate normal data. Run form a directory with a "data/" subdirectory.
#'@export
generate_data_short <- function(sigma, rho, re, start.rep, stop.rep,prefix, sample.size){
  #set.seed(2222222)
  p <- dim(normal_profiles)[1]
  K <- dim(normal_profiles)[2]

  for(rep in start.rep:stop.rep){
    file.prefix <- paste0(prefix, "_sd_", sigma, "_rho_", rho, "_re_", re, "_n", rep)
    data.file <- paste0("data/", file.prefix, "_data.RData")
    y <- matrix(nrow=p, ncol=K)
    full.data <- matrix(nrow=p, ncol=0)
    for(i in 1:K){
      my.y <- replicate(n=sample.size[i],
                        expr=normal_data(normal_profiles[,i], sigma, rho, re))
      full.data <- cbind(full.data, my.y)
    }
    R <- list("Y"=full.data, "sample.size"=sample.size)
    save(R, file=data.file)
    alternatives_normal(file.prefix = file.prefix)
  }
}
