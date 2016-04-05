#'@title Wrapper to generate binomial data
#'@description Generate binomial data. Run form a directory with a "data/" subdirectory.
#'@export
generate_data_binomial <- function(){
  set.seed(1111111)
  sample.size=c(10, 10)
  data("binom_profiles", package="jadeSims")
  p <- dim(binom_profiles)[1]
  K <- dim(binom_profiles)[2]
  for(re in c(0, 0.02, 0.05, 0.07)){
    for(rep in 1:60){
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
    }
  }
}

#'@title Wrapper to generate normal auto-regressive data
#'@description Generate normal AR  data. Run form a directory with a "data/" subdirectory.
#'@export
generate_data_ar <- function(){
  set.seed(2222222)
  sample.size=c(10, 10)
  data("normal_profiles", package="jadeSims")
  p <- dim(normal_profiles)[1]
  K <- dim(normal_profiles)[2]
  for(arsd in c(0.5, 1, 2)){
    for(rho in c(0, 0.2, 0.4)){
      for(rep in 1:60){
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
      }
    }
  }
}

#'@title Wrapper to generate normal data with random effects
#'@description Generate normal data with random effects. Run form a directory with a "data/" subdirectory.
#'@export
generate_data_re <- function(){
  set.seed(3333333)
  sample.size=c(10, 10)
  data("normal_profiles", package="jadeSims")
  p <- dim(normal_profiles)[1]
  K <- dim(normal_profiles)[2]

  sds <- c(2.18,2.12,  2.06, 2)
  res <- c(0.5, 0.707, .866, 1)
  perc <- c(5, 10, 15, 20)
  for(l in 1:4){
    for(rep in 1:60){
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
    }
  }
}



#'@title Wrapper to generate normal auto-regressive data with long profiles
#'@description Generate normal AR  data. Run form a directory with a "data/" subdirectory.
#'@export
generate_data_long <- function(arsd, rho, n.rep, sample.size=c(20, 20)){
  set.seed(2222222)
  p <- dim(long_profiles2)[1]
  K <- dim(long_profiles2)[2]

  for(rep in 1:n.rep){
    file.prefix <- paste0("ar_long_sd_", arsd, "_rho_", rho, "_n", rep)
    data.file <- paste0("data/", file.prefix, "_data.RData")
    y <- matrix(nrow=p, ncol=K)
    full.data <- matrix(nrow=p, ncol=0)
    for(i in 1:K){
      my.y <- replicate(n=sample.size[i],
                expr=normal_data(long_profiles2[,i], arsd, rho, 0))
      full.data <- cbind(full.data, my.y)
    }
    R <- list("Y"=full.data, "sample.size"=sample.size)
    save(R, file=data.file)
    alternatives_normal(file.prefix = file.prefix)
  }
}


#'@title Wrapper to generate normal data with short profiles
#'@description Generate normal data. Run form a directory with a "data/" subdirectory.
#'@export
generate_data_short <- function(sigma, rho, re, n.rep,prefix, sample.size){
  #set.seed(2222222)
  p <- dim(normal_profiles)[1]
  K <- dim(normal_profiles)[2]

  for(rep in 1:n.rep){
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


#'@title Wrapper to generate binomial data based on methylation
#'@description Generate binomial data.
#'Run form a directory with a "data/" subdirectory.
#'@export
generate_data_meth <- function(prefix, re, sample.size,
                               n.rep, rd.factor=1, seed=NULL){
  if(!is.null(seed)) set.seed(seed)

  p <- dim(binom_meth_profiles)[1]
  K <- dim(binom_meth_profiles)[2]

  for(rep in 1:n.rep){
    file.prefix <- paste0(prefix, "_re_", re, "_n", rep)
    data.file <- paste0("data/", file.prefix, "_data.RData")
    full.counts <- full.reads <-  matrix(nrow=p, ncol=0)
    for(i in 1:K){
      done <- FALSE
      while(!done){
        M <- replicate(n=sample.size[i], expr={
          rds <- binom_meth_reads[, sample(size = 1, x=1:8)]
          rds <- ceiling(rd.factor*rds)
          cbind(binomial_data(binom_meth_profiles[,i], rds, re)$y, rds)
        })
        if(all(rowSums(M[,2,])  > 0)) done <- TRUE
      }
      full.counts <- cbind(full.counts, M[,1, ])
      full.reads <- cbind(full.reads, M[,2, ])
    }
    R <- list("Y"=full.counts, "READS"=full.reads,
              "sample.size"=sample.size)
    save(R, file=data.file)
    alternatives_binomial(file.prefix=file.prefix, bs.factor=1, bs.ns=40,
                          sites=binom_meth_pos)
  }
}
