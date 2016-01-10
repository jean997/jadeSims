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

generate_data_ar <- function(){
  set.seed(2222222)
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

generate_data_re <- function(){
  set.seed(3333333)
  data("normal_profiles", package="jadeSims")
  p <- dim(normal_profiles)[1]
  K <- dim(normal_profiles)[2]

  sds <- c(2.18,2.12,  2.06, 2)
  res <- c(0.5, 0.707, .866, 1)
  perc <- c(5, 10, 15, 20)
  for(l in 1:4){
    for(rep in 1:60){
      file.prefix <- paste0("re_", perc, "_n", rep)
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

