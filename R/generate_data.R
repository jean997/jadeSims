

#'Simulate normal data
#'@description Simulate normal data
#'with auto-regressive errors or random effects
#'@param profile Mean values (length p)
#'@param sd Standard deviation
#'@param rho Auto-regression parameter
#'@param reff.sd Standard deviation of random effects
#'@return A length p vector of simulated data points.
#'@export
normal_data <- function(profile, sd, rho, reff.sd){
	#if(rho > 0 & reff.sd > 0) stop("Are you sure you wanted rho and reff.sd?")
	p <- length(profile)
	y <- profile+rnorm(1, 0, reff.sd) + filter(rnorm(p, mean=0, sd=sd), filter=rho, method="recursive")
	return(y)
}


#'Simulate binomial data
#'@description Simulate binomial data
#'by tiling reads.
#'@param profile Mean values (length p)
#'@param mean.reads Desired average number of reads per site.
#'@param reff.sd Standard deviation of random effects
#'@return A data frame with columns y and reads
#'@export
tile_binomial_data <- function(profile, mean.reads, reff.sd, mean.tile.length=30){
  p <- length(profile)
  mt <- tile_reads(nbp=p, avg_length=mean.tile.length, avg_reads=mean.reads)
  if(reff.sd > 0) profile <- shift_profile(profile, reff.sd)
  y <- rbinom(n=p, size=mt$coverage, prob=profile)
  R <- data.frame("y"=y, "reads"=mt$coverage)
  return(R)
}


tile_reads <- function(nbp, avg_length, avg_reads){
  L <- nbp+avg_length
  #number of tiles we will need
  ntiles <- L*avg_reads/avg_length
  starts <- sample(-avg_length:nbp, size=ntiles, replace=TRUE)
  lengths <- floor(rexp(n=ntiles, rate=1/avg_length))
  starts <- starts[lengths > 0]
  lengths <- lengths[lengths > 0]

  tiles <- data.frame("start"=starts, "length"=lengths)
  tiles <- tiles[order(starts),]
  #Check for total coverage
  sites_bp <- 1:nbp
  coverage_total <- rep(0, nbp)
  for(t in 1:nrow(tiles)){
    which_tile <- which( sites_bp >= tiles$start[t] & sites_bp <= tiles$start[t]+tiles$length[t])
    if(length(which_tile)==0) next
    coverage_total[which_tile] <- coverage_total[which_tile]+1
  }
  zeros <- which(coverage_total==0)
  #prevent sites with no coverage
  while(length(zeros) > 0){
    new.length <- floor(rexp(n=1, rate=1/avg_length))
    new.start <- sample((zeros[1]-new.length+1):zeros[1], size=1, replace=TRUE)
    tiles <- rbind(tiles, c(new.start, new.length))
    t <- t+1
    which_tile <- which( sites_bp >= tiles$start[t] & sites_bp <= tiles$start[t]+tiles$length[t])
    coverage_total[which_tile] <- coverage_total[which_tile]+1
    zeros <- which(coverage_total==0)
  }
  return(list("tiles"=tiles, "coverage"=coverage_total))
}

shift_profile <- function(profile, reff.sd){
	re <- rnorm(n=1, mean=0, sd=reff.sd)
	profile <- profile + re
	profile <- pmin(1, profile)
	profile <- pmax(0, profile)
	return(profile)
}

