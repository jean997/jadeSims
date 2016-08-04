
rates_by_region <- function(x, labels,
                            merge.margin=0, min.length=1){
  labI <- Intervals(get_regions(labels, merge.margin=0, min.length=1), type="Z")
  xI <- Intervals(get_regions(x, merge.margin=merge.margin, min.length=min.length), type="Z")

  if(nrow(xI)==0){
    return(c("fp.ct"=0, "fdp"=0, "tpr"=0, "tp.ct"=0,
             "tot.disc"=0))
  }

  ov_X<- distance_to_nearest(xI, labI)

  ov_L <- distance_to_nearest(labI, xI)

  tp.ct <- sum(ov_X ==0)
  acc <- rep(0, length(ov_X))
  fp.ct <- sum(ov_X > 0)

  if(any(ov_X==0)){
    ix <- which(ov_X ==0)
    for(j in ix){
      xj <- xI[j,]
      ixt <- interval_intersection(xj, labI)
      acc[j] <- sum(size(ixt))/size(xj)
    }
  }

  xTP <- xI[ov_X==0 & acc >= min.acc,]
  ov_LTP <- distance_to_nearest(labI, xTP)
  disc <- rep(0, nrow(labI))
  disc[ov_LTP==0] <- 1


  fdp <- fp.ct /(fp.ct + tp.ct)
  p <- length(x)
  tot.disc <- sum(size(xI))/p
  return(c("fp.ct"=fp.ct, "fdp"=fdp, "tpr"=sum(disc)/nrow(labI), "tp.ct"=tp.ct,
           "tot.disc"=tot.disc))
}


#'@title Get average region level tpr and fpr from and agg object
#'@param agg.obj An object produced by aggregate_sims
#'@param stat.names Which statistics to compare to
#'@param merge.margin Merge regions separated by merge.margin
#'@param min.length Only consider regions longer than min.length
#'@param min.acc Minimum accuract for a tp to count
#'@return A list of objects produced by get_tpr_fpr
#'@export
get_region_rates <- function(agg.obj,
                            stat.names, max.prop=0.5,
                             merge.margin=0, min.length=1){
  #For each stat at each p-value level - p x length(which.stats) x 5
  #For jade at each level of gamma ngamma x 5

  labels <- agg.obj$site.labels
  N <- dim(agg.obj$all.stats)[1]
  r <- length(stat.names)
  which.stats <- match(stat.names, agg.obj$names)
  p <- dim(agg.obj$all.stats)[2]

  rate.list <- list()
  ct <- 1

  tpr.list <- list()
  fdp.list <- list()
  fct.list <- list()
  prop.list <- list()
  cat("JADE\n")
  #rm.idx <- c()
  for(j in 1:N){
    cat(j, " ")
    jade.rates <-apply(agg.obj$all.sep[[j]], MARGIN=2, FUN=function(x){
      jadeSims:::rates_by_region(x, labels, merge.margin, min.length)
    })
    midx <- which(jade.rates["tot.disc",] <= max.prop)
    cat(length(midx), "\n")
    tpr.list[[j]] <- jade.rates["tpr", midx]
    fdp.list[[j]] <- jade.rates["fdp",midx]
    fct.list[[j]] <- jade.rates["fp.ct", midx]
    prop.list[[j]] <- jade.rates["tot.disc", midx]
    #if(all(fct.list[[j]]==0)) rm.idx <- c(rm.idx, j)
  }
  cat("\n")
  #if(avg.by.prop) rate.list[[ct]] <- avg_by_prop(tpr.list, fpr.list, prop.list)
   # else if(length(rm.idx) > 0) rate.list[[ct]] <- avg_by_interp(tpr.list[-rm.idx], fpr.list[-rm.idx])
  #    else rate.list[[ct]] <- avg_by_interp(tpr.list, fpr.list)
  rate.list[[ct]] <- avg_by_ct(tpr.list, fct.list)
  ct <- ct + 1

  for(ix in which.stats){
    cat(agg.obj$names[ix], "\n")
    tpr.list <- list()
    fpr.list <- list()
    prop.list <- list()
    for(j in 1:N){
      cat(j, " ")
      q <- quantile(agg.obj$all.stats[j, , ix], probs=max.prop)
      sort.p <- sort(agg.obj$all.stats[j, ,ix])
      sort.p <- sort.p[sort.p <= q]
      M <- sapply(sort.p, FUN=function(t, stats){
        x <- stats < t
        unlist(jadeSims:::rates_by_region(x, labels, merge.margin, min.length))
      }, stats=agg.obj$all.stats[j, , ix])
      midx <- which(M["tot.disc",] <= max.prop)
      cat(length(midx), " ")
      tpr.list[[j]] <- M["tpr", midx]
      fdp.list[[j]] <- M["fdp",midx]
      fct.list[[j]] <- M["fp.ct", midx]
      prop.list[[j]] <- M["tot.disc", midx]
    }
    rate.list[[ct]] <- avg_by_ct(tpr.list, fct.list)
    ct <- ct + 1
    cat("\n")
  }
  rate.list$names <- c("jade", stat.names)
  return(rate.list)
}



get_regions <- function(x, min.length=3, merge.margin=0){
  q0 <- matrix(nrow=0, ncol=2)
  if(class(x)=="logical") ind <- x
    else ind <- x==1
  y <- rle(ind)
  n <- length(y$values)
  starts <- cumsum(c(1, y$lengths[-n]))
  stops <- cumsum(y$lengths)
  #Threshold first
  if(min.length > 1){
    idx <- which(y$lengths < min.length & y$values ==TRUE)
    for(j in idx){
      ind[starts[j]:stops[j]] <- FALSE
    }
    y <- rle(ind)
    n <- length(y$values)
    starts <- cumsum(c(1, y$lengths[-n]))
    stops <- cumsum(y$lengths)
  }


  #Then Merge
  if(merge.margin > 0){
    idx <- which(y$lengths <= merge.margin & (y$values == FALSE| is.na(y$values)))
    idx <- idx[!idx ==1]
    idx <- idx[!idx ==n]
    for(j in idx){
      ind[starts[j]:stops[j]] <- TRUE
    }
    y <- rle(ind)
    n <- length(y$values)
    starts <- cumsum(c(1, y$lengths[-n]))
    stops <- cumsum(y$lengths)
  }

  idx <- which(y$values==TRUE)
  starts <- starts[idx]
  stops <- stops[idx]
  q0 <- cbind(starts, stops)

  #q0I <- Intervals(q0)
  return(q0)
}
