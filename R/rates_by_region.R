library(intervals)
rates_by_region <- function(x, labels,
                            merge.margin=0, min.length=1){
  labI <- Intervals(get_regions(labels, merge.margin=0, min.length=1))
  xI <- Intervals(get_regions(x, merge.margin=merge.margin, min.length=min.length))

  ov_X<- distance_to_nearest(xI, labI)

  ov_L <- distance_to_nearest(labI, xI)
  disc <- c(0, 0)
  disc[ov_L==0] <- 1

  tp.ct <- sum(ov_X ==0)
  acc <- rep(0, length(tp.ct))
  fp.ct <- sum(ov_X > 0)
  ct <- 1
  if(any(ov_X==0)){
    ix <- which(ov_X ==0)
    for(j in ix)
      xj <- xI[j,]
      ixt <- interval_intersection(xj, labI)
      if(nrow(ixt)==0) ixt_l <- 0
        else ixt_l <- ixt[,2]-ixt[,1] +1
      dif <- interval_difference(xj, labI)
      if(nrow(dif)==0) dif_l <- 0
        else dif_l <- dif[,2]-dif[,1]
      acc[ct] <- ixt_l/(dif_l + ixt_l)
      ct <- ct + 1
  }
  fpr <- fp.ct /(fp.ct + tp.ct)
  p <- length(x)
  tot.disc <- sum(size(xI))/p
  return(list("fp.ct"=fp.ct, "fpr"=fpr, "tpr"=sum(disc)/2, "tp.ct"=tp.ct,
              "acc"=acc, "disc"=disc, "tot.disc"=tot.disc))
}


#'@title Get average region level tpr and fpr from and agg object
#'@param agg.obj An object produced by aggregate_sims
#'@param stat.names Which statistics to compare to
#'@param merge.margin Merge regions separated by merge.margin
#'@param min.length Only consider regions longer than min.length
#'@param min.acc Not currently used
#'@return A list of objects produced by get_tpr_fpr
#'@export
get_region_rates <- function(agg.obj,
                            stat.names,
                             merge.margin=0, min.length=1, min.acc=0){
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
  fpr.list <- list()
  prop.list <- list()
  cat("JADE\n")
  for(j in 1:N){
    cat(j, " ")
    jade.rates <-apply(agg.obj$all.sep[[j]], MARGIN=2, FUN=function(x){
      unlist(rates_by_region(x, labels, merge.margin, min.length))
    })
    tpr.list[[j]] <- jade.rates["tpr",]
    fpr.list[[j]] <- jade.rates["fpr",]
    prop.list[[j]] <- jade.rates["tot.disc",]
  }
  cat("\n")
  rate.list[[ct]] <- avg_by_prop(tpr.list, fpr.list, prop.list)
  ct <- ct + 1

  for(ix in which.stats){
    cat(agg.obj$names[ix], "\n")
    tpr.list <- list()
    fpr.list <- list()
    prop.list <- list()
    for(j in 1:N){
      cat(j, " ")
      sort.p <- sort(agg.obj$all.stats[j, ,ix])
      M <- sapply(1:p, FUN=function(t, stats, sort.p){
        x <- stats < sort.p[t]
        unlist(rates_by_region(x, labels, merge.margin, min.length))
      }, stats=agg.obj$all.stats[j, , ix], sort.p=sort.p)
      tpr.list[[j]] <- M["tpr",]
      fpr.list[[j]] <- M["fpr",]
      prop.list[[j]] <- M["tot.disc",]
    }
    rate.list[[ct]] <- avg_by_prop(tpr.list, fpr.list, prop.list)
    ct <- ct + 1
    cat("\n")
  }
  return(rate.list)
}

plot_rates <- function(rate.list, cols){
  plot(rate.list[[1]]$fpr, rate.list[[1]]$tpr, type="l",
       ylim=c(0, 1), xlim=c(0, 1), col=cols[1])
  for(i in 2:length(rate.list)){
    lines(rate.list[[i]]$fpr, rate.list[[i]]$tpr, col=cols[i])
  }

  plot(rate.list[[1]]$prop, rate.list[[1]]$tpr, type="l",
        ylim=c(0, 1), xlim=c(0, 1), col=cols[1])
  for(i in 2:length(rate.list)){
    lines(rate.list[[i]]$prop, rate.list[[i]]$tpr, col=cols[i])
  }
  plot(rate.list[[1]]$prop, rate.list[[1]]$fpr, type="l",
        ylim=c(0, 1), xlim=c(0, 1), col=cols[1])
  for(i in 2:length(rate.list)){
    lines(rate.list[[i]]$prop, rate.list[[i]]$fpr, col=cols[i])
  }

}
get_regions <- function(x, min.length=3, merge.margin=0){
  q0 <- matrix(nrow=0, ncol=2)
  if(class(x)=="logical") ind <- x
    else ind <- x==1
  y <- rle(ind)
  n <- length(y$values)
  #Merge first
  starts <- cumsum(c(1, y$lengths[-n]))
  stops <- cumsum(y$lengths)
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
  #Then threshold
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
  idx <- which(y$values==TRUE)
  starts <- st.copy <- starts[idx]
  stops <- sp.copy <-  stops[idx]
  q0 <- cbind(starts, stops)

  #q0I <- Intervals(q0)
  return(q0)
}
