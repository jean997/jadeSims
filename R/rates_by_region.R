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
    ix <- which(ov ==0)
    for(j in 1:ix)
      xj <- xI[j,]
      ixt <- interval_intersection(xj, labI)
      ixt_l <- ixt[,2]-ixt[,1] +1
      dif <- interval_difference(xj, labI)
      dif_l <- dif[,2]-dif[,1]
      acc[ct] <- ixt_l/(dif_l + ixt_l)
  }
  fpr <- fp.ct /(fp.ct + tp.ct)
  return(list("fp.ct"=fp.ct, "fpr"=fpr, "tp.ct"=tp.ct, "acc"=acc, "disc"=disc))
}

get_region_rates <- function(prefix, merge.margin=0, min.length=1, min.acc=0){

}


get_regions <- function(x, min.length=3, merge.margin=0){
  q0 <- matrix(nrow=0, ncol=2)
  ind <- x==1
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
