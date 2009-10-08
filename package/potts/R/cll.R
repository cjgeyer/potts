##############################################################################
#
# calculate composite log likelihood for a potts model.
#
##############################################################################
composite.ll <- function(theta, t_cache=NULL) {
  tot <- 0
  if (is.null(t_cache)) {
    stop("!cache_t not implemented yet in composite.ll!")
  } else {
    # this could probably be re-written to be an apply, but I was
    # having trouble with the indicies
    for(j in 1:( dim(t_cache)[1] )) { # for A in script(A)
      arr <- t_cache[j,,]
      # subtract base case.
      tmp <- t(apply(arr, 1, function(r) r-arr[1,]))
      tot <- tot + (t_stat-arr[1,]) %*% theta -
        log(sum(exp(tmp %*% theta)))
    }
  }
  tot
}

##############################################################################
#
# calculate gradient of the composite log likelihood for a potts
# model.
#
##############################################################################
gr.composite.ll <- function(theta, t_cache=NULL) {
  tot <- rep(0,length(theta))
  if (is.null(t_cache)) {
    stop("!cache_t not implemented yet in composite.ll!")
  } else {
    # this could probably be cleverly rewritten to be an apply
    for(i in 1:( dim(t_cache)[1] )) {
      arr <- t_cache[i,,]
      tot <- tot + t_stat - arr[1,]
      tmp <- rowSums(apply(arr, 1, function(r) {
        rtmp <- exp( (r-arr[1,]) %*% theta)
        c( (r-arr[1,]) * rtmp, rtmp)
      }))
      tot <- tot - tmp[1:(length(tmp)-1)]/tmp[length(tmp)]
    }
  }
  tot
}
