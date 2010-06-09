##############################################################################
#
# This file contains various functions used for calculating t_cache
#
##############################################################################

##############################################################################
#
# This calculates our t_cache, it needs x, ncolor, t_stat, sizeA,
# npixel, f, and gridcache
#
# f is a function x*ncolor*A*I:-> t() .  If the index is 0, it means to
# calculate the observed for that cell.
#
# gridcache is any cache that f needs.  It may be null.
#
##############################################################################
generate_t_cache <- function(x, ncolor, t_stat, sizeA, npixel, f,
                             fapply=lapply, gridcache=NULL) {
  nIdx <- ncolor**npixel
  fapply(1:sizeA, function(a) {
    t_start <- t_stat - f(x, ncolor, a, 0, gridcache)
    arr <- array(0, dim=c(nIdx, ncolor+1-1))
    for(idx in 1:nIdx) {
      arr[idx,] <- t_start + f(x, ncolor, a, idx, gridcache)
    }
    arr
  })
}

##############################################################################
#
# This function is used in generating the grids for bigger window sizes.
# it cycles a vector through the set (1:ncolor)**windowsize
#
##############################################################################
incvector <- function(v, ncolor) {
  vorig <- v
  i <- length(v)
  v[i] <- v[i] + 1
  overflow <- FALSE
  while(i > 0 && v[i] > ncolor) {
    v[i] <- 1
    i <- i - 1
    if (i > 0)
      v[i] <- v[i] + 1
    else
      overflow <- TRUE
  }
  if (overflow)
    vorig
  else
    v
}

##############################################################################
#
# generate a grid cache
#
##############################################################################
gengridcache <- function(ncolor, v, ncol) {
  n <- ncolor**length(v)
  l <- list()
  for(i in 1:n) {
    l[[i]] <- matrix(v, ncol=ncol)
    v <- incvector(v, ncolor)
  }
  l
}
##############################################################################
#
# generate C^A for a single pixel window
#
##############################################################################
gensingleton <- function(ncolor) {
  gengridcache(ncolor, c(1), 1)
}

singleton <- function(x, ncolor, a, idx, gcache=NULL) {
  if (is.null(gcache)) {
    # singleton is so small we don't bother with cacheing
    gcache <- gensingleton(ncolor)
  }
  ncol <- dim(x)[2]
  # get the appropriate grid to put in
  i <- ceiling(a/ncol)
  j <- a - (i - 1) * ncol
  grid <- matrix(x[i,j], ncol=1)
  if (idx > 0) {
    grid <- gcache[[idx]]
  }
  calc_t(x, ncolor, grid=grid, i,j)[-1]
}

##############################################################################
#
# generate the elements of C^A for 2 pixel horizontal grids.
# overlapping and non-overlapping version
#
##############################################################################
gentwopixel <- function(ncolor) {
  gengridcache(ncolor, rep(1,2), 2)
}

twopixel <- function(x, ncolor, a, idx, gcache=NULL) {
  if (is.null(gcache)) {
    # no cache for two pixel
    gcache <- gentwopixel(ncolor)
  }
  ncol <- dim(x)[2]
  # get the appropriate grid to put in
  i <- ceiling(a/ncol)
  j <- a - (i - 1) * ncol
  j2 <- ifelse(j==ncol, 1, j+1)
  grid <- matrix(c(x[i,j], x[i,j2]), ncol=2)
  if (idx > 0) {
    grid <- gcache[[idx]]
  }
  calc_t(x, ncolor, grid=grid, i,j)[-1]
}

twopixel.nonoverlap <- function(x, ncolor, a, idx, gcache=NULL) {
  if (is.null(gcache)) {
    # no cache for two pixel
    gcache <- gentwopixel(ncolor)
  }
  ncol <- dim(x)[2]
  nc2 <- ncol/2
  i <- ceiling(a / nc2)
  j <- a - (i - 1) * nc2
  j <- (j-1)*2+1
  j2 <- j + 1
  grid <- matrix(c(x[i,j], x[i,j2]), ncol=2)
  if (idx > 0) {
    grid <- gcache[[idx]]
  }
  calc_t(x, ncolor, grid=grid, i,j)[-1]
}

##############################################################################
#
# generate C^A for 2x2 pixel windows.  Overlapping and non-overlapping
# versions.
#
#
##############################################################################
genfourpixel <- function(ncolor) {
  gengridcache(ncolor, rep(1,4), 2)
}

fourpixel <- function(x, ncolor, a, idx, gcache=NULL) {
  if (is.null(gcache)) {
    # no cache for four pixel
    gcache <- genfourpixel(ncolor)
  }
  nrow <- dim(x)[1]
  ncol <- dim(x)[2]
  # get the appropriate grid to put in
  i <- ceiling(a/ncol)
  j <- a - (i - 1) * ncol
  i2 <- ifelse(i == nrow, 1, i+1)
  j2 <- ifelse(j == ncol, 1, j+1)
  grid <- x[c(i,i2),c(j,j2)]
  if (idx > 0) {
    grid <- gcache[[idx]]
  }
  calc_t(x, ncolor, grid=grid, i,j)[-1]
}

fourpixel.nonoverlap <- function(x, ncolor, a, idx, gcache=NULL) {
  if (is.null(gcache)) {
    gcache <- genfourpixel(ncolor)
  }
  nrow <- dim(x)[1]
  ncol <- dim(x)[2]
  nr2 <- nrow/2
  nc2 <- ncol/2
  i <- ceiling(a / nc2)
  j <- a - (i - 1) * nc2
  i <- (i-1)*2+1
  j <- (j-1)*2+1
  i2 <- i+1
  j2 <- j+1
  grid <- x[c(i,i2),c(j,j2)]
  if (idx > 0) {
    grid <- gcache[[idx]]
  }
  calc_t(x, ncolor, grid=grid, i,j)[-1]
}
##############################################################################
#
# generate a list of the elements of C^A for a 3x3 window
#
##############################################################################
genthreebythree <- function(ncolor) {
  gengridcache(ncolor, rep(1,9), 3)
}

nine.cache <- list()

ninepixel.nonoverlap <- function(x, ncolor, a, idx, gcache=NULL) {
  if (is.null(gcache)) {
    # generate this windows C^A if we don' have it cached.
    if (length(nine.cache) < ncolor ||
        is.null(nine.cache[[ncolor]])) {
      nine.cache[[ncolor]] <- genthreebythree(ncolor)
    }
    gcache <- nine.cache[[ncolor]]
  }
  nrow <- dim(x)[1]
  ncol <- dim(x)[2]
  nboxrow <- ceiling(ncol/3)
  nboxcol <- ceiling(nrow/3)
  # boxes go in rows first, then columns, so determine the row number
  rnum <- ceiling(a/nboxrow)
  cnum <- a - nboxrow * (rnum-1)
  i <- (rnum - 1) * 3 + 1:3
  j <- (cnum - 1) * 3 + 1:3
  i[i > nrow] <- i[i > nrow] - nrow
  j[j > ncol] <- j[j > ncol] - ncol
  grid <- NULL 
  if (idx > 0) {
    grid <- gcache[[idx]]
  } else {
    grid <- x[i,j]
  }
  calc_t(x, ncolor, grid=grid, i[1], j[1])[-1]
}
##############################################################################
#
# Functions for 4x4 windows.
#
##############################################################################
genfourbyfour <- function(ncolor) {
  gengridcache(ncolor, rep(1,16), 4)
}

sixteen.cache <- list()

sixteenpixel.nonoverlap <- function(x, ncolor, a, idx, gcache=NULL) {
  if (is.null(gcache)) {
    if (length(sixteen.cache) < ncolor ||
        is.null(sixteen.cache[[ncolor]])) {
      sixteen.cache[[ncolor]] <<- genfourbyfour(ncolor)
    }
    gcache <- sixteen.cache[[ncolor]]
  }
  nrow <- dim(x)[1]
  ncol <- dim(x)[2]
  nboxrow <- ceiling(ncol/4)
  nboxcol <- ceiling(nrow/4)
  # boxes go rows firts, then columns, so determine the row number
  rnum <- ceiling(a/nboxrow)
  cnum <- a - nboxrow * (rnum-1)
  i <- (rnum - 1) * 4 + 1:4
  j <- (cnum - 1) * 4 + 1:4
  i[i > nrow] <- i[i > nrow] - nrow
  j[j > ncol] <- j[j > ncol] - ncol
  grid <- NULL 
  if (idx > 0) {
    grid <- gcache[[idx]]
  } else {
    grid <- x[i,j]
  }
  calc_t(x, ncolor, grid=grid, i[1], j[1])[-1]
}

