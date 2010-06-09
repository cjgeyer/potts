\name{calc_t}
\alias{calc_t}
\alias{calc_t_innergrid}
\alias{calc_t_full}
\title{Calculate Canonical Statistic for Potts Model}
\description{
  Caclulate the canonical statistic 't' for a realization of a Potts
  Model
}
\usage{
calc_t_full(x,ncolor)
calc_t_innergrid(x, ncolor, grid, i, j)
calc_t(x, ncolor, grid=NULL, i=NULL, j=NULL)
}
\arguments{
  
  \item{x}{2 dimensional matrix, elements in 1, \dots, ncolor,
           representing a Potts model}
  
  \item{ncolor}{Number of colors in this Potts Model.}

  \item{grid}{2 dimensional matrix, elements in 1, \dots, ncolor.  If
    non-NULL it is placed into x at the location x[i,j].}

  \item{i}{row to place the grid.}

  \item{j}{column to place the grid.}

	       
}
\details{
  
  Calculates the canonical statistics for a realized Potts Model.
  \code{calc_t} calls \code{calc_t_full} if \code{grid} is \code{NULL}
  and \code{calc_t_innergrid} otherwise.

  \code{calc_t_full} calculates the canonical statistics for the full image.

  \code{calc_t_innergrid} calculates the canonical statistics for the a
  window of the image, but with that window replaced by \code{grid},
  with the upper left corner of \code{grid} located at \code{x[i,j]}.
  
}
\value{
  All functions return a vector of length \code{ncolor+1}.  Elements
  1,\dots,ncolor contain the number of pixels of each color.  Element
  ncolor+1 contains the number of matching neighbor pairs for the image.

  \code{calc_t_full} returns the values for the whole image.

  \code{calc_t_innergrid} returns the value for just the selected
  window, but this includes the number of matching pairs from the
  replaced window to it's neighbors.  E.g. if \code{X} is the full
  image, and A_1 is the value of some window in the image and we want to
  know the value of \code{t( y union X\\A_1)} this would be
  \code{calc_t_full(X, ncolor) + calc_t_innergrid(X, ncolor, y, i, j) -
  calc_t_innergrid(X, ncolor, A_1, i, j)}
}
\examples{

ncolor <- 4
beta   <- log(1+sqrt(ncolor))
theta  <- c(rep(0,ncolor), beta)

nrow <- 32
ncol <- 32

x <- matrix(sample(ncolor, nrow*ncol, replace=TRUE), nrow=nrow, ncol=ncol)
foo <- packPotts(x, ncolor)
out <- potts(foo, theta, nbatch=10)
x <- unpackPotts(out$final)

t_stat <- calc_t(x, ncolor)[-1]
t_stat_inner <- calc_t(x, ncolor, matrix(1, nrow=2, ncol=2), 1, 1)

}
\keyword{misc}
