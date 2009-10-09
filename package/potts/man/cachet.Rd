\name{Cache T}
\alias{generate_t_cache}
\alias{incvector}
\alias{gengridcache}
\alias{gensingleton}
\alias{singleton}
\alias{gentwopixel}
\alias{twopixel}
\alias{twopixel.nonoverlap}
\alias{genfourpixel}
\alias{fourpixel}
\alias{fourpixel.nonoverlap}
\alias{genthreebythree}
\alias{ninepixel.nonoverlap}
\alias{genfourbyfour}
\alias{sixteenpixel.nonoverlap}
\title{Cache calculated Canonical Statistics for Potts Models.}
\description{
  Variety of functions to support cacheing of calculated canonical
  statistics for Potts Models.  There is some attempt at being 'smart'
  with when to regenerate the statistics.
}
\usage{
generate_t_cache(x, ncolor, t_stat, sizeA, npixel, f, gridcache=NULL)

incvector(v, ncolor)

gengridcache(ncolor, v, ncol)

gensingleton(ncolor)

singleton(x, ncolor, a, idx, gcache=NULL)

gentwopixel(ncolor)

twopixel(x, ncolor, a, idx, gcache=NULL)

twopixel.nonoverlap(x, ncolor, a, idx, gcache=NULL)

genfourpixel(ncolor)

fourpixel(x, ncolor, a, idx, gcache=NULL)

fourpixel.nonoverlap(x, ncolor, a, idx, gcache=NULL)

genthreebythree(ncolor)

ninepixel.nonoverlap(x, ncolor, a, idx, gcache=NULL)

genfourbyfour(ncolor)

sixteenpixel.nonoverlap(x, ncolor, a, idx, gcache=NULL)
}
\arguments{

  \item{t_stat}{numerical vector of length \code{ncolor}.  Contains the
                canonical statistic for the whole image, with one
                component dropped.}

  \item{sizeA}{numerical.  The number of elements in script(A).}

  \item{npixel}{numerical.  The number of pixels in one element of script(A).}

  \item{f}{function.  Takes arguments \code{x}, \code{ncolor}, \code{a},
           \code{idx} and \code{ncolor}.  Returns value of
           \code{t_calc_innergrid} with window A_a replaced by the
           idx'th element of C^A.}

  \item{gridcache}{list.  Optional.  If non-null, it is a list of the
                   elements of C^A.  Indexed by \code{idx}.}

  \item{v}{numerical vector of length \code{|C|^|A|}.  Elements must be
           in 1,\dots,ncolor.}

  \item{x}{2 dimensional matrix, elements in 1, \dots, ncolor,
           representing a Potts model.}

  \item{ncolor}{Number of colors in this Potts Model.}

  \item{ncol}{numerical.  Indicates number of columns to use when
    casting \code{v} as a matrix.}

  \item{a}{numerical.  Indicates which member of script(A) is being
    referenced.}

  \item{idx}{numerical.  Indicates which element of C^A is being referenced.}
  
  \item{gcache}{list.  The list contains matricies.  Elements of the
    matricies must be in 1,\dots,ncolor. Each element of the list is one
    of the elements of C^A.}
	       
}
\details{

  This set of functions is used to generate cached calculations of the
  canonical statistic of a Potts model suitable for passing into
  \code{composite.ll} or \code{gr.composite.ll}.

  \code{generate_t_cache} is the main function, and the others are
  intended to be used in conjunction with it.  See the examples below.

  Of special not are \code{incvector} and \code{gengridcache}.  Using
  these two functions you can easily generate the grid cache for any
  rectangular sized window.  This is how most of the other functions
  work.

  Functions that start with \code{gen} generate a \code{list} of the
  elements of C^A for the window size and number of colors.  Most of
  these functions call \code{gengridcache} with appropriate arguments.

  Functions that have \code{overlap} and \code{nonoverlap} versions
  generate a overlapping and nonoverlapping set of windows respectively.
  
  
  \code{singleton} is for a single pixel window (Besag or MPLE).

  \code{twopixel} does a two horizontal pixel window.

  \code{fourpixel} does a two by two pixel window.

  \code{ninepixel} does a three by three pixel window.

  \code{sixteenpixel} does a four by four pixel window.
  
}
\value{
  \code{incvector} takes a vector and returns a new vector.  These
  vectors can be regarded as elements of C^A cast into the a vector.  If
  you let v_1 be any vector in B = \{ as.vector(c) | c \\in C^A \} then the
  set of vectors \{v_1,\dots,v_\{|C|^|A|\} | v_\{j+1\} = v_j \} = B.

  Functions that start with \code{gen} return a \code{list} of the
  elements of C^A.

  The other functions (e.g. \code{twopixel}, \code{fourpixel}, \dots)
  return the result of replacing the a-th window of x with the idx-th
  element of C^A and calculating \code{calc_t_innergrid} for that window.
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
t_cache_besag <- generate_t_cache(x, ncolor, t_stat, nrow*ncol, 1,
                                  singleton)

}
\keyword{misc}
