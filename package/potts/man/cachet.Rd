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
  know the value of \code{t( y union X\A_1)} this would be
  \code{calc_t_full(X, ncolor) + calc_t_innergrid(X, ncolor, y, i, j) -
  calc_t_innergrid(X, ncolor, A_1, i, j)}
}
\examples{
}
\keyword{misc}
