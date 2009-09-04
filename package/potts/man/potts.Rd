\name{potts}
\alias{potts}
\alias{potts.raw}
\alias{potts.potts}
\title{Potts Models}
\description{
    Simulate Potts model using Swendsen-Wang algorithm.
}
\usage{
potts(obj, param, nbatch, blen = 1, nspac = 1,
    boundary = c("torus", "free", "condition"))
}
\arguments{
  \item{obj}{an R vector of class \code{"raw"} that encodes a realization
      of a Potts model, typically the output of \code{\link{packPotts}}.
      Alternatively, an object of class \code{"potts"} from a
      previous invocation of this function can be supplied,
      in which case any missing arguments
      are taken from this object.}
  \item{param}{numeric, canonical parameter vector.}
  \item{nbatch}{the number of batches.}
  \item{blen}{the length of batches.}
  \item{nspac}{the spacing of iterations that contribute to batches.}
  \item{boundary}{type of boundary conditions.
      The value of this argument can be abbreviated.}
}
\details{
Runs a Swendsen-Wang algorithm producing a Markov chain with equilibrium
distribution having the specified Potts model.
}
\value{
  an object of class \code{"potts"},
  which is a list containing at least the following components:
  \item{initial}{initial state of Markov chain.}
  \item{final}{final state of Markov chain.}
  \item{initial.seed}{value of \code{.Random.seed} before the run.}
  \item{final.seed}{value of \code{.Random.seed} after the run.}
  \item{time}{running time of Markov chain from \code{\link{system.time}}.}
  \item{param}{canonical parameter vector.}
  \item{nbatch}{the number of batches.}
  \item{blen}{the length of batches.}
  \item{nspac}{the spacing of iterations that contribute to batches.}
  \item{boundary}{the argument \code{boundary}.}
  \item{batch}{an \code{nbatch} by \code{ncolor + 1} matrix, where \code{ncolor}
      is the number of colors of the potts model; each row is the batch means
      for the canonical statistic vector for one batch of Markov chain
      iterations.}
}
\examples{
ncolor <- as.integer(4)
beta <- log(1 + sqrt(ncolor))
theta <- c(rep(0, ncolor), beta)

nrow <- 100
ncol <- 100
x <- matrix(1, nrow = nrow, ncol = ncol)
foo <- packPotts(x, ncolor)

out <- potts(foo, theta, nbatch = 10)
out$path
\dontrun{image(out$final)}
}
\keyword{misc}
