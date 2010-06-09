\name{composite.ll}
\alias{composite.ll}
\alias{gr.composite.ll}
\title{Composite Log Likelihood for Potts Models}
\description{
  Calculate Composite Log Likelihood (CLL) and the gradient of the CLL
  for Potts models.
}
\usage{
composite.ll(theta, t_stat, t_cache=NULL, fapply=lapply)
gr.composite.ll(theta, t_stat, t_cache=NULL, fapply=lapply) 
}
\arguments{

  \item{theta}{numeric canonical parameter vector.  The CLL will be
               evaluated at this point.  It is assumed that a component
               corresponding to one of the colors has been dropped.}

  \item{t_stat}{numeric, canonical statistic vector.  The value of the
                canonical statistic for the full image.  It is assumed
                that a component corresponding to one of the colors has
                been dropped.  It must be the same component as for
                \code{theta}.}

  \item{t_cache}{numeric 3 dimensional array.  Theoretically optional,
                 but required for now.  t_cache[i,j,] = t(c_j union
                 x_\{L \\ A\}), with the same component dropped as
                 \code{t_stat} and \code{theta}.  i means the ith
                 enumerated element of script(A), denoted A_i.  j
                 corresponds to the jth element of C^\{A_i\}}
	       
  \item{fapply}{function.  Expected to function as lapply does.  Usefull
    for enabling parallel processing.  E.g. use the mclapply function
    from the multicore package.}

}
\details{
  
  For the given value of \code{theta} \code{composite.ll} and
  \code{gr.composite.ll} calculate the CLL and the gradient of the CLL
  respectively for a realized Potts model represented by \code{t_stat}
  and \code{t_cache}.
  
}
\value{
  \code{composite.ll} returns CLL evaluated at \code{theta}.

  \code{gr.composite.ll} returns a numeric vector of length
  \code{length(theta)} containing the gradient of the CLL at \code{theta}.
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
t_cache_mple <- generate_t_cache(x, ncolor, t_stat, nrow*ncol, 1,
                                 singleton)
t_cache_two <- generate_t_cache(x, ncolor, t_stat, nrow*ncol/2, 2,
                                twopixel.nonoverlap)

theta.initial <- 1:ncolor

composite.ll(theta[-1], t_stat, t_cache_mple)
gr.composite.ll(theta[-1], t_stat, t_cache_mple)

\dontrun{
optim.mple <- optim(theta.initial, composite.ll, gr=gr.composite.ll,
                    t_stat, t_cache_mple, method="BFGS",
                    control=list(fnscale=-1))
optim.mple$par

optim.two <- optim(theta.initial, composite.ll, gr=gr.composite.ll,
                   t_stat, t_cache_two, method="BFGS",
                   control=list(fnscale=-1))
optim.two$par
}

\dontrun{
# or use mclapply to speed things up.
library(multicore)
optim.two <- optim(theta.initial, composite.ll, gr=gr.composite.ll,
                   t_stat, t_cache_two, mclapply, method="BFGS",
                   control=list(fnscale=-1))
optim.two$par

}

}
\keyword{misc}
