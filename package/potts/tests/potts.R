
library(potts)

set.seed(42)

ncolor <- as.integer(4)
beta <- log(1 + sqrt(ncolor))
theta <- c(rep(0, ncolor), beta)

nrow <- 100
ncol <- 100
x <- matrix(1, nrow = nrow, ncol = ncol)
foo <- packPotts(x, ncolor)

out <- potts(foo, theta, niter = 10)
out$path

out <- potts(foo, theta, niter = 10, boundary = "free")
out$path

out <- potts(foo, theta, niter = 10, boundary = "condition")
out$path

