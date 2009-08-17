
potts <- function(obj, param, niter,
    boundary = c("torus", "free", "condition"))
UseMethod("potts")

potts.potts <- function(obj, param, niter,
    boundary = c("torus", "free", "condition"))
{
    boundary <- match.arg(boundary)
    if (missing(param)) param <- obj$param
    if (missing(niter)) niter <- obj$niter
    initial <- obj$final
    .Random.seed <- obj$final.seed
    potts.raw(initial, param, niter, boundary)
}

potts.raw <- function(obj, param, niter,
    boundary = c("torus", "free", "condition"))
{
    if (! exists(".Random.seed")) runif(1)
    saveseed <- .Random.seed

    initial.info <- inspectPotts(obj)
    stopifnot(is.numeric(param))
    stopifnot(all(is.finite(param)))
    stopifnot(length(param) == initial.info$ncolor + 1)

    stopifnot(is.numeric(niter))
    niter <- as.integer(niter)
    stopifnot(niter > 0)

    boundary.code <- match(boundary, c("torus", "free", "condition"))

    out.time <- system.time(
    out <- .C("potts", final = obj, param = as.double(param),
        niter = niter, code = as.integer(boundary.code),
        path = matrix(as.integer(0), nrow = length(param), ncol = niter),
        PACKAGE = "potts")
    )
    return(structure(list(initial.seed = saveseed, final.seed = .Random.seed,
        initial = obj, final = out$final, param = param, niter = niter,
        boundary = boundary, path = t(out$path), time = out.time),
        class = "potts"))
}

