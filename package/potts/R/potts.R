
potts <- function(obj, param, nbatch, blen = 1, nspac = 1,
    boundary = c("torus", "free", "condition"))
UseMethod("potts")

potts.potts <- function(obj, param, nbatch, blen = 1, nspac = 1,
    boundary = c("torus", "free", "condition"))
{
    boundary <- match.arg(boundary)
    if (missing(param)) param <- obj$param
    if (missing(nbatch)) nbatch <- obj$nbatch
    if (missing(blen)) blen <- obj$blen
    if (missing(nspac)) nspac <- obj$nspac
    initial <- obj$final
    .Random.seed <- obj$final.seed
    potts.raw(initial, param, nbatch, blen, nspac, boundary)
}

potts.raw <- function(obj, param, nbatch, blen = 1, nspac = 1,
    boundary = c("torus", "free", "condition"))
{
    boundary <- match.arg(boundary)

    if (! exists(".Random.seed")) runif(1)
    saveseed <- .Random.seed

    initial.info <- inspectPotts(obj)
    stopifnot(is.numeric(param))
    stopifnot(all(is.finite(param)))
    if(length(param) != initial.info$ncolor + 1)
        stop("length(param) not number of colors + 1")

    stopifnot(is.numeric(nbatch))
    stopifnot(length(nbatch) == 1)
    stopifnot(nbatch == as.integer(nbatch))
    nbatch <- as.integer(nbatch)
    stopifnot(nbatch > 0)

    stopifnot(is.numeric(blen))
    stopifnot(length(blen) == 1)
    stopifnot(blen == as.integer(blen))
    blen <- as.integer(blen)
    stopifnot(blen > 0)

    stopifnot(is.numeric(nspac))
    stopifnot(length(nspac) == 1)
    stopifnot(nspac == as.integer(nspac))
    nspac <- as.integer(nspac)
    stopifnot(nspac > 0)

    boundary.code <- match(boundary, c("torus", "free", "condition"))

    out.time <- system.time(
    out <- .C("potts", final = obj, param = as.double(param),
        nbatch = nbatch, blen = blen, nspac = nspac,
        code = as.integer(boundary.code),
        batch = matrix(as.double(0), nrow = length(param), ncol = nbatch),
        PACKAGE = "potts")
    )
    return(structure(list(initial.seed = saveseed, final.seed = .Random.seed,
        initial = obj, final = out$final, param = param, nbatch = nbatch,
        blen = blen, nspac = nspac,
        boundary = boundary, batch = t(out$batch), time = out.time),
        class = "potts"))
}

