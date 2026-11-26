# Internal MkPrime functions used by rhat-ess-quartet.R, as installed in
# /nobackup/pjjg18/evotempo-dev/lib (MkPrime 0.0.0.9000, built 2026-10-16 12:29 UTC;
# source commit not recorded). Written by deparse(); for the record only.

`.AsChainMatrix` <- 
function (x) 
{
    if (is.null(dim(x))) {
        matrix(x, ncol = 1L)
    }
    else {
        as.matrix(x)
    }
}

`.Autocovariance` <- 
function (x) 
{
    n <- length(x)
    vx <- var(x)
    if (vx == 0) 
        return(rep(0, n))
    m <- nextn(n)
    m2 <- 2L * m
    yc <- x - mean(x)
    yc <- c(yc, rep.int(0, m2 - n))
    ac <- Re(fft(Mod(fft(yc))^2, inverse = TRUE))[seq_len(n)]
    ac <- ac/ac[1] * vx * (n - 1)/n
    ac
}

`.EssBulk` <- 
function (x) 
{
    x <- .AsChainMatrix(x)
    splits <- .SplitChains(x)
    .EssMultiChain(.ZScale(splits))
}

`.EssMultiChain` <- 
function (x) 
{
    x <- as.matrix(x)
    nChains <- ncol(x)
    nIter <- nrow(x)
    if (nIter < 3L || .IsConstant(x)) 
        return(NA_real_)
    acovList <- apply(x, 2, .Autocovariance, simplify = FALSE)
    acovMeans <- Reduce(`+`, acovList)/nChains
    meanVar <- acovMeans[1] * nIter/(nIter - 1)
    varPlus <- meanVar * (nIter - 1)/nIter
    if (nChains > 1L) {
        varPlus <- varPlus + var(colMeans(x))
    }
    if (varPlus <= 0 || !is.finite(varPlus)) 
        return(NA_real_)
    rhoHat <- numeric(nIter)
    t <- 0L
    rhoEven <- 1
    rhoOdd <- 1 - (meanVar - acovMeans[2])/varPlus
    rhoHat[1] <- rhoEven
    rhoHat[2] <- rhoOdd
    while (t < nIter - 5L && !is.nan(rhoEven + rhoOdd) && (rhoEven + 
        rhoOdd) > 0) {
        t <- t + 2L
        rhoEven <- 1 - (meanVar - acovMeans[t + 1L])/varPlus
        rhoOdd <- 1 - (meanVar - acovMeans[t + 2L])/varPlus
        if ((rhoEven + rhoOdd) >= 0) {
            rhoHat[t + 1L] <- rhoEven
            rhoHat[t + 2L] <- rhoOdd
        }
    }
    maxT <- t
    if (rhoEven > 0) 
        rhoHat[maxT + 1L] <- rhoEven
    t <- 0L
    while (t <= maxT - 4L) {
        t <- t + 2L
        if (rhoHat[t + 1L] + rhoHat[t + 2L] > rhoHat[t - 1L] + 
            rhoHat[t]) {
            rhoHat[t + 1L] <- (rhoHat[t - 1L] + rhoHat[t])/2
            rhoHat[t + 2L] <- rhoHat[t + 1L]
        }
    }
    tauHat <- -1 + 2 * sum(rhoHat[seq_len(maxT)]) + rhoHat[maxT + 
        1L]
    ess <- nChains * nIter/tauHat
    totalN <- nChains * nIter
    tauBound <- 1/log10(totalN)
    if (tauHat < tauBound) {
        ess <- totalN * log10(totalN)
    }
    max(ess, 1)
}

`.EssQuantile` <- 
function (x, prob) 
{
    if (.IsConstant(x)) 
        return(NA_real_)
    indicator <- x <= quantile(x, prob)
    dim(indicator) <- dim(x)
    .EssMultiChain(indicator)
}

`.EssTail` <- 
function (x) 
{
    x <- .AsChainMatrix(x)
    splits <- .SplitChains(x)
    q05Ess <- .EssQuantile(splits, 0.05)
    q95Ess <- .EssQuantile(splits, 0.95)
    min(q05Ess, q95Ess)
}

`.FoldDraws` <- 
function (x) 
{
    abs(x - median(x))
}

`.IsConstant` <- 
function (x, tol = .Machine$double.eps) 
{
    if (anyNA(x) || any(is.infinite(x))) 
        return(TRUE)
    rng <- range(x)
    abs(rng[2] - rng[1]) < tol
}

`.Rhat` <- 
function (x) 
{
    x <- .AsChainMatrix(x)
    splits <- .SplitChains(x)
    if (nrow(splits) < 4L) 
        return(NA_real_)
    bulkRhat <- .RhatClassical(.ZScale(splits))
    tailRhat <- .RhatClassical(.ZScale(.FoldDraws(splits)))
    max(bulkRhat, tailRhat)
}

`.RhatClassical` <- 
function (x) 
{
    if (.IsConstant(x)) 
        return(NA_real_)
    nIter <- nrow(x)
    chainMeans <- colMeans(x)
    chainVars <- apply(x, 2, var)
    varBetween <- nIter * var(chainMeans)
    varWithin <- mean(chainVars)
    if (varWithin == 0) 
        return(Inf)
    sqrt((varBetween/varWithin + nIter - 1)/nIter)
}

`.SplitChains` <- 
function (x) 
{
    x <- as.matrix(x)
    nIter <- nrow(x)
    half <- nIter%/%2L
    if (half < 1L) 
        return(x)
    offset <- if (nIter%%2L == 1L) 
        1L
    else 0L
    cbind(x[seq_len(half), , drop = FALSE], x[half + offset + 
        seq_len(half), , drop = FALSE])
}

`.ZScale` <- 
function (x) 
{
    r <- rank(x, ties.method = "average")
    n <- length(x)
    z <- qnorm((r - 3/8)/(n + 1/4))
    z[is.na(x)] <- NA
    if (!is.null(dim(x))) {
        dim(z) <- dim(x)
    }
    z
}
