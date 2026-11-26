MkPrime:::.Rhat <- 
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

MkPrime:::.EssBulk <- 
function (x) 
{
    x <- .AsChainMatrix(x)
    splits <- .SplitChains(x)
    .EssMultiChain(.ZScale(splits))
}

MkPrime:::.EssTail <- 
function (x) 
{
    x <- .AsChainMatrix(x)
    splits <- .SplitChains(x)
    q05Ess <- .EssQuantile(splits, 0.05)
    q95Ess <- .EssQuantile(splits, 0.95)
    min(q05Ess, q95Ess)
}

