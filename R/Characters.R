#' Mk likelihood of discrete characters
#'
#' Log-likelihood of discrete characters under the Mk (Jukes-Cantor) model
#' with `nState` states and a uniform root, by pruning on given branch lengths
#' in expected changes per character.
#' The C++ engine uses the MkPrime pruning kernels; used by [TempoMCMC()].
#'
#' Sampled ancestors are tips on zero-length branches, so their states are
#' observed at the node that bears them.
#'
#' With `characters` from [MkCharacters()], each character has its own number
#' of states and may be ordered or hold sets of states, the likelihood may be
#' conditioned on variability, and rates may vary among characters (discrete
#' Gamma with shape `alpha`); see [MkCharacters()].
#'
#' @param tree A [`TempoTree`][TempoTree()], or a tree that it accepts.
#' @param characters Integer matrix with one row per tip (in the order of the
#' tree's tips) and one column per character, giving states `0` to
#' `nState - 1`, or `NA` where unknown; or characters prepared by
#' [MkCharacters()].
#' @param lengths Numeric vector giving the length of the branch above each
#' node (tips, then internal nodes), in expected changes; ignored at the root.
#' @param nState Integer giving the number of states; ignored for
#' [MkCharacters()] input, which records the number of states per character.
#' @param alpha Numeric giving the shape of the Gamma distribution of rates
#' among characters, for [MkCharacters()] input with `nRateCat > 1`;
#' otherwise ignored.
#' @inheritParams TempoEngine
#' @returns `MkLikelihood()` returns the natural logarithm of the likelihood.
#' @examples
#' tree <- structure(list(edge = rbind(c(4, 5), c(5, 1), c(5, 2), c(4, 3)),
#'                        edge.length = c(1, 2, 2, 3), tip.label = c("a", "b", "c"),
#'                        Nnode = 2), class = "phylo")
#' chars <- matrix(c(0, 0, 1, 1, NA, 1), 3)
#' MkLikelihood(tree, chars, c(0.2, 0.2, 0.3, 0, 0.1))
#'
#' # Three states, one character ordered, a polymorphic cell, conditioned on
#' # variability, with Gamma rate variation
#' cells <- matrix(c("0", "1", "2", "{01}", "1", "0"), 3)
#' mk <- MkCharacters(cells, nState = 3, ordered = c(TRUE, FALSE),
#'                    coding = "variable", nRateCat = 4)
#' MkLikelihood(tree, mk, c(0.2, 0.2, 0.3, 0, 0.1), alpha = 0.5)
#' @family dCET inference
#' @export
MkLikelihood <- function(tree, characters, lengths, nState = 2L,
                         engine = getOption("EvoTempo.engine", "cpp"),
                         alpha = NULL) {
  tree <- TempoTree(tree)
  if (inherits(characters, "MkCharacters")) {
    return(.MkCharactersLogLik(tree, characters, lengths, alpha, engine))
  }
  nTip <- tree[["nTip"]]
  characters <- as.matrix(characters)
  if (nrow(characters) != nTip) {
    stop("`characters` needs one row per tip")
  }
  if (TempoEngine(engine) == "cpp") {
    states <- matrix(as.integer(characters), nTip)
    return(.MkLikelihoodCpp(tree, states, as.double(lengths), as.integer(nState)))
  }
  nChar <- ncol(characters)
  partial <- vector("list", nTip + tree[["nNode"]])
  for (i in seq_len(nTip)) {
    p <- matrix(0, nState, nChar)
    obs <- characters[i, ]
    p[, is.na(obs)] <- 1
    known <- which(!is.na(obs))
    p[cbind(obs[known] + 1L, known)] <- 1
    partial[[i]] <- p
  }
  logScale <- 0
  for (v in tree[["postorder"]]) {
    value <- 1
    for (kid in tree[["kids"]][v - nTip, ]) {
      ell <- lengths[[kid]]
      stay <- 1 / nState + (nState - 1) / nState * exp(-nState * ell / (nState - 1))
      move <- (1 - stay) / (nState - 1)
      L <- partial[[kid]]
      value <- value * (rep(move * colSums(L), each = nState) + (stay - move) * L)
    }
    scale <- value[cbind(max.col(t(value), "first"), seq_len(nChar))]
    partial[[v]] <- value / rep(scale, each = nState)
    logScale <- logScale + sum(log(scale))
  }
  sum(log(colSums(partial[[tree[["root"]]]]) / nState)) + logScale
}

#' Characters for the Mk model
#'
#' Prepares a matrix of discrete characters for [MkLikelihood()] and
#' [TempoMCMC()], with a number of states per character, ordered characters,
#' ambiguous cells, conditioning on variability and among-character rate
#' variation.
#'
#' Unordered characters follow the Mk (Jukes-Cantor) model with \eqn{K}
#' states (RevBayes `fnJC(K)`).
#' Ordered characters change only between adjacent states, at equal rates
#' (RevBayes `fnOrderedRateMatrix(maxState = K)`).
#' Both rate matrices are scaled to one expected change per unit of branch
#' length, and both give a uniform distribution at the root.
#' A two-state ordered character is unordered.
#' Unordered transition probabilities are computed in closed form; ordered ones
#' by uniformization with repeated squaring, which sums only non-negative terms
#' and so keeps small probabilities accurate, with \eqn{P(0) = I} exactly.
#'
#' A cell holds one state, a set of states, whose likelihood is the sum over
#' its states, or is unknown.
#' A set that covers all \eqn{K} states is unknown.
#'
#' With `coding = "variable"` the likelihood is conditioned on each
#' character being variable among the tips where it is scored (Lewis 2001),
#' as RevBayes computes `dnPhyloCTMC(coding = "variable")`, with unknown cells
#' as gaps.
#' Characters with the same number of states, ordering and missing-data
#' mask (the set of tips with a scored cell) share one correction,
#' \eqn{-n_m \log(1 - P_m)}, where \eqn{P_m} is the probability that every
#' scored tip of mask \eqn{m} has the same state, averaged over rate
#' categories, and \eqn{n_m} the number of characters with that mask.
#' The conditioning is exact only for characters without partial sets: as in
#' RevBayes, a character such as `0, 0, {01}` is kept although some of the
#' histories it sums over are constant.
#' A character whose scored cells all hold the same state or set cannot be
#' variable, and is an error (RevBayes excludes it).
#'
#' With `nRateCat` \eqn{> 1}, the branch lengths of each character are
#' multiplied by a rate \eqn{r_c}, and its likelihood is the mean over the
#' `nRateCat` equiprobable categories \eqn{c}, whose rates are those of
#' [GammaRates()] (a Gamma distribution with shape and rate \eqn{\alpha}).
#' The shape is the argument `alpha` of [MkLikelihood()], and a parameter of
#' [TempoMCMC()] (`start$alpha`, with prior `charPriors$alpha`).
#'
#' The row names of `characters`, if any (other than `1, 2, ...`), are kept
#' and must equal the tip labels of the tree that the characters are scored
#' on.
#'
#' A character's likelihood is computed without rescaling in the C++ engine,
#' so a single character with likelihood below about \eqn{10^{-308}} (log
#' likelihood below about -709; reachable with many tips at very low rates)
#' gives `-Inf` there, but not in the R engine.
#'
#' @param characters Matrix with one row per tip (in the order of the tree's
#' tips) and one column per character: either integer states `0` to
#' \eqn{K - 1}, with `NA` where unknown, or character cells holding a state
#' (`"0"`, `"12"`), a set of states (such as `"{01}"`, `"(012)"` or `"0/12"`),
#' or `"?"`, `"-"`, `""` or `NA` where unknown.
#' @param nState Integer giving the number of states \eqn{K} of each
#' character: one value, or one per character.
#' By default, one more than the largest state of the character, and at
#' least 2; RevBayes instead takes \eqn{K} from the matrix's symbols, so give
#' `nState` to match it (a message says so when ordered characters take the
#' default).
#' @param ordered Logical, one value or one per character: whether the
#' character is ordered.
#' @param coding Character string: `"all"` for no conditioning, or
#' `"variable"` to condition on variability.
#' @param nRateCat Integer giving the number of Gamma rate categories;
#' 1 for no rate variation among characters.
#' @returns `MkCharacters()` returns an object of class `MkCharacters`: a
#' list with elements `nTip`, `tipLabel` (the row names, or `NULL`), `nChar`,
#' `nState` and `ordered` (per character), `coding`, `nRateCat` and
#' `partitions`, one per combination of ordering and number of states, each
#' holding the tip partial likelihoods of its characters, followed by those of
#' the constant patterns of its missing-data masks if `coding = "variable"`.
#' @references
#' Lewis PO (2001). A likelihood approach to estimating phylogeny from
#' discrete morphological character data. Systematic Biology 50:913--925.
#' \doi{10.1080/106351501753462876}
#'
#' Yang Z (1994). Maximum likelihood phylogenetic estimation from DNA
#' sequences with variable rates over sites: approximate methods. Journal of
#' Molecular Evolution 39:306--314. \doi{10.1007/BF00160154}
#' @examples
#' cells <- matrix(c("0", "1", "2", "{01}", "1", "0", "0", "1", "?"), 3)
#' mk <- MkCharacters(cells, nState = c(3, 2, 2), ordered = c(TRUE, FALSE, FALSE),
#'                    coding = "variable", nRateCat = 4)
#' mk$nState
#' @family dCET inference
#' @export
MkCharacters <- function(characters, nState = NULL, ordered = FALSE,
                         coding = c("all", "variable"), nRateCat = 1L) {
  coding <- match.arg(coding)
  sets <- .CellSets(characters)
  nTip <- nrow(sets)
  nChar <- ncol(sets)
  if (nTip < 2L || nChar < 1L) {
    stop("`characters` needs at least two rows (tips) and one column")
  }
  tipLabel <- rownames(characters)
  if (identical(tipLabel, as.character(seq_len(nTip)))) {
    tipLabel <- NULL
  }
  maxState <- apply(sets, 2, function(col) {
    states <- unlist(col)
    if (length(states)) max(states) else 0L
  })
  if (!is.logical(ordered) || !length(ordered) %in% c(1L, nChar) || anyNA(ordered)) {
    stop("`ordered` must be TRUE or FALSE, or one value per character")
  }
  if (is.null(nState)) {
    nState <- pmax(2L, maxState + 1L)
    if (any(rep_len(ordered, nChar) & nState > 2L)) {
      message("MkCharacters(): `nState` taken from the largest observed state of ",
              "each character; give `nState` to match a matrix whose symbols ",
              "allow more states (RevBayes takes K from the symbols)")
    }
  }
  if (!is.numeric(nState) || !length(nState) %in% c(1L, nChar) ||
      any(!is.finite(nState) | nState != round(nState) | nState < 2)) {
    stop("`nState` must give one whole number >= 2, or one per character")
  }
  nState <- rep_len(as.integer(nState), nChar)
  if (any(maxState >= nState)) {
    stop("States exceed `nState` - 1 in character(s) ",
         paste(which(maxState >= nState), collapse = ", "))
  }
  ordered <- rep_len(ordered, nChar) & nState > 2L
  if (length(nRateCat) != 1L || !is.numeric(nRateCat) || !is.finite(nRateCat) ||
      nRateCat < 1 || nRateCat != round(nRateCat)) {
    stop("`nRateCat` must be a whole number >= 1")
  }
  # Partial likelihoods [state, character, tip]; a full set is unknown
  partial <- array(1, c(max(nState), nChar, nTip))
  scored <- matrix(FALSE, nTip, nChar)
  for (j in seq_len(nChar)) {
    for (i in seq_len(nTip)) {
      s <- sets[[i, j]]
      if (length(s) && length(s) < nState[[j]]) {
        partial[, j, i] <- 0
        partial[s + 1L, j, i] <- 1
        scored[i, j] <- TRUE
      }
    }
  }
  if (coding == "variable") {
    constant <- vapply(seq_len(nChar), function(j) {
      key <- vapply(sets[scored[, j], j], paste, "", collapse = ",")
      length(unique(key)) < 2L
    }, TRUE)
    if (any(constant)) {
      stop("Character(s) ", paste(which(constant), collapse = ", "),
           " cannot be variable (every scored cell holds the same state or ",
           "set); remove them for `coding = \"variable\"`")
    }
  }
  groups <- unique(data.frame(ordered = ordered, nState = nState))
  groups <- groups[order(groups[["ordered"]], groups[["nState"]]), , drop = FALSE]
  partitions <- lapply(seq_len(nrow(groups)), function(g) {
    K <- groups[["nState"]][[g]]
    ord <- groups[["ordered"]][[g]]
    cols <- which(nState == K & ordered == ord)
    tip <- partial[seq_len(K), cols, , drop = FALSE]
    count <- integer(0)
    nConst <- 0L
    if (coding == "variable") {
      key <- apply(scored[, cols, drop = FALSE], 2, function(o) {
        paste(which(o), collapse = ",")
      })
      first <- !duplicated(key)
      mask <- scored[, cols[first], drop = FALSE]
      count <- tabulate(match(key, key[first]), sum(first))
      # Under Jukes-Cantor every constant pattern has the same probability,
      # so one pattern (state 0) per mask stands for all K
      nConst <- if (ord) K else 1L
      const <- array(1, c(K, ncol(mask) * nConst, nTip))
      for (m in seq_len(ncol(mask))) {
        for (s in seq_len(nConst)) {
          const[, (m - 1L) * nConst + s, mask[, m]] <- as.numeric(seq_len(K) == s)
        }
      }
      combined <- array(0, c(K, length(cols) + ncol(mask) * nConst, nTip))
      combined[, seq_along(cols), ] <- tip
      combined[, length(cols) + seq_len(ncol(mask) * nConst), ] <- const
      tip <- combined
    }
    list(nState = K, ordered = ord, columns = cols, nData = length(cols),
         nConst = nConst, count = count, nCol = dim(tip)[[2]],
         tip = as.vector(tip), Q = .MkRateMatrix(K, ord))
  })
  structure(list(nTip = nTip, tipLabel = tipLabel, nChar = nChar, nState = nState,
                 ordered = ordered, coding = coding, nRateCat = as.integer(nRateCat),
                 partitions = partitions),
            class = "MkCharacters")
}

#' Discrete Gamma rates
#'
#' Rate multipliers of `nCat` equiprobable categories of a Gamma distribution
#' with shape and rate `alpha` (mean 1), each the mean rate within its
#' category (Yang 1994): the discretisation that RevBayes
#' `fnDiscretizeGamma(alpha, alpha, nCat)` approximates.
#' The cut points here come from R's `qgamma()`; RevBayes's chi-square
#' quantile (AS 91, relative tolerance 5e-7) moves its rates by up to a few
#' times \eqn{10^{-8}} at some shapes (e.g. 0.8 and 1.1), which changes
#' log-likelihoods by up to about \eqn{10^{-6}}.
#'
#' @param alpha Numeric giving the shape (and rate) of the Gamma
#' distribution.
#' @param nCat Integer giving the number of categories.
#' @returns `GammaRates()` returns a numeric vector of `nCat` rates with
#' mean 1.
#' @references
#' Yang Z (1994). Maximum likelihood phylogenetic estimation from DNA
#' sequences with variable rates over sites: approximate methods. Journal of
#' Molecular Evolution 39:306--314. \doi{10.1007/BF00160154}
#' @examples
#' GammaRates(0.5)
#' @family dCET inference
#' @export
GammaRates <- function(alpha, nCat = 4L) {
  if (length(alpha) != 1L || !is.numeric(alpha) || !is.finite(alpha) || alpha <= 0) {
    stop("`alpha` must be a positive number")
  }
  if (length(nCat) != 1L || !is.numeric(nCat) || !is.finite(nCat) || nCat < 1 ||
      nCat != round(nCat)) {
    stop("`nCat` must be a whole number >= 1")
  }
  if (nCat == 1L) {
    return(1)
  }
  cut <- stats::qgamma(seq_len(nCat - 1L) / nCat, shape = alpha, rate = alpha)
  diff(c(0, stats::pgamma(cut * alpha, alpha + 1), 1)) * nCat
}

# Character matrix -> matrix (list) of integer state sets, empty if unknown
.CellSets <- function(characters) {
  characters <- as.matrix(characters)
  if (is.numeric(characters) || is.logical(characters)) {
    if (any(!is.na(characters) &
            (characters < 0 | characters != round(characters)))) {
      stop("Numeric `characters` must be whole numbers >= 0, or NA")
    }
    cells <- lapply(as.vector(characters), function(x) {
      if (is.na(x)) integer(0) else as.integer(x)
    })
  } else {
    cells <- lapply(as.vector(characters), function(x) {
      if (is.na(x)) {
        return(integer(0))
      }
      x <- gsub("\\s", "", x)
      if (x %in% c("", "?", "-")) {
        return(integer(0))
      }
      if (grepl("^[0-9]+$", x)) {
        return(as.integer(x))  # one state, possibly of several digits
      }
      if (grepl("^[0-9]+(/[0-9]+)+$", x)) {
        return(sort(unique(as.integer(strsplit(x, "/")[[1]]))))  # "0/12"
      }
      if (!grepl("^[\\{\\(][0-9]+[\\}\\)]$", x)) {
        stop("Cannot read character cell \"", x, "\"")
      }
      sort(unique(as.integer(strsplit(gsub("[^0-9]", "", x), "")[[1]])))  # "{01}"
    })
  }
  matrix(cells, nrow(characters), ncol(characters))
}

# Rate matrix of K states, mean rate 1 under the uniform distribution.
# Ordered: adjacent states only, at rate r, mean rate 2 (K - 1) r / K = 1
.MkRateMatrix <- function(K, ordered) {
  Q <- matrix(0, K, K)
  if (ordered) {
    r <- K / (2 * (K - 1))
    Q[cbind(seq_len(K - 1), 2:K)] <- r
    Q[cbind(2:K, seq_len(K - 1))] <- r
  } else {
    Q[] <- 1 / (K - 1)
    diag(Q) <- 0
  }
  diag(Q) <- -rowSums(Q)
  Q
}

# Transition probabilities exp(Q t) by uniformization and squaring (#98,
# R98-13): with mu = max_i -Q_ii and B = I + Q / mu >= 0, the argument
# x = mu t is halved s times until y = x / 2^s <= 1/2, exp(Q t / 2^s) =
# exp(-y) sum_{n <= K + 17} y^n / n! B^n (truncation below 1e-19 relative to
# any term that the sum contains), and the result is squared s times.  Every
# step adds and multiplies non-negative numbers, so small probabilities keep
# their relative accuracy, and P(0) = I exactly.  Mirrored step by step by
# UniformizedTransition() in src/characters.cpp.
.MkTransition <- function(Q, t) {
  K <- nrow(Q)
  if (!(t > 0)) {
    return(diag(K))
  }
  mu <- max(-diag(Q))
  x <- mu * t
  s <- if (x > 0.5) as.integer(ceiling(log2(x / 0.5))) else 0L
  y <- x / 2^s
  B <- diag(K) + Q / mu
  term <- diag(K)
  total <- diag(K)
  for (n in seq_len(K + 17L)) {
    term <- (term %*% B) * (y / n)
    total <- total + term
  }
  P <- exp(-y) * total
  for (i in seq_len(s)) {
    P <- P %*% P
  }
  P
}

# log(1 - exp(lp)) for lp < 0, without cancellation (Maechler 2012)
.Log1mExp <- function(lp) {
  ifelse(lp > -log(2), log(-expm1(lp)), log1p(-exp(lp)))
}

# Log-likelihood of MkCharacters on branch lengths, with rate categories at
# shape `alpha` (or the given `rates`); the Mkv correction per mask if
# coding = "variable".  -Inf if a constant pattern's probability reaches 1
.MkCharactersLogLik <- function(tree, characters, lengths, alpha,
                                engine = getOption("EvoTempo.engine", "cpp"),
                                rates = NULL) {
  if (characters[["nTip"]] != tree[["nTip"]]) {
    stop("`characters` needs one row per tip")
  }
  labels <- characters[["tipLabel"]]
  if (!is.null(labels) && !is.null(tree[["tipLabel"]]) &&
      !identical(labels, tree[["tipLabel"]])) {
    stop("The row names of `characters` do not match the tree's tip labels ",
         "(in the tree's order)")
  }
  if (is.null(rates)) {
    rates <- if (characters[["nRateCat"]] > 1L) {
      if (is.null(alpha)) {
        stop("`alpha` is needed with `nRateCat` > 1")
      }
      GammaRates(alpha, characters[["nRateCat"]])
    } else {
      1
    }
  }
  variable <- characters[["coding"]] == "variable"
  if (TempoEngine(engine) == "cpp") {
    return(.MkCharactersCpp(tree, characters[["partitions"]], as.double(lengths),
                            as.double(rates), variable))
  }
  total <- 0
  for (g in characters[["partitions"]]) {
    # Site likelihoods (log), characters x categories
    logLik <- vapply(rates, function(r) .MkColumnLogLik(tree, g, lengths * r),
                     numeric(g[["nCol"]]))
    logLik <- matrix(logLik, g[["nCol"]])
    top <- apply(logLik, 1, max)
    siteLik <- log(rowMeans(exp(logLik - top))) + top
    # Likelihood 0 in every category (as the C++ engine's sum): -Inf, not NaN
    siteLik[top == -Inf] <- -Inf
    data <- seq_len(g[["nData"]])
    total <- total + sum(siteLik[data])
    if (variable) {
      # log P_m: of state 0 times K (unordered), or of every state (ordered)
      logConst <- siteLik[-data]
      logMask <- if (g[["ordered"]]) {
        m <- matrix(logConst, g[["nConst"]])
        mTop <- apply(m, 2, max)
        log(colSums(exp(m - rep(mTop, each = nrow(m))))) + mTop
      } else {
        log(g[["nState"]]) + logConst
      }
      if (any(!(logMask < 0))) {
        return(-Inf)
      }
      total <- total - sum(g[["count"]] * .Log1mExp(logMask))
    }
  }
  total
}

# Pure-R pruning of one partition's columns, with per-node rescaling:
# log site likelihoods with a uniform root.  Unordered: the closed-form
# Jukes-Cantor probabilities of MkPrime's JCRates; ordered: .MkTransition()
.MkColumnLogLik <- function(tree, g, lengths) {
  nTip <- tree[["nTip"]]
  K <- g[["nState"]]
  tip <- array(g[["tip"]], c(K, g[["nCol"]], nTip))
  part <- vector("list", nTip + tree[["nNode"]])
  for (i in seq_len(nTip)) {
    part[[i]] <- matrix(tip[, , i], K)
  }
  logScale <- numeric(g[["nCol"]])
  for (v in tree[["postorder"]]) {
    a <- 1
    for (kid in tree[["kids"]][v - nTip, ]) {
      ell <- lengths[[kid]]
      L <- part[[kid]]
      a <- a * if (g[["ordered"]]) {
        .MkTransition(g[["Q"]], ell) %*% L
      } else {
        negExpm1 <- -expm1(-K * ell / (K - 1))
        pSame <- 1 / K + (1 - 1 / K) * (1 - negExpm1)
        pDiff <- negExpm1 / K
        rep(pDiff * colSums(L), each = K) + (pSame - pDiff) * L
      }
    }
    mx <- apply(a, 2, max)
    # A column with likelihood 0 (e.g. changes on branches of length ~0, at a
    # Gamma category rate near 0) stays 0, with log scale -Inf, not NaN
    part[[v]] <- a / rep(ifelse(mx > 0, mx, 1), each = K)
    logScale <- logScale + log(mx)
  }
  log(colSums(part[[tree[["root"]]]]) / K) + logScale
}

# Number of characters in a matrix or MkCharacters
.NCharacters <- function(characters) {
  if (inherits(characters, "MkCharacters")) {
    characters[["nChar"]]
  } else {
    NCOL(characters)
  }
}

# Tempo-weighted branch lengths c * sum_k exp(gamma x_k) dwell_k, per node
# (with an origin, the root's is the stem's, which MkLikelihood() ignores)
.CharacterLengths <- function(tree, x, paths, rate, gamma) {
  if (.UseCpp()) {
    return(.CharacterLengthsCpp(tree, x, paths, rate, gamma))
  }
  age <- tree[["age"]]
  out <- numeric(tree[["nTip"]] + tree[["nNode"]])
  for (v in .PathNodes(tree)) {
    int <- .Intervals(paths[[v]], .TopAge(tree, v), age[[v]])
    out[v] <- rate * sum(exp(gamma * x[int[["bin"]]]) * (int[["hi"]] - int[["lo"]]))
  }
  out
}
