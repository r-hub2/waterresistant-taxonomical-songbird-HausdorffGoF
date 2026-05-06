# =============================================================================
# Private helpers
# =============================================================================

.validate_psi <- function(psi) {
  if (!is.numeric(psi) || length(psi) != 2)
    stop("`scale_psi` must be a numeric vector of length 2.")
  psi <- sort(psi)
  if (psi[1] <= 0 || psi[2] >= 1)
    stop("Both elements of `scale_psi` must be strictly between 0 and 1.")
  if (psi[1] == psi[2])
    stop("The two elements of `scale_psi` must be distinct.")
  psi  # returns c(psi_low, psi_high)
}

# Numerically invert a CDF at a single probability p: solve CDF(u) = p.
#   pdf supplied:  Newton-Raphson (derivative = pdf), mirroring H_stat_1s_1d.
#   pdf absent:    uniroot with an automatically expanded bracket.
.invert_cdf <- function(CDF, pdf = NULL, p, tol = 1e-10, max.init = 1000) {
  if (!is.null(pdf)) {
    g  <- function(u) CDF(u) - p   # g'(u) = pdf(u)
    u  <- 0
    h  <- g(u) / pdf(u)            # initial step before entering the loop
    num.init <- 0
    while (abs(h) > tol) {
      h  <- g(u) / pdf(u)
      u  <- u - h
      num.init <- num.init + 1
      if (max.init < num.init)
        stop("Maximum number of iterations reached in .invert_cdf.")
    }
    return(u)
  }
  uniroot(function(u) CDF(u) - p, interval = c(0, 1),
          extendInt = "yes", tol = tol)$root
}

# Compute sigma* for the one-sample case (Proposition 13, DJK 2026):
#   sigma* = (psi[2] - psi[1]) / (F^{-1}(psi[2]) - F^{-1}(psi[1]))
# Priority: CDFinverse > Newton-Raphson (pdf supplied) > uniroot (CDF only).
# max.init is forwarded to .invert_cdf and used only on the Newton-Raphson path.
.sigma_from_NullDist <- function(y, psi, tol = 1e-10, max.init = 1000) {
  q_inv <- if (!is.null(y$CDFinverse)) {
    y$CDFinverse
  } else {
    function(p) .invert_cdf(y$CDF, y$pdf, p, tol = tol, max.init = max.init)
  }
  (psi[2] - psi[1]) / (q_inv(psi[2]) - q_inv(psi[1]))
}

# Build F_sigma(t) = F(t/sigma) as a new NullDist, capturing sigma by value.
# Only functions that exist in y are scaled; absent ones remain NULL.
.scale_NullDist <- function(y, sigma) {
  CDF_s <- local({
    s <- sigma; cdf <- y$CDF
    function(t) cdf(t / s)
  })
  
  CDFinverse_s <- if (!is.null(y$CDFinverse)) {
    local({
      s <- sigma; qf <- y$CDFinverse
      function(p) {
        result        <- s * qf(p)
        result[p < 0] <- -Inf
        result[p > 1] <-  Inf
        result
      }
    })
  } else NULL
  
  pdf_s <- if (!is.null(y$pdf)) {
    local({ s <- sigma; dens <- y$pdf; function(t) dens(t / s) / s })
  } else NULL
  
  distribution(CDF = CDF_s, CDFinverse = CDFinverse_s, pdf = pdf_s)
}

# Monte Carlo estimate of sigma* for the two-sample case (Eq. 49, DJK 2025).
# Works on a single coordinate vector pair (x_col, y_col).
.sigma_2s_1d <- function(x_col, y_col, psi, nperms) {
  m <- length(x_col);  n <- length(y_col)
  Z <- c(x_col, y_col)
  mean(replicate(nperms, {
    idx <- sample.int(m + n, m, replace = FALSE)
    Xd  <- Z[idx];  Yd <- Z[-idx]
    Qx  <- quantile(Xd, psi, names = FALSE)
    Qy  <- quantile(Yd, psi, names = FALSE)
    max((psi[2] - psi[1]) / diff(Qx),
        (psi[2] - psi[1]) / diff(Qy))
  }))
}

# Append scaling metadata to an htest object and tag the method string.
.attach_scale_info <- function(result, sigma, psi) {
  result$sigma     <- sigma
  result$scale_psi <- psi
  result$method    <- paste0(result$method, " (scale-tuned)")
  result
}


# =============================================================================
# distribution()
# =============================================================================

distribution <- function(CDF, CDFinverse = NULL, pdf = NULL) {
  if (!is.function(CDF))
    stop("`CDF` must be a function.")
  if (!is.null(CDFinverse) && !is.function(CDFinverse))
    stop("`CDFinverse` must be a function or NULL.")
  if (!is.null(pdf) && !is.function(pdf))
    stop("`pdf` must be a function or NULL.")
  
  CDF <- Vectorize(CDF)
  if (!is.null(pdf))
    pdf <- Vectorize(pdf)
  
  if (!is.null(CDFinverse)) {
    raw_CDFinverse <- Vectorize(CDFinverse)
    CDFinverse <- function(p) {
      result        <- raw_CDFinverse(p)
      result[p < 0] <- -Inf
      result[p > 1] <-  Inf
      result
    }
  }
  
  structure(list(CDF = CDF, CDFinverse = CDFinverse, pdf = pdf),
            class = "NullDist")
}

print.NullDist <- function(x, ...) {
  cat("<NullDist>\n")
  cat("  CDF        :", if (!is.null(x$CDF))        "supplied" else "---", "\n")
  cat("  CDFinverse :", if (!is.null(x$CDFinverse)) "supplied"
      else "NULL (root-finding used)", "\n")
  cat("  pdf        :", if (!is.null(x$pdf))        "supplied"
      else "NULL (uniroot used)", "\n")
  invisible(x)
}


# =============================================================================
# Hausdorff_test()
# =============================================================================
# method = "default"  ->  one-sample:           exact p-value
#                         two-sample univariate: Monte Carlo permutation
#                         two-sample bivariate:  Monte Carlo permutation
# method = "exact"    ->  one-sample:           exact p-value
#                         two-sample univariate: full permutation enumeration
#                         two-sample bivariate:  not yet supported (warning + mc)
# method = "mc"       ->  one-sample:           Monte Carlo bootstrap
#                         two-sample:           Monte Carlo permutation
#
# scale_psi  NULL (no scaling) or length-2 numeric c(psi1, psi2) with
#            0 < psi2 < psi1 < 1. When supplied, sigma* is computed and the
#            data (and null distribution) are scaled before the test. Scaling
#            metadata is attached to the returned htest object.
#
# scale_nperms  number of Monte Carlo replications used to estimate sigma* in
#               the two-sample paths. Ignored in the one-sample path.

Hausdorff_test <- function(x, y,
                           method       = "default",
                           nboots       = 2000,
                           tol          = 1e-10,
                           scale_psi    = NULL,
                           scale_nperms = 1000,
                           max.init     = 1000,
                           ...) {
  method <- match.arg(method, choices = c("default", "exact", "mc"))
  UseMethod("Hausdorff_test", y)
}

# ── One-sample: NullDist ──────────────────────────────────────────────────────

Hausdorff_test.NullDist <- function(x, y,
                                    method       = "default",
                                    nboots       = 2000,
                                    tol          = 1e-10,
                                    scale_psi    = NULL,
                                    scale_nperms = 1000,
                                    max.init     = 1000,
                                    ...) {
  method    <- match.arg(method, choices = c("default", "exact", "mc"))
  data_name <- deparse(substitute(x))   # capture before x is overwritten
  
  sigma <- NULL
  if (!is.null(scale_psi)) {
    psi   <- .validate_psi(scale_psi)
    sigma <- .sigma_from_NullDist(y, psi, tol = tol, max.init = max.init)
    x     <- sigma * x
    y     <- .scale_NullDist(y, sigma)
  }
  
  result <- if (method %in% c("default", "exact")) {
    H_test_1s_1d(x          = x,
                 CDF        = y$CDF,
                 CDFinverse = y$CDFinverse,
                 pdf        = y$pdf,
                 tol        = tol)
  } else {
    .Hausdorff_1s_mc(x          = x,
                     CDF        = y$CDF,
                     CDFinverse = y$CDFinverse,
                     nboots     = nboots,
                     tol        = tol)
  }
  
  result$data.name <- data_name   # restore original variable name
  if (!is.null(sigma))
    result <- .attach_scale_info(result, sigma, psi)
  result
}

# ── One-sample: bare CDF function shorthand ───────────────────────────────────

Hausdorff_test.function <- function(x, y,
                                    method       = "default",
                                    nboots       = 2000,
                                    tol          = 1e-10,
                                    scale_psi    = NULL,
                                    scale_nperms = 1000,
                                    max.init     = 1000,
                                    ...) {
  method <- match.arg(method, choices = c("default", "exact", "mc"))
  # Promotes bare CDF to NullDist; .sigma_from_NullDist will error informatively
  # if scale_psi is set but CDFinverse is absent (which it will be here).
  Hausdorff_test.NullDist(x,
                          y            = distribution(CDF = y),
                          method       = method,
                          nboots       = nboots,
                          tol          = tol,
                          scale_psi    = scale_psi,
                          scale_nperms = scale_nperms,
                          max.init     = max.init)
}

# ── Two-sample univariate ─────────────────────────────────────────────────────

Hausdorff_test.numeric <- function(x, y,
                                   method       = "default",
                                   nboots       = 2000,
                                   tol          = 1e-10,
                                   scale_psi    = NULL,
                                   scale_nperms = 1000,
                                   ...) {
  if (!is.numeric(x))
    stop("`x` must be a numeric vector.")
  method    <- match.arg(method, choices = c("default", "exact", "mc"))
  data_name <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  
  sigma <- NULL
  if (!is.null(scale_psi)) {
    psi   <- .validate_psi(scale_psi)
    sigma <- .sigma_2s_1d(x, y, psi, scale_nperms)
    x     <- sigma * x
    y     <- sigma * y
  }
  
  use_exact <- (method == "exact")
  result    <- H_test_2s_1d(x1 = x, x2 = y, nboots = nboots, Exact = use_exact)
  
  result$data.name <- data_name   # restore original variable names
  if (!is.null(sigma))
    result <- .attach_scale_info(result, sigma, psi)
  result
}

# ── Two-sample bivariate: matrix ──────────────────────────────────────────────

Hausdorff_test.matrix <- function(x, y,
                                  method       = "default",
                                  nboots       = 2000,
                                  tol          = 1e-6,
                                  scale_psi    = NULL,
                                  scale_nperms = 1000,
                                  max.init     = 1000,
                                  invariant    = FALSE,
                                  ...) {
  method    <- match.arg(method, choices = c("default", "exact", "mc"))
  data_name <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  
  sigma <- NULL
  if (!is.null(scale_psi)) {
    psi   <- .validate_psi(scale_psi)
    # Column-wise sigma*: apply the univariate formula to each coordinate.
    sigma <- vapply(1:2, function(j)
      .sigma_2s_1d(x[, j], y[, j], psi, scale_nperms), numeric(1))
    x <- sweep(x, 2, sigma, `*`)
    y <- sweep(y, 2, sigma, `*`)
  }
  
  result <- H_test_2s_2d(x = x, y = y, nboots = nboots,
                         Exact = (method == "exact"),
                         invariant = invariant, tol = tol)
  
  result$data.name <- data_name   # restore original variable names
  if (!is.null(sigma))
    result <- .attach_scale_info(result, sigma, psi)
  result
}

# ── Two-sample bivariate: list ────────────────────────────────────────────────

Hausdorff_test.list <- function(x, y,
                                method       = "default",
                                nboots       = 2000,
                                tol          = 1e-6,
                                scale_psi    = NULL,
                                scale_nperms = 1000,
                                max.init     = 1000,
                                invariant    = FALSE,
                                ...) {
  method <- match.arg(method, choices = c("default", "exact", "mc"))
  x <- do.call(cbind, x)
  y <- do.call(cbind, y)
  Hausdorff_test.matrix(x, y,
                        method       = method,
                        nboots       = nboots,
                        tol          = tol,
                        scale_psi    = scale_psi,
                        scale_nperms = scale_nperms,
                        invariant    = invariant,
                        ...)
}


# =============================================================================
# Hausdorff_stat()
# =============================================================================

Hausdorff_stat <- function(x, y, tol = 1e-10, ...) {
  UseMethod("Hausdorff_stat", y)
}

Hausdorff_stat.numeric <- function(x, y, tol = 1e-10, ...) {
  H_stat_2s_1d_tr(a = x, b = y)
}

Hausdorff_stat.NullDist <- function(x, y, tol = 1e-10, ...) {
  H_stat_1s_1d(x  = x,
               CDF = y$CDF,
               pdf = y$pdf,
               tol = tol)
}

Hausdorff_stat.function <- function(x, y, tol = 1e-10, ...) {
  H_stat_1s_1d(x  = x,
               CDF = y,
               tol = tol)
}

Hausdorff_stat.matrix <- function(x, y, tol = 1e-6, invariant = FALSE, ...) {
  if (!invariant) {
    H_stat_2s_2d(x = x, y = y, tol = tol)
  } else {
    max(
      H_stat_2s_2d(x =  x,                    y =  y,                    tol = tol),
      H_stat_2s_2d(x = cbind( x[,1], -x[,2]), y = cbind( y[,1], -y[,2]), tol = tol),
      H_stat_2s_2d(x = cbind(-x[,1],  x[,2]), y = cbind(-y[,1],  y[,2]), tol = tol),
      H_stat_2s_2d(x = cbind(-x[,1], -x[,2]), y = cbind(-y[,1], -y[,2]), tol = tol)
    )
  }
}

Hausdorff_stat.list <- function(x, y, tol = 1e-6, invariant = FALSE, ...) {
  x <- do.call(cbind, x)
  y <- do.call(cbind, y)
  Hausdorff_stat.matrix(x, y, tol = tol, invariant = invariant, ...)
}


# =============================================================================
# Private Monte Carlo one-sample bootstrap (not exported)
# =============================================================================
# Called by Hausdorff_test.NullDist when method = "mc".
# Generates nboots bootstrap samples from F via the quantile transform
# (CDFinverse if supplied, otherwise uniroot inversion of CDF) and returns
# the proportion of bootstrap H statistics exceeding the observed value.

.Hausdorff_1s_mc <- function(x, CDF, CDFinverse = NULL,
                             nboots = 2000, tol = 1e-6, max.init = 1000) {
  n         <- length(x)
  test.stat <- H_stat_1s_1d(x, CDF, tol = tol, max.init = max.init)
  
  if (is.null(CDFinverse)) {
    CDFinverse <- function(p) {
      uniroot(function(t) CDF(t) - p, c(0, 1),
              extendInt = "yes", tol = tol)$root
    }
  }
  
  bigger <- 0L
  for (i in seq_len(nboots)) {
    x0     <- CDFinverse(runif(n))
    bigger <- bigger +
      (H_stat_1s_1d(x0, CDF, tol = tol, max.init = max.init) > test.stat)
  }
  
  names(test.stat) <- "H"
  result <- list(
    p.value     = bigger / nboots,
    method      = "One-sample Hausdorff test",
    statistic   = test.stat,
    alternative = "two-sided",
    data.name   = deparse(substitute(x))
  )
  class(result) <- "htest"
  result
}