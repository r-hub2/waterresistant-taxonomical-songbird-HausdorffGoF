\name{Hausdorff_stat}
\alias{Hausdorff_stat}
\alias{Hausdorff_stat.numeric}
\alias{Hausdorff_stat.NullDist}
\alias{Hausdorff_stat.function}
\alias{Hausdorff_stat.matrix}
\alias{Hausdorff_stat.list}
\title{Unified one-sample and two-sample Hausdorff test statistic}

\description{
An S3 generic that computes the Hausdorff (\eqn{H}{}) test statistic
without a p-value, dispatching on the class of \code{y} in the same way as
\code{\link{Hausdorff_test}}:

\describe{
\item{\code{y} numeric}{Two-sample statistic \eqn{\mathcal{H}_{m,n} = H(F^{c}_{m}, G^{c}_{n})}{}
via \code{\link{H_stat_2s_1d_tr}}.}
\item{\code{y} \code{NullDist}}{One-sample statistic \eqn{H(F^{c}, F^{c}_{n})}{} via
\code{\link{H_stat_1s_1d}}.}
\item{\code{y} function}{One-sample shorthand; the bare cdf is passed directly to
\code{\link{H_stat_1s_1d}} with
\code{pdf = NULL}.}
\item{\code{y} matrix or list}{Two-sample bivariate statistic \eqn{H(F^{c}_{m}, G^{c}_{n})}{} via
\code{\link{H_stat_2s_2d}}.  Both \code{x}
and \code{y} must be two-column matrices (or lists of two numeric vectors)
of sizes \eqn{m}{} and \eqn{n}{}.  When \code{invariant = TRUE} the
order-invariant statistic
\eqn{H^{\circ}_{m,n} = \max_{1 \leq i \leq 4} H(F^{c}_{m,i}, G^{c}_{n,i})}{}
is returned by maximising over the \eqn{2^{2} = 4}{} sign-flip orderings of
the two coordinates (Equation (17) of Dimitrova, Jia, and Kaishev 2026b).
List input is converted to a two-column matrix and the same path is followed.}
}
}

\usage{
Hausdorff_stat(x, y, tol = 1e-10, ...)

\method{Hausdorff_stat}{numeric}(x, y, tol = 1e-10, ...)

\method{Hausdorff_stat}{NullDist}(x, y, tol = 1e-10, ...)

\method{Hausdorff_stat}{function}(x, y, tol = 1e-10, ...)

\method{Hausdorff_stat}{matrix}(x, y, tol = 1e-6, invariant = FALSE, ...)

\method{Hausdorff_stat}{list}(x, y, tol = 1e-6, invariant = FALSE, ...)
}

\arguments{
\item{x}{a numeric vector (univariate case) or a two-column numeric matrix or list
of two numeric vectors (bivariate case).}
\item{y}{one of: (i) a numeric vector (two-sample univariate), (ii) a
\code{\link{distribution}} object (one-sample, full
control), (iii) a bare vectorised \R function for the null cdf
(one-sample shorthand), or (iv) a two-column numeric matrix or list of two
numeric vectors (two-sample bivariate).}
\item{tol}{a numeric value giving the tolerance for the root-finding step (one-sample
path, default \code{1e-10}) or for locating omnidirectional jump vertices
(bivariate path, default \code{1e-6}).}
\item{invariant}{logical; bivariate path only.  When \code{TRUE} the order-invariant
statistic \eqn{H^{\circ}_{m,n}}{} is returned by maximising over the four
sign-flip orderings of the two coordinates.  Defaults to \code{FALSE}.}
\item{...}{further arguments (currently unused).}
}

\value{
A single numeric value: the observed Hausdorff statistic \eqn{\hat{q}}{}.
}

\seealso{
\code{\link{distribution}},
\code{\link{Hausdorff_test}},
\code{\link{H_stat_1s_1d}},
\code{\link{H_stat_2s_1d_tr}},
\code{\link{H_stat_2s_2d}}.
}

\examples{
## --- Univariate two-sample statistic ----------------------------------------
set.seed(1)
Hausdorff_stat(rnorm(40), rnorm(40))

## --- One-sample statistic: full distribution object -------------------------
set.seed(2)
x <- rexp(50, rate = 1)
null_exp <- distribution(CDF = function(t) pexp(t, rate = 1),
                         pdf = function(t) dexp(t, rate = 1))
Hausdorff_stat(x, null_exp)

## --- One-sample statistic: bare CDF shorthand -------------------------------
set.seed(3)
Hausdorff_stat(rnorm(60), pnorm)

## --- Bivariate two-sample statistic (standard ordering) ---------------------
set.seed(4)
x <- matrix(rnorm(100), ncol = 2)
y <- matrix(rnorm(100), ncol = 2)
Hausdorff_stat(x, y)

## --- Bivariate two-sample statistic (order-invariant) -----------------------
Hausdorff_stat(x, y, invariant = TRUE)

## --- Bivariate two-sample statistic: list input -----------------------------
set.seed(5)
x3 <- list(rnorm(30), rnorm(30))
y3 <- list(rnorm(30), rnorm(30))
Hausdorff_stat(x3, y3)
}
