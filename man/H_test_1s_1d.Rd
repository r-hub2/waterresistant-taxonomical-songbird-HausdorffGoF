\name{H_test_1s_1d}
\alias{H_test_1s_1d}
\title{One-sample Hausdorff goodness-of-fit test}

\description{
A thin wrapper that first computes the one-sample Hausdorff test statistic
\eqn{\hat{q} = H(F^{c}, F^{c}_{n})}{} via
\code{\link{H_stat_1s_1d}}, then delegates
the p-value computation to
\code{\link{H_test_c_cdf}}.  The split into
two functions allows the p-value to be computed independently when
\eqn{\hat{q}}{} is already known.
}

\usage{
H_test_1s_1d(x, CDF, CDFinverse = NULL, pdf = NULL, tol = 1e-10,
             max.init = 1000)
}

\arguments{
\item{x}{a numeric vector of data sample values \eqn{\{x_{1}, \ldots, x_{n}\}}{}.}
\item{CDF}{a vectorised \R function for the null cdf \eqn{F(x)}{}.}
\item{CDFinverse}{(optional) a vectorised \R function for \eqn{F^{-1}(p)}{}, extended so
that \eqn{F^{-1}(p) = -\infty}{} for \eqn{p < 0}{} and
\eqn{F^{-1}(p) = +\infty}{} for \eqn{p > 1}{}.  When supplied, the
boundary values \eqn{a_{i}}{} and \eqn{b_{i}}{} are computed directly in
\code{\link{H_test_c_cdf}}.
Defaults to \code{NULL}.}
\item{pdf}{(optional) a vectorised \R function for the null density \eqn{f(x)}{}.
When supplied, Newton--Raphson is used in
\code{\link{H_stat_1s_1d}} to solve the
equation \eqn{x + F(x) = \lambda}{}.
Defaults to \code{NULL}.}
\item{tol}{tolerance for root-finding, passed to both
\code{\link{H_stat_1s_1d}} and
\code{\link{H_test_c_cdf}}.
Defaults to \code{1e-10}.}
\item{max.init}{maximum number of Newton--Raphson iterations, passed to
\code{\link{H_stat_1s_1d}}.
Used only when \code{pdf} is supplied. Defaults to \code{1000}.}
}

\value{
An object of class \code{"htest"} with the following components:
\describe{
\item{\code{statistic}}{the observed Hausdorff statistic
\eqn{\hat{q} = H(F^{c}, F^{c}_{n})}{}, named \code{"H"}.}
\item{\code{p.value}}{the exact p-value \eqn{P(H(F^{c}, F^{c}_{n}) > \hat{q})}{}.}
\item{\code{method}}{the character string \code{"One-sample Hausdorff test"}.}
\item{\code{alternative}}{the character string \code{"two-sided"}.}
\item{\code{data.name}}{a character string giving the name of the data object.}
}
}

\seealso{
\code{\link{H_stat_1s_1d}},
\code{\link{H_test_c_cdf}},
\code{\link{Hausdorff_test}}.
}

\examples{
set.seed(1)
x <- rexp(50, rate = 1)
H_test_1s_1d(x,
             CDF        = function(t) pexp(t, rate = 1),
             CDFinverse = function(p) qexp(p, rate = 1),
             pdf        = function(t) dexp(t, rate = 1))
}
