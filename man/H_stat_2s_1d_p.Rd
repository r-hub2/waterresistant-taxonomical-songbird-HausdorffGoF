\name{H_stat_2s_1d_p}
\alias{H_stat_2s_1d_p}
\title{C++ implementation of the two-sample Hausdorff distance via the projection method}

\description{
A C++ implementation (via \pkg{Rcpp}) of the projection method for computing the
Hausdorff distance \eqn{H(F^{c}_{m}, G^{c}_{n})}{} between the planar curves of
the empirical cdfs of two independent univariate samples
\eqn{\{x_{1}, \ldots, x_{m}\}}{} and \eqn{\{y_{1}, \ldots, y_{n}\}}{}.
See Appendix E of Dimitrova, Jia, and Kaishev (2026b) for a numerical accuracy
and timing comparison with
\code{\link{H_stat_2s_1d_tr}}.
}

\details{
\code{\link{H_stat_2s_1d_p}} implements the
\bold{projection method} of Section 3.1 of Dimitrova, Jia, and Kaishev (2026b)
and proceeds in five steps identical to those described for the univariate case
in the bivariate function
\code{\link{H_stat_2s_2d}}, specialised to
\eqn{k = 1}{}.  The overall computation is \eqn{O(m + n)}{} in the number of
distinct observations.
}

\usage{
H_stat_2s_1d_p(a, b)
}

\arguments{
\item{a}{a numeric vector of data sample values \eqn{\{x_{1}, \ldots, x_{m}\}}{}.}
\item{b}{a numeric vector of data sample values \eqn{\{y_{1}, \ldots, y_{n}\}}{}.}
}

\value{
A single numeric value: the Hausdorff distance \eqn{H(F^{c}_{m}, G^{c}_{n})}{}.
}

\references{
Dimitrina S. Dimitrova, Yun Jia, Vladimir K. Kaishev (2026b). ``On a Two-Sample
Multivariate Goodness-of-Fit Test Based on the Hausdorff Metric''.
\emph{Submitted}.
}

\seealso{
\code{\link{H_stat_2s_1d_tr}} for the
faster transformation-method implementation.
}

\examples{
set.seed(1)
a <- rnorm(50)
b <- rnorm(50)
H_stat_2s_1d_p(a, b)

## Verify agreement with the transformation-method C++ implementation
H_stat_2s_1d_tr(a, b)
}
