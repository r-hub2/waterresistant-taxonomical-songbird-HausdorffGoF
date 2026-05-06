\name{H_stat_2s_1d_tr}
\alias{H_stat_2s_1d_tr}
\title{C++ implementation of the two-sample Hausdorff distance via the transformation method}

\description{
A C++ implementation (via \pkg{Rcpp}) of the transformation method for computing
the Hausdorff distance \eqn{H(F^{c}_{m}, G^{c}_{n})}{} between the planar curves
of the empirical cdfs of two independent univariate samples
\eqn{\{x_{1}, \ldots, x_{m}\}}{} and \eqn{\{y_{1}, \ldots, y_{n}\}}{}.  This is
the faster of the two exported C++ routines; see Appendix E of Dimitrova, Jia, and
Kaishev (2026b) for timing comparisons.
}

\details{
\code{\link{H_stat_2s_1d_tr}} implements
the \bold{transformation method} of Section 3.1 of Dimitrova, Jia, and Kaishev
(2026b). The plane is rotated by \eqn{\pi/4}{} via the linear map
\eqn{\mathrm{Tr}(x, z) = (x - z,\; x + z)}{} (Equation (35), ibid.), after
which the Hausdorff distance equals half the supremum of the absolute vertical
difference between the two rotated staircase curves (Equation (36), ibid.). The
algorithm proceeds in four steps.

\bold{Step 1 (signed-increment encoding).}
For each sample, the sorted distinct values and their relative jump heights are
extracted.  Both samples are shifted by subtracting the joint minimum.  Each
staircase curve is encoded as an interleaved signed-increment vector.

\bold{Step 2 (rotated coordinates via cumulative sums).}
The \eqn{\pi/4}{} rotation is performed implicitly: the \eqn{z'}{}-coordinate
(\eqn{x + z}{}) and the \eqn{-x'}{}-coordinate (\eqn{z - x}{}) of every
staircase vertex are derived from the encoded vectors via two cumulative sums.

\bold{Step 3 (greedy pointer sweep).}
A single left-to-right pass through both encoded curves simultaneously identifies,
for each vertex of one curve, the segment of the other curve whose
\eqn{\lambda}{}-projection interval contains it.

\bold{Step 4 (distance formula and return value).}
At each vertex, the candidate Hausdorff distance is computed in closed form from
the parity of the current pointer position.  The Hausdorff distance is
\deqn{H(F^{c}_{m}, G^{c}_{n}) = \frac{1}{2} \max_{i}\, |H_{\mathrm{temp}}[i]|.}{}
The entire computation is \eqn{O(m + n)}{} in the number of distinct observations.
}

\usage{
H_stat_2s_1d_tr(a, b)
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
\code{\link{H_stat_2s_1d_p}} for the C++
projection-method implementation (same result, slower).
}

\examples{
set.seed(1)
a <- rnorm(50)
b <- rnorm(50)
H_stat_2s_1d_tr(a, b)

## Verify agreement with the projection-method C++ implementation
H_stat_2s_1d_p(a, b)
}
