\name{H_test_c_cdf}
\alias{H_test_c_cdf}
\title{Exact complementary cdf of the one-sample Hausdorff statistic}

\description{
Computes the exact complementary cdf
\eqn{P\!\left(H(F^{c},F_n^{c})>q\right)}{}
of the one-sample Hausdorff goodness-of-fit statistic for a sample of size
\eqn{n}{}. For a continuous null cdf \eqn{F}{}, the rectangle-probability
representation of Theorem 4 of Dimitrova, Jia, and Kaishev (2026a) shows that
\deqn{
P\!\left(H(F^{c},F_n^{c})>q\right)
=
1-P(a_i\le U_{(i)}\le b_i,\;1\le i\le n),
}{
P(H(F^c,F_n^c)>q)=1-P(a_i<=U_(i)<=b_i,1<=i<=n),
}
where \eqn{U_{(i)}}{} are the order statistics of an i.i.d.\ \eqn{U(0,1)}{}
sample and
\deqn{
a_i = F\!\left(F^{-1}\!\left(\frac{i}{n}-q\right)-q\right),
\qquad
b_i = F\!\left(F^{-1}\!\left(\frac{i-1}{n}+q\right)+q\right).
}{
a_i = F(F^(-1)(i/n-q)-q), b_i = F(F^(-1)((i-1)/n+q)+q).
}
The function evaluates this rectangle-probability representation numerically
via \code{KSgeneral::ks_c_cdf_Rcpp()} after constructing the boundary
vectors required by that routine.

This is the low-level p-value engine used by
\code{\link{H_test_1s_1d}} and is useful when
the observed Hausdorff statistic \eqn{q}{} has already been computed.
}

\details{
For \eqn{q \in (0,1)}{}, the function constructs two boundary vectors
\code{f_a} and \code{f_b} of length \eqn{n}{}.  If \code{CDFinverse} is
available, these are computed directly from the closed-form expressions for
\eqn{a_i}{} and \eqn{b_i}{}.  Otherwise, the code solves
\eqn{F(x) = i/n - q}{} and \eqn{F(x) = (i-1)/n + q}{} numerically and maps
the resulting roots through \eqn{F(x - q)}{} and \eqn{F(x + q)}{},
respectively.  When \eqn{i/n - q < 0}{} the entry \eqn{a_i = 0}{} is assigned
directly; when \eqn{(i-1)/n + q > 1}{} the entry \eqn{b_i = 1}{} is assigned
directly, avoiding unnecessary root-finding.

The two boundary vectors are written to a temporary file
\file{Boundary_Crossing_Time.txt}, then passed to
\code{KSgeneral::ks_c_cdf_Rcpp(n)}, which evaluates the required
rectangle probability in \eqn{O(n^{2} \log n)}{} time.  The file is deleted
after use.

The trivial boundary cases \code{q >= 1} (returns 1) and \code{q <= 0}
(returns 0) are handled analytically before any computation.
}

\usage{
H_test_c_cdf(q, n, CDF, CDFinverse = NULL, pdf = NULL, tol = 1e-10,
             max.init = 1000)
}

\arguments{
\item{q}{a single numeric value giving the observed Hausdorff statistic.}
\item{n}{a positive integer, the sample size.}
\item{CDF}{a vectorised \R function implementing the null cdf \eqn{F(x)}{}.}
\item{CDFinverse}{(optional) a vectorised \R function implementing \eqn{F^{-1}(p)}{}.
When supplied, the boundary values are computed directly:
\deqn{
a_i = F\!\left(F^{-1}\!\left(\frac{i}{n}-q\right)-q\right),
\qquad
b_i = F\!\left(F^{-1}\!\left(\frac{i-1}{n}+q\right)+q\right),
}{
a_i = F(F^(-1)(i/n-q)-q), b_i = F(F^(-1)((i-1)/n+q)+q).
}
where \eqn{F^{-1}}{} is understood in the extended sense with
\eqn{F^{-1}(y) = -\infty}{} for \eqn{y \le 0}{} and
\eqn{F^{-1}(y) = +\infty}{} for \eqn{y \ge 1}{}.  Defaults to \code{NULL}.}
\item{pdf}{(optional) a vectorised \R function for the null density \eqn{f(x)}{}.
Used only when \code{CDFinverse} is \code{NULL}: if \code{pdf} is
supplied, the required quantiles are obtained by Newton--Raphson iteration;
otherwise by \code{uniroot} with \code{extendInt = "yes"}.
Defaults to \code{NULL}.}
\item{tol}{numerical tolerance used in the root-finding step. Defaults to
\code{1e-10}.}
\item{max.init}{maximum number of Newton--Raphson iterations allowed at each boundary
point. Defaults to \code{1000}.}
}

\value{
A single numeric value equal to the exact complementary cdf
\eqn{P(H(F^{c},F_n^{c})>q)}{}.
}

\references{
Dimitrina S. Dimitrova, Yun Jia, Vladimir K. Kaishev (2026a).
``On a One Sample Goodness-of-fit Test Based on the Hausdorff Metric''.
\emph{Submitted}.

Dimitrina S. Dimitrova, Vladimir K. Kaishev, Senren Tan (2020).
``Computing the Kolmogorov-Smirnov Distribution When the Underlying CDF is
Purely Discrete, Mixed or Continuous''. \emph{Journal of Statistical
Software}, \bold{95}(10), 1--42. doi:10.18637/jss.v095.i10.
}

\seealso{
\code{\link{H_test_1s_1d}},
\code{\link{H_stat_1s_1d}},
\code{\link{Hausdorff_test}}.
}

\examples{
set.seed(1)
x <- rexp(50, rate = 1)
q <- H_stat_1s_1d(x   = x,
                  CDF = function(t) pexp(t, rate = 1),
                  pdf = function(t) dexp(t, rate = 1))
H_test_c_cdf(q          = q,
             n          = length(x),
             CDF        = function(t) pexp(t, rate = 1),
             CDFinverse = function(p) qexp(p, rate = 1),
             pdf        = function(t) dexp(t, rate = 1))
}
