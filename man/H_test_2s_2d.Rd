\name{H_test_2s_2d}
\alias{H_test_2s_2d}
\title{Two-sample bivariate Hausdorff goodness-of-fit test with Monte Carlo
permutation p-values}

\description{
Tests the null hypothesis \eqn{H_{0}: F(\bm{x}) = G(\bm{x})}{} for all
\eqn{\bm{x} \in \mathbb{R}^{2}}{} against the two-sided alternative
\eqn{H_{1}: F(\bm{x}) \neq G(\bm{x})}{} for at least one \eqn{\bm{x}}{},
using the Hausdorff distance between the bivariate empirical cdfs of two
independent samples.  The p-value is obtained by Monte Carlo permutation.
The test statistic is computed by \code{\link{H_stat_2s_2d}}.
}

\details{
Given two independent random samples
\eqn{\bm{X} = \{(X_{1,1}, X_{1,2}), \ldots, (X_{m,1}, X_{m,2})\}}{} and
\eqn{\bm{Y} = \{(Y_{1,1}, Y_{1,2}), \ldots, (Y_{n,1}, Y_{n,2})\}}{} from
unknown bivariate cdfs \eqn{F}{} and \eqn{G}{} respectively, the two-sample
Hausdorff (\eqn{H}{}) test statistic is
\eqn{\hat{q} = \mathcal{H}_{m,n} = H(F^{c}_{m}, G^{c}_{n})}{},
computed by \code{\link{H_stat_2s_2d}}.  By Theorem 2.5 of Dimitrova, Jia,
and Kaishev (2026b), this equals the side length of the largest axis-aligned
hypercube that can be inscribed in the region between the two bivariate
empirical cdf surfaces \eqn{F^{c}_{m}}{} and \eqn{G^{c}_{n}}{}.

\bold{Permutation p-value.}  Because the null distribution of
\eqn{\mathcal{H}_{m,n}}{} depends on the unknown \eqn{F}{} and \eqn{G}{},
p-values are obtained by permutation (Proposition 3.5 of Dimitrova, Jia, and
Kaishev 2026b).  When \code{Exact = TRUE} all \eqn{\frac{(m+n)!}{m!\,n!}}{}
splits of the pooled sample are enumerated; if the count exceeds
\eqn{2.147 \times 10^{9}}{} the function falls back to Monte Carlo
automatically.  When \code{Exact = FALSE} the pooled sample is randomly split
for \code{nboots} replications and the p-value is the proportion of
permutation statistics \eqn{\geq \hat{q}}{}, floored at
\eqn{1/(2 \times \code{nboots})}{}.  The asymptotic equivalence between the
permutation p-value and the true p-value of \eqn{\mathcal{H}_{m,n}}{}, and
the control of type I error, follow from Theorems 3.2--3.4 of Dimitrova, Jia,
and Kaishev (2026b); see \code{\link{H_test_2s_1d}} for the full statement.

\bold{Order-invariant statistic.}  When \code{invariant = TRUE} the
order-invariant statistic
\eqn{H^{\circ}_{m,n} = \max_{1 \leq i \leq 4} H(F^{c}_{m,i}, G^{c}_{n,i})}{}
(Equation (17) of Dimitrova, Jia, and Kaishev 2026b) is used for both the
observed value and every permutation replicate.  It maximises over the
\eqn{2^{2} = 4}{} sign-flip orderings of the two coordinates, making the
test invariant to the labelling of coordinate axes.

\bold{Scale tuning.}  For scale-tuned testing pass \code{scale_psi} to
\code{\link{Hausdorff_test}}, which estimates the column-wise optimal scale
vector \eqn{(\sigma^{*}_{1}, \sigma^{*}_{2})}{} and calls this function
internally on the already-scaled samples.
}

\usage{
H_test_2s_2d(x, y, nboots = 2000, Exact = FALSE, invariant = FALSE, tol = 1e-6)
}

\arguments{
\item{x}{a two-column numeric matrix representing the first bivariate sample
of size \eqn{m}{}.}
\item{y}{a two-column numeric matrix representing the second bivariate sample
of size \eqn{n}{}.}
\item{nboots}{a positive integer giving the number of Monte Carlo permutation
replications, used only when \code{Exact = FALSE}.  Defaults to \code{2000}.}
\item{Exact}{a logical value indicating whether an exact permutation p-value
should be computed by enumerating all \eqn{\frac{(m+n)!}{m!\,n!}}{} splits of
the pooled sample.  If \eqn{\frac{(m+n)!}{m!\,n!} > 2.147483647 \times 10^{9}}{}
the function automatically switches to Monte Carlo and emits a message.
Defaults to \code{FALSE}.}
\item{invariant}{logical.  When \code{TRUE}, the order-invariant statistic
\eqn{H^{\circ}_{m,n}}{} is used for the observed value and all permutation
replicates, maximising over the four sign-flip orderings of the two
coordinates (Equation (17) of Dimitrova, Jia, and Kaishev 2026b).
Defaults to \code{FALSE}.}
\item{tol}{a numeric value giving the tolerance for locating omnidirectional
jump vertices in the bivariate empirical cdfs, passed to
\code{\link{H_stat_2s_2d}}.  Defaults to \code{1e-6}.}
}

\value{
An object of class \code{"htest"} with the following components:
\describe{
\item{\code{statistic}}{the observed Hausdorff statistic
\eqn{\hat{q} = \mathcal{H}_{m,n} = H(F^{c}_{m}, G^{c}_{n})}{} (or
\eqn{H^{\circ}_{m,n}}{} when \code{invariant = TRUE}), named \code{"H"}.}
\item{\code{p.value}}{the Monte Carlo permutation p-value, i.e. the proportion
of permutation statistics greater than or equal to \eqn{\hat{q}}{}.  If
zero, replaced by \eqn{1/(2 \times \code{nboots})}{}.}
\item{\code{method}}{a character string indicating the procedure:
\code{"Two-sample bivariate Hausdorff test (Exact)"} or
\code{"Two-sample bivariate Hausdorff test (Monte Carlo permutation)"},
with the suffix \code{", order-invariant"} appended when
\code{invariant = TRUE}.}
\item{\code{alternative}}{the character string \code{"two-sided"}.}
\item{\code{data.name}}{a character string giving the names of the two data
objects.}
}
}

\references{
Dimitrina S. Dimitrova, Yun Jia, Vladimir K. Kaishev (2026a). ``On a One
Sample Goodness-of-fit Test Based on the Hausdorff Metric''. \emph{Submitted}.

Dimitrina S. Dimitrova, Yun Jia, Vladimir K. Kaishev (2026b). ``On a Two-Sample
Multivariate Goodness-of-Fit Test Based on the Hausdorff Metric''.
\emph{Submitted to the Annals of Statistics}.
}

\seealso{
\code{\link{H_stat_2s_2d}},
\code{\link{Hausdorff_test}},
\code{\link{H_test_2s_1d}}.
}

\examples{\donttest{
## Bivariate H test: both samples from N(0,I) (null is true)
set.seed(1)
xm <- matrix(rnorm(100), ncol = 2)
ym <- matrix(rnorm(100), ncol = 2)
H_test_2s_2d(xm, ym, nboots = 1000)

## Bivariate H test: samples from N(0,I) and N(0.5,I) (null is false)
set.seed(2)
xm2 <- matrix(rnorm(100), ncol = 2)
ym2 <- matrix(rnorm(100, mean = 0.5), ncol = 2)
H_test_2s_2d(xm2, ym2, nboots = 1000)

## Exact permutation test for small samples
set.seed(3)
H_test_2s_2d(matrix(rnorm(16), ncol = 2),
             matrix(rnorm(16), ncol = 2), Exact = TRUE)

## Order-invariant statistic
set.seed(4)
H_test_2s_2d(xm, ym, nboots = 1000, invariant = TRUE)

## Scale-tuned test via the unified wrapper (recommended)
set.seed(4)
H_test_2s_2d(xm, ym, nboots = 1000)
Hausdorff_test(xm, ym, nboots = 1000,
               scale_psi = c(0.99, 0.95), scale_nperms = 500)
}}
