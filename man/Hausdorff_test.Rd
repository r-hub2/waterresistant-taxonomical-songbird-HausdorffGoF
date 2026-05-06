\name{Hausdorff_test}
\alias{Hausdorff_test}
\alias{Hausdorff_test.NullDist}
\alias{Hausdorff_test.function}
\alias{Hausdorff_test.numeric}
\alias{Hausdorff_test.matrix}
\alias{Hausdorff_test.list}
\title{Unified one-sample and two-sample Hausdorff goodness-of-fit test}

\description{
An S3 generic that performs the Hausdorff (\eqn{H}{}) goodness-of-fit test
and returns an object of class \code{"htest"}, dispatching on the class of
\code{y}.

\bold{Dispatch paths.}
\describe{
\item{\code{y} \code{NullDist}}{One-sample test.
\code{\link{H_test_1s_1d}} is called when
\code{method} is \code{"default"} or \code{"exact"}, giving the exact
rectangle-probability p-value.  When \code{method = "mc"}, a Monte Carlo
bootstrap p-value is computed instead: \code{nboots} samples of size
\eqn{n}{} are drawn from \eqn{F}{} via the quantile transform and the
proportion of bootstrap \eqn{H}{} statistics exceeding the observed value
is returned.}
\item{\code{y} function}{One-sample shorthand.  The bare CDF is promoted to a \code{NullDist} with
\code{CDFinverse = NULL} and \code{pdf = NULL}, then dispatched as above.
Activating \code{scale_psi} with a bare function is not supported because
\eqn{\sigma^{*}}{} requires a quantile function; an informative error
directs the user to \code{\link{distribution}}.}
\item{\code{y} numeric}{Two-sample univariate test via
\code{\link{H_test_2s_1d}}.
\code{"default"} and \code{"mc"} use Monte Carlo permutation;
\code{"exact"} enumerates all \eqn{\frac{(m+n)!}{m!\,n!}}{} permutations.}
\item{\code{y} matrix or list}{Two-sample bivariate test via
\code{\link{H_test_2s_2d}}.  Monte Carlo permutation is used by default;
\code{method = "exact"} enumerates all \eqn{\frac{(m+n)!}{m!\,n!}}{} splits
(automatically falling back to Monte Carlo if the count exceeds
\eqn{2.147 \times 10^{9}}{}).  List input is converted to a two-column matrix
before proceeding.}
}

\bold{Resolution of} \code{method} \bold{by path.}
\describe{
\item{\code{"default"}}{One-sample: exact rectangle-probability p-value.
Two-sample (univariate and bivariate): Monte Carlo permutation.}
\item{\code{"exact"}}{One-sample: exact rectangle-probability p-value (same as \code{"default"}).
Two-sample univariate: full enumeration of all permutations.
Two-sample bivariate: full enumeration of all permutations (falls back to
Monte Carlo automatically if sample sizes are too large).}
\item{\code{"mc"}}{One-sample: Monte Carlo bootstrap p-value.  Bootstrap samples are drawn
from \eqn{F}{} via \code{CDFinverse} if supplied in the \code{NullDist}
object, otherwise by \code{uniroot} inversion of \code{CDF}.
Two-sample (univariate and bivariate): Monte Carlo permutation.}
}

\bold{Scale tuning} (\code{scale_psi})\bold{.}
When \code{scale_psi} is supplied, \eqn{\sigma^{*}}{} is computed before
the test, the data are rescaled, and the test is run on the scaled inputs.

\emph{One-sample} (\code{y} is \code{NullDist}): \eqn{\sigma^{*}}{} uses
the closed-form formula of Proposition 13 of Dimitrova, Jia, and Kaishev
(2026a),
\deqn{\sigma^{*} = \frac{\psi_{1} - \psi_{2}}{F^{-1}(\psi_{1}) - F^{-1}(\psi_{2})}.}{}
The quantile \eqn{F^{-1}}{} is evaluated by priority: \code{CDFinverse}
(if supplied in the \code{NullDist} object) \eqn{>}{} Newton--Raphson (if
\code{pdf} is supplied) \eqn{>}{} \code{uniroot} applied to \code{CDF}.
The Newton--Raphson iteration is bounded by \code{max.init} steps.  The
sample is replaced by \eqn{\sigma^{*} x}{} and the null distribution by
\eqn{F_{\sigma^{*}}(t) = F(t/\sigma^{*})}{}.  Only functions that exist in
the original \code{NullDist} object are scaled; absent ones remain
\code{NULL}.

\emph{Two-sample univariate} (\code{y} is numeric): \eqn{\sigma^{*}}{} is
estimated from the pooled sample by averaging over \code{scale_nperms}
random splits (Equation (49) of Dimitrova, Jia, and Kaishev 2026b),
\deqn{\sigma^{*} \approx \frac{1}{B}\sum_{b=1}^{B}
\max\!\left(
  \frac{\psi_{1}-\psi_{2}}{Q_{X^{\dagger}_{m}}^{(b)}(\psi_{1}) -
    Q_{X^{\dagger}_{m}}^{(b)}(\psi_{2})},\;
  \frac{\psi_{1}-\psi_{2}}{Q_{Y^{\dagger}_{n}}^{(b)}(\psi_{1}) -
    Q_{Y^{\dagger}_{n}}^{(b)}(\psi_{2})}
\right).}{}
Both samples are scaled by the scalar \eqn{\sigma^{*}}{}.

\emph{Two-sample bivariate} (\code{y} is matrix or list): the univariate
formula is applied column-wise, yielding a length-2 vector
\eqn{(\sigma^{*}_{1}, \sigma^{*}_{2})}{}.  Column \eqn{j}{} of each matrix
is scaled by \eqn{\sigma^{*}_{j}}{}.

In all cases, \code{scale_psi} and \eqn{\sigma^{*}}{} are attached to the
returned \code{"htest"} object as \code{$scale_psi} and \code{$sigma},
and the string \code{"(scale-tuned)"} is appended to \code{$method}.
}

\usage{
Hausdorff_test(x, y, method = "default", nboots = 2000, tol = 1e-10,
               scale_psi = NULL, scale_nperms = 1000, max.init = 1000, ...)

\method{Hausdorff_test}{NullDist}(x, y, method = "default", nboots = 2000,
                                   tol = 1e-10, scale_psi = NULL,
                                   scale_nperms = 1000, max.init = 1000, ...)

\method{Hausdorff_test}{function}(x, y, method = "default", nboots = 2000,
                                   tol = 1e-10, scale_psi = NULL,
                                   scale_nperms = 1000, max.init = 1000, ...)

\method{Hausdorff_test}{numeric}(x, y, method = "default", nboots = 2000,
                                  tol = 1e-10, scale_psi = NULL,
                                  scale_nperms = 1000, ...)

\method{Hausdorff_test}{matrix}(x, y, method = "default", nboots = 2000,
                                 tol = 1e-6, scale_psi = NULL,
                                 scale_nperms = 1000, max.init = 1000,
                                 invariant = FALSE, ...)

\method{Hausdorff_test}{list}(x, y, method = "default", nboots = 2000,
                               tol = 1e-6, scale_psi = NULL,
                               scale_nperms = 1000, max.init = 1000,
                               invariant = FALSE, ...)
}

\arguments{
\item{x}{a numeric vector (one-sample or two-sample univariate), or a two-column
numeric matrix or list of two numeric vectors (bivariate).}
\item{y}{one of: (i) a \code{\link{distribution}} object, (ii) a
bare CDF function, (iii) a numeric vector, or (iv) a two-column numeric
matrix or list.  See Dispatch paths above.}
\item{method}{\code{"default"}, \code{"exact"}, or \code{"mc"}.  See the resolution
table in the Description.  Defaults to \code{"default"}.}
\item{nboots}{a positive integer: number of Monte Carlo replications for the test
p-value (one-sample bootstrap or two-sample permutation).
Defaults to \code{2000}.}
\item{tol}{numeric tolerance for root-finding in the one-sample statistic and p-value
computation (default \code{1e-10}), or for locating omnidirectional jump
vertices in the bivariate statistic (default \code{1e-6}).
Also used when \code{scale_psi} is supplied and Newton--Raphson inversion
is needed to compute \eqn{\sigma^{*}}{}.}
\item{scale_psi}{\code{NULL} (no scaling, default) or a length-2 numeric vector
\code{c(psi1, psi2)} with \eqn{0 < \psi_{2} < \psi_{1} < 1}{}.  The two
values may be supplied in any order; they are sorted internally.
Recommended choices (Dimitrova, Jia, and Kaishev 2026a, 2026b):
\describe{
\item{Right tail}{\code{c(0.99, 0.95)}.}
\item{Body}{\code{c(0.70, 0.40)} or \code{c(0.60, 0.30)}.}
\item{Left tail}{\code{c(0.05, 0.01)}.}
}}
\item{scale_nperms}{a positive integer: number of Monte Carlo splits used to estimate
\eqn{\sigma^{*}}{} in the two-sample paths.  Ignored when
\code{scale_psi = NULL} or in the one-sample path.  Defaults to
\code{1000}.}
\item{max.init}{a positive integer: maximum number of Newton--Raphson iterations when
inverting \eqn{F}{} numerically to compute \eqn{\sigma^{*}}{} in the
one-sample path.  Ignored when \code{CDFinverse} is supplied in the
\code{NullDist} object, or when \code{scale_psi = NULL}.
Defaults to \code{1000}.}
\item{invariant}{logical; bivariate path only.  When \code{TRUE}, the order-invariant
statistic \eqn{H^{\circ}_{m,n}}{} is used throughout (observed value and
all permutation replicates), maximising over the \eqn{2^{2} = 4}{}
sign-flip orderings of the two coordinates (Equation (17) of Dimitrova,
Jia, and Kaishev 2026b).  Defaults to \code{FALSE}.}
\item{...}{further arguments (currently unused).}
}

\value{
An object of class \code{"htest"} with components:
\describe{
\item{\code{statistic}}{the observed Hausdorff statistic \eqn{\hat{q}}{} (computed on the scaled
data when \code{scale_psi} is supplied), named \code{"H"}.}
\item{\code{p.value}}{the p-value computed according to the active path and \code{method}.}
\item{\code{method}}{a character string identifying the procedure.  The suffix
\code{"(scale-tuned)"} is appended when \code{scale_psi} is supplied.}
\item{\code{alternative}}{the character string \code{"two-sided"}.}
\item{\code{data.name}}{a character string giving the names of the data objects.}
\item{\code{sigma}}{(only when \code{scale_psi} is supplied) the computed \eqn{\sigma^{*}}{}:
a scalar for the one-sample and two-sample univariate paths; a length-2
vector for the bivariate path.}
\item{\code{scale_psi}}{(only when \code{scale_psi} is supplied) the validated sorted
\code{c(psi_low, psi_high)} vector used to compute \eqn{\sigma^{*}}{}.}
}
}

\references{
Dimitrina S. Dimitrova, Yun Jia, Vladimir K. Kaishev (2026a). ``On a One
Sample Goodness-of-fit Test Based on the Hausdorff Metric''.
\emph{Submitted}.

Dimitrina S. Dimitrova, Yun Jia, Vladimir K. Kaishev (2026b). ``On a
Two-Sample Multivariate Goodness-of-Fit Test Based on the Hausdorff
Metric''. \emph{Submitted}.

Dimitrina S. Dimitrova, Vladimir K. Kaishev, Senren Tan (2020).
``Computing the Kolmogorov-Smirnov Distribution When the Underlying CDF is
Purely Discrete, Mixed or Continuous''. \emph{Journal of Statistical
Software}, \bold{95}(10): 1--42. doi:10.18637/jss.v095.i10.
}

\seealso{
\code{\link{distribution}},
\code{\link{Hausdorff_stat}},
\code{\link{H_test_1s_1d}},
\code{\link{H_test_c_cdf}},
\code{\link{H_test_2s_1d}},
\code{\link{H_stat_2s_2d}}.
}

\examples{
## ---- One-sample, no scaling ------------------------------------------------
set.seed(1)
x      <- rexp(50, rate = 1)
null_e <- distribution(CDF        = function(t) pexp(t, rate = 1),
                       CDFinverse = function(p) qexp(p, rate = 1),
                       pdf        = function(t) dexp(t, rate = 1))
Hausdorff_test(x, null_e)

## ---- One-sample, scale-tuned: right tail -----------------------------------
res <- Hausdorff_test(x, null_e, scale_psi = c(0.99, 0.95))
res$method    # "... (scale-tuned)"
res$sigma
res$scale_psi


\donttest{
## ---- One-sample, Monte Carlo p-value, scale-tuned: body --------------------
Hausdorff_test(x, null_e, method = "mc", nboots = 1000,
               scale_psi = c(0.70, 0.40))
}

## ---- One-sample, only CDF and pdf supplied (Newton-Raphson sigma*) ---------
null_cdf_only <- distribution(CDF = function(t) pexp(t, rate = 1),
                              pdf = function(t) dexp(t, rate = 1))
Hausdorff_test(x, null_cdf_only, scale_psi = c(0.99, 0.95))

## ---- Two-sample univariate, no scaling -------------------------------------
set.seed(2)
x1 <- rnorm(40);  x2 <- rnorm(40)
Hausdorff_test(x1, x2)


## ---- Two-sample univariate, scale-tuned ------------------------------------
res2 <- Hausdorff_test(x1, x2, scale_psi = c(0.99, 0.95), scale_nperms = 500)
res2$sigma

## ---- Two-sample univariate, exact permutation, small samples ---------------
set.seed(3)
Hausdorff_test(rnorm(8), rnorm(8), method = "exact")

\donttest{
## ---- Two-sample bivariate, no scaling --------------------------------------
set.seed(4)
xm <- matrix(rnorm(100), ncol = 2)
ym <- matrix(rnorm(100, mean = 0.5), ncol = 2)
Hausdorff_test(xm, ym, nboots = 1000)

## ---- Two-sample bivariate, scale-tuned (column-wise sigma*) ----------------
res3 <- Hausdorff_test(xm, ym, nboots = 1000,
                       scale_psi = c(0.70, 0.40), scale_nperms = 500)
res3$sigma   # length-2 vector, one sigma* per coordinate

## ---- Two-sample bivariate, order-invariant statistic -----------------------
Hausdorff_test(xm, ym, nboots = 1000, invariant = TRUE)
}}
