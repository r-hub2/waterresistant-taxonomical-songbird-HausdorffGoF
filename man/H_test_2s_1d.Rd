\name{H_test_2s_1d}
\alias{H_test_2s_1d}
\title{Two-sample Hausdorff goodness-of-fit test with exact or Monte Carlo permutation p-values}

\description{
Tests the null hypothesis \eqn{H_{0}: F(x) = G(x)}{} for all \eqn{x}{} against
the two-sided alternative \eqn{H_{1}: F(x) \neq G(x)}{} for at least one
\eqn{x}{}, using the Hausdorff distance between the empirical cdfs of two
independent univariate samples.  The p-value is obtained either by exact
enumeration of all \eqn{\frac{(m+n)!}{m!n!}}{} permutations of the pooled
sample or by Monte Carlo permutation.  The test statistic is computed by
\code{\link{H_stat_2s_1d_tr}}.
}

\details{
Given two independent random samples
\eqn{\{x_{1}, \ldots, x_{m}\}}{} and \eqn{\{y_{1}, \ldots, y_{n}\}}{} from
unknown univariate cdfs \eqn{F}{} and \eqn{G}{} respectively, the two-sample
Hausdorff (\eqn{H}{}) test statistic is
\eqn{\hat{q} = \mathcal{H}_{m,n} = H(F^{c}_{m}, G^{c}_{n})}{},
computed by
\code{\link{H_stat_2s_1d_tr}}.  By
Theorem 2.5 of Dimitrova, Jia, and Kaishev (2026b), this equals the side length
of the largest axis-aligned square that can be inscribed in the region between
the two empirical cdf staircase curves \eqn{F^{c}_{m}}{} and \eqn{G^{c}_{n}}{}.

\bold{Permutation p-value.} Because the distribution of \eqn{\mathcal{H}_{m,n}}{}
depends on the unknown \eqn{F}{} and \eqn{G}{}, p-values are computed via a
permutation argument (Proposition 3.5 of Dimitrova, Jia, and Kaishev 2026b). Let
\eqn{Z_{m+n} = \{x_{1},\ldots,x_{m},y_{1},\ldots,y_{n}\}}{} be the pooled
sample. The permutation Hausdorff statistic
\eqn{\mathcal{H}^{\dagger}_{m,n}}{} is defined on all
\eqn{C = \frac{(m+n)!}{m!n!}}{} splits of \eqn{Z_{m+n}}{} into groups of
sizes \eqn{m}{} and \eqn{n}{}, and its conditional distribution given
\eqn{Z_{m+n}}{} is:
\deqn{P(\mathcal{H}^{\dagger}_{m,n} > q \mid Z_{m+n})
  = \frac{1}{C} \sum_{i=1}^{C}
    \mathbf{1}\!\left\{H(F^{c\dagger}_{m}(\pi_{i}), G^{c\dagger}_{n}(\pi_{i})) > q\right\},}{}
where \eqn{F^{c\dagger}_{m}(\pi_{i})}{} and \eqn{G^{c\dagger}_{n}(\pi_{i})}{}
are the empirical cdf planar curves of the \eqn{i}{}-th permuted sub-samples.

When \code{Exact = TRUE} all \eqn{C}{} splits are enumerated.  When
\code{Exact = FALSE}, the p-value is estimated by Monte Carlo: in each of
\code{nboots} replications the pooled sample is randomly split into groups of
sizes \eqn{m}{} and \eqn{n}{} and \eqn{H}{} is recomputed; the p-value is the
proportion of replications for which \eqn{H \geq \hat{q}}{}.  If the estimated
p-value is zero it is replaced by \eqn{1/(2 \times \code{nboots})}{}.

\bold{Asymptotic equivalence (Theorems 3.2--3.4).} Under the null hypothesis,
or under fixed or contiguous alternatives (with bounded densities and
\eqn{m/(m+n) \to \eta \in (0,1)}{}), the permutation p-value is asymptotically
equivalent to the true p-value of \eqn{\mathcal{H}_{m,n}}{}.  The permutation
test controls the type I error at the nominal level for all \eqn{k \geq 1}{},
and the power of the two tests are asymptotically equal.

\bold{Scale tuning.} The power of \eqn{\mathcal{H}_{m,n}}{} is not
scale-invariant.  The optimal scale \eqn{\sigma^{*}}{} can be estimated and
applied automatically by passing \code{scale_psi} to the unified wrapper
\code{\link{Hausdorff_test}}, which calls this function
internally on the already-scaled samples.
}

\usage{
H_test_2s_1d(x1, x2, nboots = 2000, Exact = FALSE)
}

\arguments{
\item{x1}{a numeric vector of data sample values \eqn{\{x_{1}, \ldots, x_{m}\}}{}.}
\item{x2}{a numeric vector of data sample values \eqn{\{y_{1}, \ldots, y_{n}\}}{}.}
\item{nboots}{a positive integer giving the number of Monte Carlo permutation
replications, used only when \code{Exact = FALSE}. Defaults to \code{2000}.}
\item{Exact}{a logical value indicating whether an exact permutation p-value should be
computed by enumerating all \eqn{\frac{(m+n)!}{m!n!}}{} splits of the
pooled sample.  If \eqn{\frac{(m+n)!}{m!n!} > 2.147483647 \times 10^{9}}{},
the function automatically switches to Monte Carlo and emits a message.
Defaults to \code{FALSE}.}
}

\value{
An object of class \code{"htest"} with the following components:
\describe{
\item{\code{statistic}}{the observed Hausdorff statistic
\eqn{\hat{q} = \mathcal{H}_{m,n} = H(F^{c}_{m}, G^{c}_{n})}{}, named
\code{"H"}.}
\item{\code{p.value}}{the exact or Monte Carlo permutation p-value, i.e. the proportion of
permutation statistics
\eqn{H(F^{c\dagger}_{m}(\pi_{i}), G^{c\dagger}_{n}(\pi_{i}))}{} that are
greater than or equal to \eqn{\hat{q}}{}.  If the Monte Carlo estimate is
zero it is replaced by \eqn{1/(2 \times \code{nboots})}{}.}
\item{\code{method}}{a character string indicating the method:
\code{"Two-sample Hausdorff Test (Exact)"} when all permutations are
enumerated, or \code{"Two-sample Hausdorff Test (Monte Carlo)"} otherwise.}
\item{\code{alternative}}{the character string \code{"two-sided"}.}
\item{\code{data.name}}{a character string giving the names of the two data objects.}
}
}

\references{
Dimitrina S. Dimitrova, Yun Jia, Vladimir K. Kaishev (2026a). ``On a One
Sample Goodness-of-fit Test Based on the Hausdorff Metric''. \emph{Submitted}.

Dimitrina S. Dimitrova, Yun Jia, Vladimir K. Kaishev (2026b). ``On a Two-Sample
Multivariate Goodness-of-Fit Test Based on the Hausdorff Metric''.
\emph{Submitted to the Annals of Statistics}.

Dimitrina S. Dimitrova, Vladimir K. Kaishev, Senren Tan (2020). ``Computing
the Kolmogorov-Smirnov Distribution When the Underlying CDF is Purely
Discrete, Mixed or Continuous''. \emph{Journal of Statistical Software},
\bold{95}(10): 1--42. doi:10.18637/jss.v095.i10.
}

\seealso{
\code{\link{H_stat_2s_1d_tr}},
\code{\link{Hausdorff_test}}.
}

\examples{
## Two-sample H test: both samples from N(0,1) (null is true)
set.seed(1)
x1 <- rnorm(30)
x2 <- rnorm(30)
H_test_2s_1d(x1, x2, nboots = 1000)

## Two-sample H test: samples from N(0,1) and N(0.5,1) (null is false)
set.seed(2)
x3 <- rnorm(30)
x4 <- rnorm(30, mean = 0.5)
H_test_2s_1d(x3, x4, nboots = 1000)

## Exact permutation test for small samples
set.seed(3)
a <- rnorm(8)
b <- rnorm(8)
H_test_2s_1d(a, b, Exact = TRUE)

## Scale-tuned test via the unified wrapper (recommended)
set.seed(4)
x_exp <- rexp(60, rate = 2)
y_exp <- rexp(60, rate = 3)
Hausdorff_test(x_exp, y_exp, scale_psi = c(0.99, 0.95), scale_nperms = 500)
}
