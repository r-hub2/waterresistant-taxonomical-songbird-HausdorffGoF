\name{HausdorffGoF-package}
\alias{HausdorffGoF-package}
\alias{HausdorffGoF}
\title{One- and Two-Sample Hausdorff Goodness-of-Fit Test}

\description{
This package computes the test statistic and p-values of the one-sample and
two-sample Hausdorff (\eqn{H}{}) goodness-of-fit test, proposed by Dimitrova,
Jia, and Kaishev (2026a) and Dimitrova, Jia, and Kaishev (2026b) respectively.
Exploiting the scale dependence of \eqn{H}{}, it is shown that the sensitivity
of the \eqn{H}{} test can be controlled, making it left-tail, central body or
right-tail sensitive and thereby allowing its power to be optimized for the
problem at hand.

Simulation studies further demonstrate that, in terms of statistical power, the
\eqn{H}{} test substantially outperforms classical goodness-of-fit tests such as
the Kolmogorov-Smirnov, Cramer-von Mises, and Anderson-Darling tests. The
improvement is particularly pronounced when discrepancies in the right tail of
the distribution are of primary interest, as is often the case in finance and
insurance, economics, the natural sciences, and extreme value theory. These
results make the \eqn{H}{} test a compelling and effective alternative to
existing testing procedures.
 
In the \bold{one-sample} setting, given a random sample
\eqn{{X_{1}, \ldots, X_{n}}}{} of size \eqn{n}{} with empirical cdf
\eqn{F_{n}(x)}{} and a prespecified null cdf \eqn{F(x)}{}, the \eqn{H}{}
statistic is defined as the Hausdorff distance between the planar curves (completed graphs)
\eqn{F^{c}}{} and \eqn{F^{c}_{n}}{} of \eqn{F(x)}{} and \eqn{F_{n}(x)}{}
(completed by vertical segments at jump discontinuities). Based on Lemma 2.1 of Dimitrova, Jia, and Kaishev (2026a), the latter can be expressed as

\deqn{H(F^{c}, F^{c}_{n}) =  \sup_{A \in F^{c}_{n}} \inf_{B \in F^{c}} \rho_{1}(A,B)= \sup_{A \in F^{c}} \inf_{B \in F^{c}_{n}} \rho_{1}(A,B) ,}{}
 where, \eqn{\rho_1(\cdot, \cdot)}{} is the Chebyshev
(\eqn{\ell_{\infty}}{}) i.e., \eqn{\rho_{1}(A,B) = \max{|x_{A}-x_{B}|, |z_{A}-z_{B}|}}{}.

In the \bold{two-sample} setting, the null hypothesis is \eqn{H_{0}: F(x) = G(x)}{} for all \eqn{x \in R^{k}}{}, against the 
alternative \eqn{H_{1}: F(x) \neq G(x)}{} for at least one \eqn{x}{}, where \eqn{F(x)}{} and \eqn{G(x)}{} are unknown \eqn{k}{}-dimensional cdfs. 
The null hypothesis \eqn{H_{0}}{} is tested using two independent random samples
\eqn{{X_{1}, \ldots, X_{m}}}{} and \eqn{{Y_{1}, \ldots, Y_{n}}}{} of sizes
\eqn{m}{} and \eqn{n}{} from \eqn{F(x)}{} and \eqn{G(x)}{}
respectively. 

The two-sample \eqn{H}{} statistic is then defined as the Hausdorff
distance between the completed graphs \eqn{F^{c}_{m}}{} and \eqn{G^{c}_{n}}{} of
the empirical cdfs \eqn{F_{m}(x)}{} and \eqn{G_{n}(x)}{}, again under the
Chebyshev distance \eqn{\rho_{1}}{}. Similarly to the one-sample case, using Lemma 2.4 of Dimitrova, Jia, and Kaishev (2026b), it can be expressed as
\deqn{\mathcal{H}_{m,n} = H(F^{c}_{m}, G^{c}_{n}) = \sup_{A \in F^{c}_{m}} \inf_{B \in G^{c}_{n}} \rho_{1}(A,B) = \sup_{A \in G^{c}_{n}} \inf_{B \in F^{c}_{m}} \rho_{1}(A,B).}{}

A geometric interpretation given by Theorem 2.5 of Dimitrova, Jia, and Kaishev
(2026b) is that the \eqn{H}{} statistic equals the side length of the largest
square (hypercube) that can be inscribed in the region between the two curves (hypersurfaces).
In the univariate case this is equivalent to the L\'evy metric between
\eqn{F}{} and \eqn{F_{n}}{} or \eqn{F_{m}}{} and \eqn{G_{n}}{} (see Remark 6 of Dimitrova, Jia, and Kaishev 2026b).

The package \pkg{HausdorffGoF} implements efficient algorithms to compute the
\eqn{H}{} statistic in the one-sample (Algorithm 1 of Dimitrova, Jia, and
Kaishev 2026a), two-sample univariate (Section 3.1 of Dimitrova, Jia, and
Kaishev 2026b) and two-sample bivariate (Appendix A of Dimitrova, Jia, and
Kaishev 2026b) cases. Exact one-sample p-values are provided via a
rectangle-probability representation (Theorem 4 of Dimitrova, Jia, and Kaishev
2026a), evaluated in \eqn{O(n^{2}\log(n))}{} time using the Exact-KS-FFT method
of Dimitrova, Kaishev, and Tan (2020) via
\code{\link[KSgeneral]{ks_c_cdf_Rcpp}} from the package
\pkg{KSgeneral}. Two-sample p-values are obtained by applying exact or Monte Carlo
permutation of the \eqn{H}{} test (see functions \code{\link{H_test_2s_1d}} and \code{\link{H_test_2s_2d}}).

A key feature of the \eqn{H}{} test is that its p-values are \emph{not}
scale-invariant. This property is used to make the test
the left-tail or body or right-tail sensitive and correspondingly optimize the power of the test. This is achieved by rescaling the data with a scale parameter \eqn{\sigma > 0}{} that is selected as follows.

In the \bold{one-sample} case the scaled statistic is
\eqn{\mathcal{H}_{n}(\sigma) = H(F^{c}_{\sigma}, F^{c}_{n,\sigma})}{}, where
\eqn{F_{\sigma}(x) = F(x/\sigma)}{} and \eqn{F_{n,\sigma}(x) = F_{n}(x/\sigma)}{}.
An explicit scale-tuning formula
\deqn{\sigma^{*} = \frac{\psi_{1} - \psi_{2}}{F^{-1}(\psi_{1}) - F^{-1}(\psi_{2})},}{}
depending only on the null distribution \eqn{F}{} and two probability levels
\eqn{\psi_{1} > \psi_{2} \in (0,1)}{}, allows the power to be focused on the
right tail (e.g. \eqn{\psi_{1} = 0.99, \psi_{2} = 0.95}{}), the left tail
(e.g. \eqn{\psi_{1} = 0.05, \psi_{2} = 0.01}{}), or the body
(e.g. \eqn{\psi_{1} = 0.70, \psi_{2} = 0.40}{}) of the distribution. Proposition
13 of Dimitrova, Jia, and Kaishev (2026a) shows that
\eqn{\mathcal{H}_{n}(\sigma^{*})}{} and its p-value are then invariant under any
further affine rescaling of the data.

In the \bold{two-sample} case no null distribution is assumed, so \eqn{\sigma^{*}}{}
cannot be computed from \eqn{F}{} directly. Instead, the scale is estimated from
the pooled sample via the permutation-based rule
\deqn{\sigma^{*} = \mathrm{E}_{R_{m+n}}\!\left[
  \max\!\left(
    \frac{\psi_{1}-\psi_{2}}{Q_{X^{\dagger}_{m}}(\psi_{1}) - Q_{X^{\dagger}_{m}}(\psi_{2})},\;
    \frac{\psi_{1}-\psi_{2}}{Q_{Y^{\dagger}_{n}}(\psi_{1}) - Q_{Y^{\dagger}_{n}}(\psi_{2})}
  \right) \;\Bigg|\; Z_{m+n}
\right],}{}
where \eqn{Q_{X^{\dagger}_{m}}(\psi)}{} and \eqn{Q_{Y^{\dagger}_{n}}(\psi)}{}
are the \eqn{\psi}{}-quantiles of permuted sub-samples \eqn{X^{\dagger}_{m}}{}
and \eqn{Y^{\dagger}_{n}}{} of the pooled sample
\eqn{Z_{m+n} = \{x_{1},\ldots,x_{m},y_{1},\ldots,y_{n}\}}{} (Equation (49) of
Dimitrova, Jia, and Kaishev 2026b). Theorem 4.4 of Dimitrova, Jia, and Kaishev
(2026b) shows that \eqn{\sigma^{*}}{} converges in probability to the one-sample
rule applied to the pooled distribution \eqn{E = \eta F + (1-\eta)G}{}, and
that the critical value of the scale-tuned permutation test is asymptotically
equivalent to that of the original \eqn{H}{} test. Recommended choices of
\eqn{(\psi_{1}, \psi_{2})}{} are the same as in the one-sample case: right tail
\eqn{(0.99,\,0.95)}{}, body \eqn{(0.70,\,0.40)}{} or \eqn{(0.60,\,0.30)}{},
left tail \eqn{(0.05,\,0.01)}{}.

Scale tuning is integrated directly into
\code{\link{Hausdorff_test}} via the \code{scale_psi}
argument, which accepts the pair \eqn{(\psi_{1}, \psi_{2})}{} and computes
\eqn{\sigma^{*}}{} automatically before running the test.

The package provides the following exported functions.  For the
\bold{one-sample} \eqn{H}{} test:
\code{\link{H_stat_1s_1d}} computes the test
statistic, \code{\link{H_test_1s_1d}} computes
the exact p-value, and
\code{\link{H_test_c_cdf}} is the low-level
rectangle-probability engine used by \code{H_test_1s_1d} and useful when
\eqn{\hat{q}}{} is already available.
For the \bold{two-sample} \eqn{H}{} test:
\code{\link{H_test_2s_1d}} performs the
univariate permutation test,
\code{\link{H_stat_2s_1d_tr}} computes
the univariate two-sample statistic via the C++ transformation method,
\code{\link{H_stat_2s_1d_p}} computes
the univariate two-sample statistic via the C++ projection method,
\code{\link{H_test_2s_2d}} performs the bivariate permutation test, and
\code{\link{H_stat_2s_2d}} computes the
Hausdorff distance between two bivariate empirical cdfs.
The \bold{unified interface}
\code{\link{Hausdorff_test}} and
\code{\link{Hausdorff_stat}} are S3 generics that dispatch
to the appropriate underlying function based on the type of input, covering
one-sample, two-sample univariate, and two-sample bivariate cases in a single
call. Null distributions are specified via
\code{\link{distribution}}, which constructs a
\code{"NullDist"} object bundling the cdf, quantile function, and density.
}

\details{
\bold{One-sample \eqn{H}{} test}:

The algorithm to compute the \eqn{H}{} statistic is based on Lemmas 2 and 3 of
Dimitrova, Jia, and Kaishev (2026a). It identifies the locally farthest vertices
\eqn{\{B_{l}\}}{} of the planar curve \eqn{F^{c}_{n}}{} from the null curve
\eqn{F^{c}}{}, then for each \eqn{B_{l}}{} finds the intersection point
\eqn{E_{l}}{} of the line \eqn{\{(x,z) : x + z = x_{B_{l}} + z_{B_{l}}\}}{}
with \eqn{F^{c}}{}, and returns \eqn{\max_{l} \rho_{1}(B_{l}, E_{l})}{}.
This is implemented in
\code{\link{H_stat_1s_1d}}.

The exact p-value \eqn{P(H(F^{c}, F^{c}_{n}) > q)}{} is computed by
\code{\link{H_test_1s_1d}} using the
rectangle-probability representation of Theorem 4 of Dimitrova, Jia, and
Kaishev (2026a):
\deqn{P(H(F^{c}, F^{c}_{n}) > q) = 1 - P(a_{i} \leq U_{(i)} \leq b_{i},\; 1 \leq i \leq n),}{}
where \eqn{U_{(1)} \leq \cdots \leq U_{(n)}}{} are the order statistics of
\eqn{n}{} i.i.d. Uniform\eqn{(0,1)}{} random variables, and
\deqn{a_{i} = F\!\left(F^{-1}\!\!\left(\frac{i}{n} - q\right) - q\right), \quad b_{i} = F\!\left(F^{-1}\!\!\left(\frac{i-1}{n} + q\right) + q\right).}{}
This rectangle probability is evaluated in \eqn{O(n^{2}\log(n))}{} time via
\code{KSgeneral::ks_c_cdf_Rcpp}, with the boundary construction handled by
\code{\link{H_test_c_cdf}}.  A Monte Carlo
alternative (bootstrap p-value) is available by passing \code{method = "mc"}
to \code{\link{Hausdorff_test}}.

\bold{Two-sample \eqn{H}{} test}:

The two-sample \eqn{H}{} statistic
\eqn{\mathcal{H}_{m,n} = H(F^{c}_{m}, G^{c}_{n})}{} is computed by
\code{\link{H_stat_2s_1d_tr}} using the
transformation method of Section 3.1 of Dimitrova, Jia, and Kaishev (2026b).
The two staircase curves \eqn{F^{c}_{m}}{} and \eqn{G^{c}_{n}}{} are first
modified so that one lies entirely above the other (Lemma 2.13), then the plane
is rotated by \eqn{\pi/4}{} so that \eqn{\mathcal{H}_{m,n}}{} becomes the
supremum of half the vertical difference between the rotated curves (Section
3.1, ibid.), attained over a finite set of transformed vertices. The bivariate
version \code{\link{H_stat_2s_2d}} implements
the projection method of Appendix A of Dimitrova, Jia, and Kaishev (2026b) via
an internal C++ routine.

Because the distribution of \eqn{\mathcal{H}_{m,n}}{} depends on the unknown
underlying distributions \eqn{F}{} and \eqn{G}{}, p-values are obtained by the
permutation approach of Proposition 3.5 of Dimitrova, Jia, and Kaishev (2026b).
The permutation Hausdorff statistic \eqn{\mathcal{H}^{\dagger}_{m,n}}{} is
defined by randomly splitting the pooled sample \eqn{Z_{m+n}}{} into two groups
of sizes \eqn{m}{} and \eqn{n}{}, computing \eqn{H}{} on each split, and
averaging over all \eqn{C = \frac{(m+n)!}{m!n!}}{} splits:
\deqn{P(\mathcal{H}^{\dagger}_{m,n} > q \mid Z_{m+n}) = \frac{1}{C} \sum_{i=1}^{C} \mathbf{1}\!\left\{H(F^{c\dagger}_{m}(\pi_{i}), G^{c\dagger}_{n}(\pi_{i})) > q\right\}.}{}
Theorems 3.2--3.4 of Dimitrova, Jia, and Kaishev (2026b) establish that this
permutation p-value is asymptotically equivalent to the true p-value of
\eqn{\mathcal{H}_{m,n}}{} under the null and under fixed or contiguous
alternatives, and that the permutation test controls the type I error at the
nominal level for any \eqn{k \geq 1}{}. This is implemented in
\code{\link{H_test_2s_1d}}, which performs exact
enumeration when \eqn{C \leq 2.147 \times 10^{9}}{} and Monte Carlo permutation
otherwise.

The asymptotic null distribution of the normalised statistic
\eqn{\sqrt{mn/(m+n)}\,\mathcal{H}^{\dagger}_{m,n}}{} is given by Theorem 3.6
of Dimitrova, Jia, and Kaishev (2026b) and involves the supremum of a Brownian
bridge weighted by the density of the pooled distribution.

\bold{Unified interface}:

\code{\link{Hausdorff_test}} and
\code{\link{Hausdorff_stat}} are S3 generics that dispatch
on the class of their second argument \code{y}: a
\code{\link{distribution}} object triggers the one-sample path;
a numeric vector triggers the two-sample univariate path; a matrix or list
triggers the two-sample bivariate path.  \code{Hausdorff_test} also accepts
\code{scale_psi} to activate scale tuning in a single call, computing
\eqn{\sigma^{*}}{} automatically and returning the augmented \code{"htest"}
object with \code{$sigma} and \code{$scale_psi} components.
}

\author{
Dimitrina S. Dimitrova \email{D.Dimitrova@city.ac.uk},
Yun Jia \email{yunjia2019@gmail.com},
Vladimir K. Kaishev \email{Vladimir.Kaishev.1@city.ac.uk}

Maintainer: Yun Jia \email{yunjia2019@gmail.com}
}

\references{
Dimitrina S. Dimitrova, Yun Jia, Vladimir K. Kaishev (2026a). ``On a One Sample
Goodness-of-fit Test Based on the Hausdorff Metric''. \emph{Submitted}.

Dimitrina S. Dimitrova, Yun Jia, Vladimir K. Kaishev (2026b). ``On a Two-Sample
Multivariate Goodness-of-Fit Test Based on the Hausdorff Metric''.
\emph{Submitted}.

Dimitrina S. Dimitrova, Vladimir K. Kaishev, Senren Tan (2020). ``Computing the
Kolmogorov-Smirnov Distribution When the Underlying CDF is Purely Discrete,
Mixed or Continuous''. \emph{Journal of Statistical Software}, \bold{95}(10):
1--42. doi:10.18637/jss.v095.i10.
}

