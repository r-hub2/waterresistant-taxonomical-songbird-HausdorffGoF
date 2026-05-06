\name{H_stat_2s_2d}
\alias{H_stat_2s_2d}
\title{Computes the Hausdorff distance between two bivariate empirical cdfs}

\description{
Computes the Hausdorff distance between the bivariate empirical cdfs of two
independent two-dimensional data samples
\eqn{\bm{X} = \{(X_{1,1}, X_{1,2}), \ldots, (X_{m,1}, X_{m,2})\}}{} and
\eqn{\bm{Y} = \{(Y_{1,1}, Y_{1,2}), \ldots, (Y_{n,1}, Y_{n,2})\}}{}.
}

\details{
\code{\link{H_stat_2s_2d}} computes the Hausdorff
distance \eqn{H(F^{c}_{m}, G^{c}_{n})}{} between the three-dimensional
staircase surfaces corresponding to the bivariate empirical cdfs
\eqn{F_{m}(\bm{x})}{} and \eqn{G_{n}(\bm{x})}{} of \code{x} and \code{y},
under the Chebyshev distance
\eqn{\rho_{1}(A,B) = \max_{1 \leq i \leq 3}\{|w^{(i)}_{A} - w^{(i)}_{B}|\}}{}.
It implements the projection method of Appendix A of Dimitrova, Jia, and
Kaishev (2026b), which generalises the univariate algorithm of Section 3.1 to
\eqn{k = 2}{}.
 
\bold{Bivariate empirical cdf computation.}
The evaluation of the bivariate empirical cdfs \eqn{F_{m}}{} and \eqn{G_{n}}{}
at the required grid points is performed using the fast divide-and-conquer
algorithm of Langren\'e and Warin (2021).  For \eqn{N}{} data points in
\eqn{d}{} dimensions this algorithm achieves
\eqn{O\!\left(N\log(N)^{(d-1) \vee 1}\right)}{} complexity, compared to the
naive \eqn{O(N^{2})}{}, by recursively splitting the dataset along alternating
coordinate axes and accumulating counts in each sub-problem.  In the bivariate
case (\eqn{d = 2}{}) this gives \eqn{O(N \log^2 N)}{} complexity.
 
\bold{Projection method (Appendix A, Dimitrova, Jia, and Kaishev 2026b).}
The algorithm proceeds in five steps.
 
\bold{Step 1 (omnidirectional jump vertices of \eqn{F_{m}}{}).}
A point \eqn{\alpha \in A_{x^{(1)}} \times A_{x^{(2)}}}{} is an
\emph{omnidirectional jump} of \eqn{F_{m}}{} if and only if
\eqn{F_{m}(\alpha) \neq F_{m}(\alpha - \varepsilon_{0} e_{i})}{} for both
\eqn{i = 1}{} and \eqn{i = 2}{} (Definition 2.11), where
\eqn{e_{1} = (1,0)}{}, \eqn{e_{2} = (0,1)}{}, and
\eqn{\varepsilon_{0} = \frac{1}{2}\min_{i; j \neq l} |x^{(i)}_{j} -
x^{(i)}_{l}|}{}. Each omnidirectional jump \eqn{\alpha_{i}}{} gives rise to
two vertices of the planar curve \eqn{F^{c}_{m}}{} (Equation (19)):
\deqn{A_{2i-1} = (\alpha_{i},\, F_{m}(\alpha_{i} - \varepsilon)),
\qquad
A_{2i} = (\alpha_{i},\, F_{m}(\alpha_{i})),
\qquad i = 1, \ldots, \upsilon,}{}
where \eqn{\varepsilon = (\varepsilon_{0}, \varepsilon_{0})}{}.
 
\bold{Step 2 (truncation boundary and additional vertices).}
The region between the bivariate planar curves is unbounded.  Lemma 2.3 is
applied to truncate both curves at
\eqn{M = 1 + \max\{x^{(i)}_{j}, y^{(i)}_{l}\}}{}, ensuring that
\eqn{H(\hat{F}^{c}_{m}(M), \hat{G}^{c}_{n}(M)) = H(F^{c}_{m}, G^{c}_{n})}{}.
 
\bold{Step 3 (projection of \eqn{F^{c}_{m}}{} vertices).}
Each vertex \eqn{A_{l}}{} of \eqn{\hat{F}^{c}_{m}(M)}{} is projected onto the
\eqn{x^{(1)}Ox^{(2)}}{} plane along the direction
\eqn{E_{0} = (1, 1, -1) \in {R}^{3}}{}.
 
\bold{Step 4 (projection partition for \eqn{G^{c}_{n}}{}).}
The vertices and sides of \eqn{G^{c}_{n}}{} are projected in the same
direction \eqn{E_{0}}{}, forming a partition of \eqn{{R}^{2}}{} whose regions
determine which coordinate of \eqn{A_{l}}{} drives the distance to
\eqn{G^{c}_{n}}{}.
 
\bold{Step 5 (nearest-point distances via Lemma A.3).}
For each locally farthest vertex \eqn{A_{l} \in \mathcal{V}_{loc}}{}, the
nearest point on \eqn{G^{c}_{n}}{} is found via the Pareto frontier of
dominated projected vertices, and the closed-form distance formula of Lemma
A.3 is applied.  The Hausdorff distance is
\deqn{H(F^{c}_{m}, G^{c}_{n}) = \max_{A_{l} \in \mathcal{V}_{loc}}
  \inf_{B \in G^{c}_{n}} \rho_{1}(A_{l}, B).}{}
The combinatorial search is executed by the internal C++ routine
\code{hsearch_Rcpp} (via \pkg{Rcpp}).
 
\bold{Component-wise orders and \eqn{H^{\circ}_{m,n}}{}.}
The function computes the statistic under the standard component-wise order
\eqn{\preceq_{1}}{} (\eqn{\bm{x} \preceq_{1} \bm{y}}{} iff
\eqn{x^{(1)} \leq y^{(1)}}{} and \eqn{x^{(2)} \leq y^{(2)}}{}).  The
order-invariant statistic
\eqn{H^{\circ}_{m,n} = \max_{1 \leq i \leq 4} H(F^{c}_{m,i}, G^{c}_{n,i})}{}
(Equation (17) of Dimitrova, Jia, and Kaishev 2026b) is available directly
via \code{\link{Hausdorff_stat}} with
\code{invariant = TRUE}, which calls this function over the four sign-flip
orderings and takes the maximum.
 
The inputs \code{x} and \code{y} must each be two-column numeric matrices;
otherwise the function stops with an error.  To pass list input, use
\code{\link{Hausdorff_stat}}, which converts lists to matrices before
calling this function.
}

\usage{
H_stat_2s_2d(x, y, tol = 1e-6)
}

\arguments{
\item{x}{a two-column numeric matrix representing the first bivariate sample
of size \eqn{m}{}.  For list input use \code{\link{Hausdorff_stat}}.}
\item{y}{a two-column numeric matrix representing the second bivariate sample
of size \eqn{n}{}.  For list input use \code{\link{Hausdorff_stat}}.}
\item{tol}{a numeric value giving the tolerance for locating omnidirectional jump
vertices in the bivariate empirical cdfs (see Definition 2.11 of
Dimitrova, Jia, and Kaishev 2026b). Defaults to \code{1e-6}.}
}

\value{
A single numeric value: the Hausdorff distance \eqn{H(F^{c}_{m}, G^{c}_{n})}{}
between the bivariate empirical cdfs of \code{x} and \code{y}.
}

\references{
Dimitrina S. Dimitrova, Yun Jia, Vladimir K. Kaishev (2026b). ``On a Two-Sample
Multivariate Goodness-of-Fit Test Based on the Hausdorff Metric''.
\emph{Submitted to the Annals of Statistics}.
 
Nicolas Langren\'e, Xavier Warin (2021). ``Fast Multivariate Empirical
Cumulative Distribution Function with Connection to Kernel Density Estimation''.
\emph{Computational Statistics & Data Analysis}, \bold{162}, 107267.
doi:10.1016/j.csda.2021.107267.
}

\seealso{
\code{\link{Hausdorff_stat}},
\code{\link{Hausdorff_test}}.
}

\examples{
## Hausdorff distance between two bivariate samples from the same distribution
set.seed(1)
x <- matrix(rnorm(100), ncol = 2)
y <- matrix(rnorm(100), ncol = 2)
H_stat_2s_2d(x, y)
 
## Hausdorff distance between two bivariate samples from different distributions
set.seed(2)
x2 <- matrix(rnorm(100), ncol = 2)
y2 <- matrix(rnorm(100, mean = 0.5), ncol = 2)
H_stat_2s_2d(x2, y2)
 
## List input: use Hausdorff_stat() which handles the conversion
set.seed(3)
x3 <- list(rnorm(30), rnorm(30))
y3 <- list(rnorm(30), rnorm(30))
Hausdorff_stat(x3, y3)
 
## Order-invariant statistic via the unified wrapper
Hausdorff_stat(x, y, invariant = TRUE)
}
