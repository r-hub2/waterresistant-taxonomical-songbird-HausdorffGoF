\name{H_stat_1s_1d}
\alias{H_stat_1s_1d}
\title{Computes the one-sample Hausdorff test statistic for a given sample and a prespecified null distribution}

\description{
Computes the Hausdorff distance \eqn{H(F^{c}, F^{c}_{n})}{} between the planar
curve \eqn{F^{c}}{} of a prespecified null cdf \eqn{F(x)}{} and the planar curve
\eqn{F^{c}_{n}}{} of the empirical cdf based on a given sample
\eqn{\{x_{1}, \ldots, x_{n}\}}{}.
}

\details{
Given a random sample \eqn{\{x_{1}, \ldots, x_{n}\}}{} with empirical cdf
\eqn{F_{n}(x)}{} and a prespecified null cdf \eqn{F(x)}{},
\code{\link{H_stat_1s_1d}} computes the Hausdorff
test statistic
\deqn{H(F^{c}, F^{c}_{n}) = \max_{l} \rho_{1}(B_{l}, E_{l}),}{}
where \eqn{\{B_{l}\}}{} are the locally farthest vertices of the planar curve
\eqn{F^{c}_{n}}{} from the null curve \eqn{F^{c}}{}, \eqn{E_{l}}{} is the
intersection of the line \eqn{\{(x,z): x + z = x_{B_{l}} + z_{B_{l}}\}}{} with
\eqn{F^{c}}{}, and
\eqn{\rho_{1}(A,B) = \max\{|x_{A}-x_{B}|, |z_{A}-z_{B}|\}}{}.
This implements Algorithm 1 (Lemmas 2 and 3) of Dimitrova, Jia, and Kaishev
(2026a).

The function first identifies the set \eqn{\mathcal{B}_{loc}}{} of locally
farthest vertices, then for each \eqn{B_{l} \in \mathcal{B}_{loc}}{} solves the
equation \eqn{x + F(x) = x_{B_{l}} + z_{B_{l}}}{} to find \eqn{E_{l}}{}, and
finally returns the maximum of \eqn{|x^{*}_{l} - x_{B_{l}}|}{}.  When
\code{pdf} is supplied the equation is solved via Newton--Raphson; otherwise
\code{uniroot} is used.
}

\usage{
H_stat_1s_1d(x, CDF, pdf = NULL, tol = 1e-10, max.init = 1000)
}

\arguments{
\item{x}{a numeric vector of data sample values \eqn{\{x_{1}, \ldots, x_{n}\}}{}.}
\item{CDF}{a vectorised \R function for the null cdf \eqn{F(x)}{}.}
\item{pdf}{(optional) a vectorised \R function for the null density \eqn{f(x)}{}.
When supplied, Newton--Raphson iterations are used to solve the equation
\eqn{x + F(x) = \lambda}{} in Algorithm 1, which is faster and more
accurate than bisection via \code{uniroot} for smooth densities.
Defaults to \code{NULL}.}
\item{tol}{a numeric value giving the tolerance for the root-finding step.
Defaults to \code{1e-10}.}
\item{max.init}{maximum number of Newton--Raphson iterations.
Used only when a density \code{pdf} is supplied. Defaults to \code{1000}.}
}

\value{
A single numeric value: the observed Hausdorff statistic
\eqn{\hat{q} = H(F^{c}, F^{c}_{n})}{}.
}

\references{
Dimitrina S. Dimitrova, Yun Jia, Vladimir K. Kaishev (2026a). ``On a One Sample
Goodness-of-fit Test Based on the Hausdorff Metric''. \emph{Submitted}.
}

\seealso{
\code{\link{H_test_1s_1d}},
\code{\link{Hausdorff_stat}}.
}

\examples{
## Compute the H statistic for a sample from Exp(1) against the Exp(1) null
set.seed(1)
x <- rexp(50, rate = 1)
H_stat <- H_stat_1s_1d(x   = x,
                       CDF = function(t) pexp(t, rate = 1),
                       pdf = function(t) dexp(t, rate = 1))
H_stat

## Compute the H statistic for a sample from N(0,1) against the N(0,1) null
set.seed(2)
y <- rnorm(100)
H_stat2 <- H_stat_1s_1d(x   = y,
                         CDF = pnorm,
                         pdf = dnorm)
H_stat2
}
