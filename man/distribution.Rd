\name{distribution}
\alias{distribution}
\alias{print.NullDist}
\title{Specify a null distribution for the one-sample Hausdorff test}

\description{
Creates an object of class \code{"NullDist"} that bundles the null
cumulative distribution function \eqn{F(x)}{} together with its optional
quantile function \eqn{F^{-1}(p)}{} and density \eqn{f(x)}{}.
The resulting object is passed as the second argument \code{y} to
\code{\link{Hausdorff_test}} or
\code{\link{Hausdorff_stat}} to trigger the one-sample
dispatch path.
}

\usage{
distribution(CDF, CDFinverse = NULL, pdf = NULL)
}

\arguments{
\item{CDF}{an \R function for the null cdf \eqn{F(x)}{}.  This argument is
mandatory.  The function is passed through \code{Vectorize} internally, so
it does not need to operate on vectors natively.}
\item{CDFinverse}{(optional) an \R function for the null quantile function
\eqn{F^{-1}(p)}{}.  The function is passed through \code{Vectorize}
internally.  It is additionally extended with the standard
probability-space boundary convention used by
\code{\link{H_test_c_cdf}} when evaluating
the boundary values \eqn{a_{i}}{} and \eqn{b_{i}}{} (Theorem 4 of
Dimitrova, Jia, and Kaishev 2026a):
\deqn{
F^{-1}(p) =
\left\{
\begin{array}{ll}
-\infty, & \mbox{if } p < 0, \\
+\infty, & \mbox{if } p > 1, \\
F^{-1}(p), & \mbox{otherwise.}
\end{array}
\right.
}{
F^(-1)(p) = -Inf if p < 0; +Inf if p > 1; F^(-1)(p) otherwise.
}
Only in-range values of \code{p} are passed to the user-supplied function;
out-of-range values are assigned the boundary constants directly, so no
warnings are produced.  When supplied, boundary values in
\code{\link{H_test_c_cdf}} are computed
directly rather than via root-finding, giving higher numerical precision.
Defaults to \code{NULL}.}
\item{pdf}{(optional) an \R function for the null density \eqn{f(x)}{}.  The
function is passed through \code{Vectorize} internally.  When supplied,
Newton--Raphson iterations are used in
\code{\link{H_stat_1s_1d}} to solve the
equation \eqn{x + F(x) = \lambda}{}, which is faster and more accurate
than \code{uniroot} for smooth densities.  Also used to compute
\eqn{\sigma^{*}}{} numerically when \code{CDFinverse} is absent and
\code{scale_psi} is supplied to
\code{\link{Hausdorff_test}}.  Defaults to \code{NULL}.}
}

\value{
An object of class \code{"NullDist"}: a named list with elements
\code{CDF}, \code{CDFinverse}, and \code{pdf}, all stored in their
vectorised and (for \code{CDFinverse}) boundary-extended forms.  A
\code{print} method is provided.
}

\references{
Dimitrina S. Dimitrova, Yun Jia, Vladimir K. Kaishev (2026a). ``On a One Sample
Goodness-of-fit Test Based on the Hausdorff Metric''. \emph{Submitted}.
}

\seealso{
\code{\link{Hausdorff_test}},
\code{\link{Hausdorff_stat}}.
}

\examples{
## Standard normal null (all three functions supplied)
null_norm <- distribution(CDF        = pnorm,
                          CDFinverse = qnorm,
                          pdf        = dnorm)
null_norm

## Exp(2) null -- only CDF supplied (root-finding used internally)
null_exp2 <- distribution(CDF = function(x) pexp(x, rate = 2))
null_exp2
}
