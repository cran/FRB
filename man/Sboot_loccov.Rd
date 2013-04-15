\name{Sboot_loccov}
\alias{Sboot_loccov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fast and Robust Bootstrap for S-estimates of location/covariance}
\description{
  Calculates bootstrapped S-estimates using the Fast and Robust Bootstrap
method. 
}
\usage{
Sboot_loccov(Y, R = 999, ests = Sest_loccov(Y))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{ matrix or data frame. }
  \item{R}{ number of bootstrap samples. Default is \code{R=999}. }
  \item{ests}{ original S-estimates as returned by \code{\link{Sest_loccov}}(). }
}
\details{
This function is called by \code{\link{FRBpcaS}} and \code{\link{FRBhotellingS}}, it is typically not to be used on its own. 
It requires the S-estimates of multivariate location and scatter/shape
(the result of \code{\link{Sest_loccov}} applied on \code{Y}), supplied through the argument \code{ests}. 
If \code{ests} is not provided, \code{\link{Sest_loccov}} calls the implementation of the multivariate S-estimates in package \pkg{rrcov} of Todorov and Filzmoser (2009) with default arguments.

For multivariate data the fast and robust bootstrap was developed by Salibian-Barrera, Van Aelst and Willems (2006).

The value \code{centered} gives a matrix with \code{R} columns and \eqn{p+p*p} rows (\eqn{p} is the number of variables in \code{Y}), 
containing the recalculated estimates of the S-location and -covariance. Each column represents a different bootstrap sample. 
The first \eqn{p} rows are the location estimates and the next \eqn{p*p} rows are the covariance estimates (vectorized). The estimates
are centered by the original estimates, which are also returned through \code{Sest}.
}
\value{
A list containing:
  \item{centered}{recalculated estimates of location and covariance (centered by original estimates)}
  \item{Sest}{original estimates of location and covariance}
}
\references{ 
\itemize{
\item M. Salibian-Barrera, S. Van Aelst and G. Willems (2006) PCA based on multivariate MM-estimators
with fast and robust bootstrap. \emph{Journal of the American Statistical Association}, \bold{101}, 1198--1211.
\item M. Salibian-Barrera, S. Van Aelst and G. Willems (2008) Fast and robust 
bootstrap. \emph{Statistical Methods and Applications}, \bold{17}, 41--71. 
\item V. Todorov and P. Filzmoser (2009) An object-oriented framework for robust multivariate analysis. \emph{Journal of Statistical Software}, \bold{32}, 1--47. URL \url{http://www.jstatsoft.org/v32/i03/}.
\item S. Van Aelst and G. Willems (2013). Fast and robust bootstrap for multivariate inference: The R package FRB. \emph{Journal of Statistical Software}, \bold{53}(3), 1--32.
        URL: \url{http://www.jstatsoft.org/v53/i03/}.
        }
}
\author{ Gert Willems, Ella Roelant and Stefan Van Aelst}
%\note{ ~~further notes~~ 

% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{FRBpcaS}}, \code{\link{FRBhotellingS}},    \code{\link{MMboot_loccov}} }
\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.

Y <- matrix(rnorm(50*5), ncol=5)
Sests <- Sest_loccov(Y, bdp = 0.25) 
bootresult <- Sboot_loccov(Y, R = 1000, ests = Sests)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
