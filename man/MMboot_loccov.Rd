\name{MMboot_loccov}
\alias{MMboot_loccov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fast and Robust Bootstrap for MM-estimates of Location and Covariance }
\description{
  Calculates bootstrapped MM-estimates of multivariate location and scatter using the Fast and Robust Bootstrap
method. }
\usage{
MMboot_loccov(Y, R = 999, ests = MMest_loccov(Y))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{ matrix or data frame. }
  \item{R}{ number of bootstrap samples. Default is \code{R=999}. }
  \item{ests}{ original MM-estimates as returned by \code{\link{MMest_loccov}}(). }
}
\details{
This function is called by \code{\link{FRBpcaMM}} and \code{\link{FRBhotellingMM}}, it is typically not to be used on its own. 
It requires the MM-estimates of multivariate location and scatter/shape
(the result of \code{\link{MMest_loccov}} applied on \code{Y}), supplied through the argument \code{ests}. 
If \code{ests} is not provided, \code{\link{MMest_loccov}} calls the implementation of the multivariate MM-estimates in package \pkg{rrcov} of Todorov and Filzmoser (2009) with default arguments. 

For multivariate data the fast and robust bootstrap was developed by Salibian-Barrera, Van Aelst and Willems (2006).

The value \code{centered} gives a matrix with \code{R} columns and \eqn{2*(p+p*p)} rows (\eqn{p} is the number of variables in \code{Y}), 
containing the recalculated estimates of the MM-location, MM-shape, S-covariance and S-location. 
Each column represents a different bootstrap sample. 
The first \eqn{p} rows are the MM-location estimates, the next \eqn{p*p} rows are the MM-shape estimates (vectorized). Then the next 
\eqn{p*p} rows are the S-covariance estimates (vectorized) and the final \eqn{p} rows are the S-location estimates. 
The estimates are centered by the original estimates, which are also returned through \code{MMest} in vectorized form.
}
\value{
A list containing:
  \item{centered}{recalculated MM- and S-estimates of location and scatter (centered by original estimates), see Details}
  \item{MMest}{original MM- and S-estimates of location and scatter, see Details}
}
\references{
\itemize{ 
\item M. Salibian-Barrera, S. Van Aelst and G. Willems (2006) PCA based on multivariate MM-estimators
with fast and robust bootstrap. \emph{Journal of the American Statistical Association}, \bold{101}, 1198--1211.
\item M. Salibian-Barrera, S. Van Aelst and G. Willems (2008) Fast and robust 
bootstrap. \emph{Statistical Methods and Applications}, \bold{17}, 41--71. 
regression. 
\item V. Todorov and P. Filzmoser (2009) An Object-Oriented Framework for Robust Multivariate Analysis. \emph{Journal of Statistical Software}, \bold{32}, 1--47. URL \url{http://www.jstatsoft.org/v32/i03/}.
}
}
\author{ Gert Willems, Ella Roelant and Stefan Van Aelst}
%\note{ ~~further notes~~ 

% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{FRBpcaMM}}, \code{\link{FRBhotellingMM}},  \code{\link{Sboot_loccov}} }
\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.

Y <- matrix(rnorm(50*5), ncol=5)
MMests <- MMest_loccov(Y) 
bootresult <- MMboot_loccov(Y, R = 1000, ests = MMests)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
