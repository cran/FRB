\name{Sboot_loccov}
\alias{Sboot_loccov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fast and Robust Bootstrap for S-estimates of location/covariance}
\description{
  Calculates bootstrapped S-estimates using the Fast and Robust Bootstrap
method. 
}
\usage{
Sboot_loccov(Y, R, ests = Sest_loccov(Y))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{ matrix or data frame }
  \item{R}{ number of bootstrap samples }
  \item{ests}{ original S-estimates as returned by \code{\link{Sest_loccov}}() }
}
\details{
This function is called by \code{\link{FRBpcaS}} and \code{\link{FRBhotellingS}}, it is typically not to be used on its own. 
It requires the result of \code{\link{Sest_loccov}} applied on \code{Y}, supplied through the argument \code{ests}. 
If \code{ests} is not provided, \code{\link{Sest_loccov}} will be called with default arguments. 

The fast and robust bootstrap was first developed by Salibian-Barrera and Zamar (2002) for univariate regression MM-estimators. 

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
\item M. Salibian-Barrera, S. Van Aelst and G. Willems (2008) Fast and robust 
bootstrap. \emph{Statistical Methods and Applications}, \bold{17}, 41-71. 
\item M. Salibian-Barrera, R.H. Zamar (2002) Bootstrapping robust estimates of 
regression. \emph{The Annals of Statistics}, \bold{30}, 556-582.
}
}
\author{ Gert Willems and Ella Roelant }
%\note{ ~~further notes~~ 

% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{FRBpcaS}}, \code{\link{FRBhotellingS}},  \code{\link{Sest_loccov}},  \code{\link{MMboot_loccov}} }
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
