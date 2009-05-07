\name{MMboot_loccov}
\alias{MMboot_loccov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fast and Robust Bootstrap for MM-estimates of Location and Covariance }
\description{
  Calculates bootstrapped MM-estimates of multivariate location and scatter using the Fast and Robust Bootstrap
method. }
\usage{
MMboot_loccov(Y, R, ests = MMest_loccov(Y))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{ matrix or data frame }
  \item{R}{ number of bootstrap samples }
  \item{ests}{ original MM-estimates as returned by \code{\link{MMest_loccov}}() }
}
\details{
This function is called by \code{\link{FRBpcaMM}} and \code{\link{FRBhotellingMM}}, it is typically not to be used on its own. 
It requires the result of \code{\link{MMest_loccov}} applied on \code{Y}, supplied through the argument \code{ests}. 
If \code{ests} is not provided, \code{\link{MMest_loccov}} will be called with default arguments. 

The fast and robust bootstrap was first developed by Salibian-Barrera and Zamar (2002) for univariate regression MM-estimators. 

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
\seealso{ \code{\link{FRBpcaMM}}, \code{\link{FRBhotellingMM}},  \code{\link{MMest_loccov}},  \code{\link{Sboot_loccov}} }
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
