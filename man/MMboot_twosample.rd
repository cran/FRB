\name{MMboot_twosample}
\alias{MMboot_twosample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fast and Robust Bootstrap for Two-Sample MM-estimates of Location and Covariance} 
\description{
  Calculates bootstrapped two sample MM-estimates using the Fast and Robust Bootstrap
method. 

}
\usage{
MMboot_twosample(X, groups, R, ests = MMest_twosample(X, groups))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ matrix of data frame }
  \item{groups}{ vector of 1's and 2's, indicating group numbers }
  \item{R}{ number of bootstrap samples }
  \item{ests}{ original MM-estimates as returned by \code{\link{MMest_twosample}}() }
}
\details{
 This function is called by \code{\link{FRBhotellingMM}}, it is typically not to be used on its own. 
It requires the result of \code{\link{MMest_twosample}} applied on \code{X}, supplied through the argument \code{ests}. 
If \code{ests} is not provided, \code{\link{MMest_twosample}} will be called with default arguments. 

The fast and robust bootstrap was first developed by Salibian-Barrera and Zamar (2002) for univariate regression MM-estimators. 

The value \code{centered} gives a matrix with \code{R} columns and \eqn{2*(2*p+p*p)} rows (\eqn{p} is the number of variables in \code{X}), 
containing the recalculated estimates of the MM-locations, MM-shape, S-covariance and S-locations. 
Each column represents a different bootstrap sample. 
The first \eqn{p} rows are the MM-location estimates of the first sample, the next \eqn{p} rows are the MM-location estimates of the second sample, 
the next \eqn{p*p} rows are the common MM-shape estimates (vectorized). Then the next 
\eqn{p*p} rows are the common S-covariance estimates (vectorized), the next \eqn{p} are the S-location estimates of the first sample, 
the final \eqn{p} rows are the S-location estimates of the second sample. 
The estimates are centered by the original estimates, which are also returned through \code{MMest} in vectorized form.

}
\value{
  A list containing:
  \item{centered}{recalculated two sample MM- and S-estimates of location and scatter (centered by original estimates), see Details}
  \item{MMest}{original two sample MM- and S-estimates of location and scatter, see Details}

}
\references{
\itemize{ 
\item M. Salibian-Barrera, S. Van Aelst and G. Willems (2008) Fast and robust 
bootstrap. \emph{Statistical Methods and Applications}, \bold{17}, 41-71. 
\item M. Salibian-Barrera, R.H. Zamar (2002) Bootstrapping robust estimates of 
regression. \emph{The Annals of Statistics}, \bold{30}, 556-582.
}
}
\author{ Ella Roelant and Gert Willems }
%\note{ ~~further notes~~ 

% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ See Also \code{\link{FRBhotellingMM}},  \code{\link{MMest_twosample}},  \code{\link{Sboot_twosample}} }
 
\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
Y1 <- matrix(rnorm(50*5), ncol=5)
Y2 <- matrix(rnorm(50*5), ncol=5)
Ybig <- rbind(Y1,Y2)
grp <- c(rep(1,50),rep(2,50))
MMests <- MMest_twosample(Ybig, grp)
bootresult <- MMboot_twosample(Ybig, grp, R=1000, ests=MMests)
}  
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
