\name{Sboottwosample}
\alias{Sboottwosample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fast and Robust Bootstrap for Two-Sample S-estimates of Location and Covariance}
\description{
  Calculates bootstrapped two-sample S-estimates using the Fast and Robust Bootstrap
method. 

}
\usage{
Sboottwosample(X, groups, R, ests = twosampleS(X, groups))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ matrix or data frame }
  \item{groups}{ vector of 1's and 2's, indicating group numbers }
 \item{R}{ number of bootstrap samples }
  \item{ests}{ original two-sample S-estimates as returned by \code{\link{twosampleS}}() }
}
\details{This function is called by \code{\link{FRBhotellingS}}, it is typically not to be used on its own. 
It requires the result of \code{\link{twosampleS}} applied on \code{X}, supplied through the argument \code{ests}. 
If \code{ests} is not provided, \code{\link{twosampleS}} will be called with default arguments. 

The fast and robust bootstrap was first developed by Salibian-Barrera and Zamar (2002) for univariate regression MM-estimators. 

The value \code{centered} gives a matrix with \code{R} columns and \eqn{2*p+p*p} rows (\eqn{p} is the number of variables in \code{X}), 
containing the recalculated estimates of the S-location for the first and second center and common S-covariance. Each column represents 
a different bootstrap sample. 
The first \eqn{p} rows are the location estimates of the first center, the next \eqn{p} rows are the location
estimates of the second center and the last \eqn{p*p} rows are the common covariance estimates (vectorized). The estimates
are centered by the original estimates, which are also returned through \code{Sest}.

 
}
\value{
  A list containing:
  \item{centered}{recalculated estimates of location of first and second center
 and covariance (centered by original estimates)}
  \item{Sest}{original estimates of first and second center and common covariance}

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
\seealso{ \code{\link{FRBhotellingS}} and \code{\link{twosampleS}} }
\examples{
Y1 <- matrix(rnorm(50*5), ncol=5)
Y2 <- matrix(rnorm(50*5), ncol=5)
Ybig <- rbind(Y1,Y2)
grp <- c(rep(1,50),rep(2,50))
Sests <- twosampleS(Ybig, grp, bdp=0.25)
bootresult <- Sboottwosample(Ybig,grp,R=1000,ests=Sests)

%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
