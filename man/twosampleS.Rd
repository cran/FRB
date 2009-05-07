\name{twosampleS}
\alias{twosampleS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Two Sample S-Estimates of Location and Covariance }
\description{
  Computes two sample S-estimates of location and 
common covariance
}
\usage{
twosampleS(X, groups, bdp = 0.5, control=Scontrol(...), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ matrix or data frame }
  \item{groups}{ vector of 1's and 2's, indicating group numbers }
  \item{bdp}{ required breakdown point of the two sample S-estimate. Should have \eqn{0 < } \code{bdp} \eqn{\le 0.5}, the default is 0.5}
  \item{control}{a list with control parameters for tuning the computing algorithm, see \code{\link{Scontrol}}().}
  \item{...}{allows for specifying control parameters directly instead of via \code{control}}
}
\details{
 This function is called by \code{\link{FRBhotellingS}}. 
The algorithm is a multivariate version of the fast-S algorithm introduced by Salibian-Barrera and Yohai (2006). 
See \code{\link{Scontrol}} for the adjustable tuning parameters of this algorithm.
 
The function both returns the covariance estimate \code{Sigma} and shape estimate \code{Gamma} (which has determinant equal to 1). 
The \code{scale} is determined by \eqn{det(Sigma)^{1/2/p}}, with \eqn{p} the number of variables.

}
\value{
A list containing:
  \item{Mu1 }{S-estimate of first center}
  \item{Mu2 }{S-estimate of second center}
  \item{Sigma }{S-estimate of commmon covariance}
  \item{Gamma }{S-estimate of common shape}
  \item{scale }{S-estimate of scale (univariate)}
  \item{b,c}{tuning parameters used in Tukey biweight loss function, as determined by \code{bdp}}  
}
\references{ 
\itemize{
\item  X. He and W.K. Fung (2000) High breakdown estimation for multiple populations with 
applications to discriminant analysis. \emph{Journal of Multivariate Analysis}, \bold{72}, 151-162.
\item M. Salibian-Barrera and V. Yohai (2006) A fast algorithm for S-regression estimates. 
\emph{Journal of Computational and Graphical Statistics}, \bold{15}, 414-427. 
}
}
\author{ Ella Roelant and Gert Willems }
%\note{ ~~further notes~~ 

% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{twosampleMM}}, \code{\link{FRBhotellingS}},  \code{\link{Sboottwosample}}, \code{\link{Scontrol}} }
\examples{
Y1 <- matrix(rnorm(50*5), ncol=5)
Y2 <- matrix(rnorm(50*5), ncol=5)
Ybig <- rbind(Y1,Y2)
grp <- c(rep(1,50),rep(2,50))
Sests <- twosampleS(Ybig, grp, bdp=0.25)
  
# S-estimate of first center:
Sests$Mu1
# S-estimate of second center:
Sests$Mu1
# S-estimate of common covariance:
Sests$Sigma


%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
