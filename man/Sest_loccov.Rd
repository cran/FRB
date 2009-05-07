\name{Sest_loccov}
\alias{Sest_loccov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ S-estimates of location/covariance}
\description{
  Computes S-estimates of multivariate location and covariance using the fast-S algorithm
}
\usage{
Sest_loccov(Y, bdp = 0.5, control=Scontrol(...), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{ matrix or data frame }
  \item{bdp}{ required breakdown point of the S-estimate. Should have \eqn{0 < } \code{bdp} \eqn{\le 0.5}, the default is 0.5.}
  \item{control}{a list with control parameters for tuning the computing algorithm, see \code{\link{Scontrol}}().}
  \item{...}{allows for specifying control parameters directly instead of via \code{control}}
}
\details{
  This function is called by \code{\link{FRBpcaS}} and \code{\link{FRBhotellingS}}. 
  
  Multivariate S-estimates were introduced by Davies (1987). The algorithm used here is a multivariate version of the
  fast-S algorithm introduced by Salibian-Barrera and Yohai (2006). 
  See \code{\link{Scontrol}} for the adjustable tuning parameters of this algorithm.
  
  The function both returns the covariance estimate \code{Sigma} and shape estimate \code{Gamma} (which has determinant equal to 1). 
  The \code{scale} is determined by \eqn{det(Sigma)^{1/2/p}}, with \eqn{p} the number of variables.
}
\value{
A list containing:
  \item{Mu}{S-estimate of location}
  \item{Gamma}{S-estimate of shape}
  \item{Sigma}{S-estimate of covariance}
  \item{scale}{S-estimate of scale (univariate)}
  \item{b,c}{tuning parameters used in Tukey biweight loss function, as determined by \code{bdp}}
}
\references{ 
\itemize{
\item P.L. Davies (1987). Asymptotic behavior of S-estimates of multivariate location parameters and dispersion
matrices. \emph{The Annals of Statistics}, \bold{15}, 1269-1292.
\item M. Salibian-Barrera and V. Yohai (2006) A fast algorithm for S-regression estimates. 
\emph{Journal of Computational and Graphical Statistics}, \bold{15}, 414-427. 
}
}
\author{ Gert Willems and Ella Roelant }
%\note{ ~~further notes~~ 

% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{FRBpcaS}}, \code{\link{FRBhotellingS}}, \code{\link{Sboot_loccov}}, \code{\link{MMest_loccov}}, \code{\link{Scontrol}}  }
\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
Y <- matrix(rnorm(50*5), ncol=5)
Sests <- Sest_loccov(Y, bdp = 0.25) 
# S-estimate of location:
Sests$Mu
# S-estimate of covariance:
Sests$Sigma

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
