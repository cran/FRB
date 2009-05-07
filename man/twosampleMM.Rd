\name{twosampleMM}
\alias{twosampleMM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Two Sample MM-Estimates of Location and Covariance  }
\description{
 Computes two-sample MM-estimates of multivariate location and common covariance, 
using initial two-sample S-estimates.
}

\usage{
twosampleMM(X, groups, control=MMcontrol(...), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ matrix or data frame }
  \item{groups}{ vector of 1's and 2's, indicating group numbers }
  \item{control}{a list with control parameters for tuning the MM-estimate and its computing algorithm, 
  see \code{\link{MMcontrol}}().}
  \item{...}{allows for specifying control parameters directly instead of via \code{control}}
}
\details{This function is called by \code{\link{FRBhotellingMM}}
  
  The two-sample MM-estimates are defined by first computing a two-sample S-estimate of location for each sample and common covariance, 
  then fixing its scale component and re-estimating the location vectors and shape by a more efficient M-estimate (see Tatsuoka and Tyler (2000)). 
  Tukey's biweight is used for the loss functions. By default, the first loss function (in the two-sample S-estimate) is tuned 
  in order to obtain 50\% breakdown point. The default tuning of the second loss function (M-estimate) ensures 95\% efficiency at 
  the normal model. This tuning can be changed via argument \code{control} if desired.
  %When interested in location estimates, \code{shapeEff}
 %should be \code{FALSE} (the default), in which case the particular efficiency is that of the location estimates. When interest lies in the covariance or 
  %shape part, it makes sense to set \code{shapeEff=TRUE}, in which case the shape efficiency is considered instead.
      
  The computation of the two-sample S-estimate is performed by a call to \code{\link{twosampleS}}, which uses a fast-S-type
  algorithm. Its tuning parameters can be changed via the \code{control} argument.
  
  Apart from the MM-location estimates \code{Mu1} and \code{Mu2}, the function returns both the common MM-covariance \code{Sigma} and
  common MM-shape estimate \code{Gamma} (which has determinant equal to 1). 
  Additionally, the S-estimates are returned as well (their Gaussian efficiency is usually lower than the MM-estimates but they may 
  have a lower bias). 
}
\value{
A list containing:
  \item{Mu1 }{MM-estimate of first center}
  \item{Mu2 }{MM-estimate of second center}
  \item{Sigma }{MM-estimate of covariance}
  \item{Gamma }{MM-estimate of shape}
  \item{SMu1 }{S-estimate of first center}
  \item{SMu2 }{S-estimate of second center}
  \item{SSigma }{S-estimate of covariance}
  \item{SGamma }{S-estimate of shape}
  \item{scale }{S-estimate of scale (univariate)}
  \item{c0,b,c1}{tuning parameters of the loss functions (depend on control parameters \code{bdp} and \code{eff})}
}
\references{ 
\itemize{
\item K.S. Tatsuoka and D.E. Tyler (2000). The uniqueness of S and M-functionals under non-elliptical distributions.
\emph{The Annals of Statistics}, \bold{28}, 1219-1243 }
}
\author{ Ella Roelant and Gert Willems }
%\note{ ~~further notes~~ 

% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{twosampleS}}, \code{\link{FRBhotellingMM}},  \code{\link{MMboottwosample}}, \code{\link{MMcontrol}} }
\examples{
Y1 <- matrix(rnorm(50*5), ncol=5)
Y2 <- matrix(rnorm(50*5), ncol=5)
Ybig <- rbind(Y1,Y2)
grp <- c(rep(1,50),rep(2,50))
MMests <- twosampleMM(Ybig, grp)

# MM-estimate of first center:
MMests$Mu1
# MM-estimate of second center:
MMests$Mu1
# MM-estimate of common covariance:
MMests$Sigma
#initial S-estimate of first center:
MMests$SMu1
#initial S-estimate of second center:
MMests$SMu2
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
