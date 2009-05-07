\name{MMest_multireg}
\alias{MMest_multireg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{MM-Estimates for Multivariate Regression}
\description{
  Computes MM-Estimates of multivariate regression, using initial S-estimates 
}
\usage{
MMest_multireg(X, Y, control=MMcontrol(...), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ a matrix or data frame containing the explanatory variables (possibly including intercept). }
  \item{Y}{ a matrix or data frame containing the response variables. }
  \item{control}{a list with control parameters for tuning the MM-estimate and its computing algorithm, 
  see \code{\link{MMcontrol}}().}
  \item{...}{allows for specifying control parameters directly instead of via \code{control}}
}
\details{
  This function is called by \code{\link{FRBmultiregMM}}.
  
  The MM-estimates are defined by first computing S-estimates of regression, then fixing the scale component of the error covariance
  estimate, and finally re-estimating the regression coefficients and the shape part of the error covariance by more efficient 
  M-estimates (see Tatsuoka and Tyler (2000) for MM-estimates in the special case of location/scatter estimation, and Van Aelst and 
  Willems (2005) for S-estimates of multivariate regression). Tukey's biweight is used for 
  the loss functions. By default, the first loss function (in the S-estimates) is tuned in order to obtain 50\% breakdown point. 
  The default tuning of the second loss function (M-estimates) ensures 95\% efficiency at the normal model for the coefficient estimates. 
  This tuning can be changed via argument \code{control} if desired.
  
  The computation of the S-estimates is performed by a call to \code{\link{Sest_multireg}}, which uses the fast-S algorithm. 
  See \code{\link{MMcontrol}}() to see or change the tuning parameters for this algorithm.
  
  Apart from the MM-estimate of the regression coefficients \code{Beta}, the function returns both the MM-estimate of the error 
  covariance \code{Sigma} and the corresponding shape estimate \code{Gamma} (which has determinant equal to 1). 
  Additionally, the initial S-estimates are returned as well (their Gaussian efficiency is usually lower than the MM-estimates but they may 
  have a lower bias). 
}
\value{
A list containing:
  \item{Beta }{MM-estimate of the regression coefficient matrix}
  \item{Sigma }{MM-estimate of the error covariance matrix}
  \item{Gamma }{MM-estimate of the error shape matrix}
  \item{SBeta }{S-estimate of the regression coefficient matrix}
  \item{SSigma}{S-estimate of the error covariance matrix}
  \item{SGamma}{S-estimate of the error shape matrix}
  \item{scale }{S-estimate of scale (univariate)}
  \item{c0,b,c1}{tuning parameters of the loss functions (depend on control parameters \code{bdp} and \code{eff})}
}
\references{ 
\itemize{
\item K.S. Tatsuoka and D.E. Tyler (2000). The uniqueness of S and M-functionals under non-elliptical distributions.
\emph{The Annals of Statistics}, \bold{28}, 1219-1243. 
\item  S. Van Aelst and G. Willems (2005). Multivariate regression S-estimators for robust estimation and
inference. \emph{Statistica Sinica}, \bold{15}, 981-1001.
}
}
\author{ Gert Willems and Ella Roelant }

%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{  \code{\link{FRBmultiregMM}}, \code{\link{MMboot_multireg}},  \code{\link{Sest_multireg}}, \code{\link{MMcontrol}} }
\examples{
data(schooldata)
school.x <- data.matrix(schooldata[,1:5])
school.y <- data.matrix(schooldata[,6:8])

# compute 25\% breakdown S-estimates
MMres <- MMest_multireg(school.x,school.y)
# the MM-estimate of the regression coefficient matrix:
MMres$Beta
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
