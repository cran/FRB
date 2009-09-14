\name{Sest_multireg}
\alias{Sest_multireg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{S-Estimates for Multivariate Regression}
\description{
  Computes S-Estimates of multivariate regression based on Tukey's biweight 
function using the fast-S algorithm. 
}
\usage{
Sest_multireg(X, Y, bdp = 0.5, control=Scontrol(...), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ a matrix or data frame containing the explanatory variables (possibly including intercept). }
  \item{Y}{ a matrix or data frame containing the response variables. }
  \item{bdp}{ required breakdown point. Should have \eqn{0 < } \code{bdp} \eqn{\le 0.5}, the default is 0.5.}
  \item{control}{a list with control parameters for tuning the computing algorithm, see \code{\link{Scontrol}}().}
  \item{...}{allows for specifying control parameters directly instead of via \code{control}}
}
\details{
  This function is called by \code{\link{FRBmultiregS}}. 
  
  S-estimates for multivariate regression were discussed in Van Aelst and Willems (2005). The algorithm used here is a multivariate 
  version of the fast-S algorithm introduced by Salibian-Barrera and Yohai (2006). 
  See \code{\link{Scontrol}} for the adjustable tuning parameters of this algorithm. 
  
  Apart from the regression coefficients \code{Beta}, the function both returns the error covariance matrix estimate \code{Sigma} and 
  the corresponding shape estimate \code{Gamma} (which has determinant equal to 1). 
  The \code{scale} is determined by \eqn{det(Sigma)^{1/2/q}}, with \eqn{q} the number of response variables.
}
\value{
  \item{Beta}{S-estimate of the regression coefficient matrix}
  \item{Gamma}{S-estimate of the error shape matrix}
  \item{Sigma}{S-estimate of the error covariance matrix}
  \item{scale}{S-estimate of the error scale (univariate)}
  \item{b,c}{tuning parameters used in Tukey biweight loss function, as determined by \code{bdp}}
  \item{w}{implicit weights corresponding to the S-estimates (i.e. final weights in the RWLS procedure at the end of the fast-S algorithm)}
  \item{outFlag}{outlier flags: 1 if the robust distance of the residual exceeds the .975 quantile of (the square root of)
  the chi-square distribution with degrees of freedom equal to the dimension of the responses; 0 otherwise}
}
\references{ 
\itemize{
\item M. Salibian-Barrera and V. Yohai (2006) A fast algorithm for S-regression estimates. 
\emph{Journal of Computational and Graphical Statistics}, \bold{15}, 414-427. 
\item  S. Van Aelst and G. Willems (2005). Multivariate regression S-estimators for robust estimation and
inference. \emph{Statistica Sinica}, \bold{15}, 981-1001
}
}
\author{ Gert Willems and Ella Roelant }
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{FRBmultiregS}}, \code{\link{Sboot_multireg}},  \code{\link{MMest_multireg}}, \code{\link{Scontrol}} }
\examples{
data(schooldata)
school.x <- data.matrix(schooldata[,1:5])
school.y <- data.matrix(schooldata[,6:8])

# compute 25\% breakdown S-estimates
Sres <- Sest_multireg(school.x,school.y, bdp=0.25)
# the regression coefficients:
Sres$Beta

n <- nrow(school.x)
par(mfrow=c(2,1))
# the estimates can be considered as weighted least squares estimates with the 
# following implicit weights
plot(1:n, Sres$w)
# Sres$outFlag tells which points are outliers based on whether or not their 
# robust distance exceeds the .975 chi-square cut-off:
plot(1:n, Sres$outFlag)
# (see also the diagnostic plot in plotDiag())

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
