\name{MMest_multireg}
\alias{MMest_multireg}
\alias{MMest_multireg.default}
\alias{MMest_multireg.formula}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{MM-Estimates for Multivariate Regression}
\description{
  Computes MM-Estimates of multivariate regression, using initial S-estimates 
}
\usage{
\method{MMest_multireg}{formula}(formula, data=NULL, ...)

\method{MMest_multireg}{default}(X, Y, int = TRUE, control=MMcontrol(...),
na.action=na.omit, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{ an object of class \code{\link{formula}}; a symbolic description of the model to be fit.}
  \item{data}{ data frame from which variables specified in formula are to be taken.}
  \item{X}{ a matrix or data frame containing the explanatory variables (possibly including intercept). }
  \item{Y}{ a matrix or data frame containing the response variables. }
  \item{int}{ logical: if \code{TRUE} an intercept term is added to the model (unless it is already present in \code{X})}  
  \item{control}{a list with control parameters for tuning the MM-estimate and its computing algorithm, 
  see \code{\link{MMcontrol}}().}
  \item{na.action}{ a function which indicates what should happen when the data contain NAs. Defaults to \code{\link{na.omit}}.}
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
  The desired efficiency can be changed via argument \code{control}.
  
  The computation of the S-estimates is performed by a call to \code{\link{Sest_multireg}}, which uses the fast-S algorithm. 
  See \code{\link{MMcontrol}}() to see or change the tuning parameters for this algorithm. The M-estimate part is computed
  through iteratively reweighted least squares (RWLS).
  
  Apart from the MM-estimate of the regression coefficients, the function returns both the MM-estimate of the error 
  covariance \code{Sigma} and the corresponding shape estimate \code{Gamma} (which has determinant equal to 1). 
  Additionally, the initial S-estimates are returned as well (their Gaussian efficiency is usually lower than the MM-estimates but they may 
  have a lower bias). 

The returned object inherits from class \code{mlm} such that the standard \code{\link{coef}}, \code{\link{residuals}}, \code{\link{fitted}} and \code{\link[=predict.lm]{predict}} functions can be used.
}
\value{
  An object of class \code{FRBmultireg} which extends class \code{mlm} and contains at least the following components:
  \item{coefficients}{ MM-estimates of the regression coefficients}
%  \item{intercept}{ estimate for the intercept}
  \item{residuals}{ the residuals, that is response minus fitted values}
  \item{fitted.values}{ the fitted values.}
%  \item{Beta }{MM-estimate of the regression coefficient matrix}
  \item{Sigma }{MM-estimate of the error covariance matrix}
  \item{Gamma }{MM-estimate of the error shape matrix}
  \item{scale}{ S-estimate of the size of the multivariate errors}
  \item{weights}{implicit weights corresponding to the MM-estimates (i.e. final weights in the RWLS procedure)}
  \item{outFlag}{outlier flags: 1 if the robust distance of the residual exceeds the .975 quantile of (the square root of)
  the chi-square distribution with degrees of freedom equal to the dimension of the responses; 0 otherwise}  
  \item{c0,b,c1}{tuning parameters of the loss functions (depend on control parameters \code{bdp} and \code{eff})}
  \item{method}{ a list with following components: \code{est} = character string indicating that GS-estimates were used,\code{bdp} = a copy of the \code{bdp} argument, \code{eff} a copy of the \code{eff}  argument }
  \item{control}{ a copy of the \code{control} argument }
  \item{SBeta }{S-estimate of the regression coefficient matrix}
  \item{SSigma}{S-estimate of the error covariance matrix}
  \item{SGamma}{S-estimate of the error shape matrix}
}
\references{ 
\itemize{
\item K.S. Tatsuoka and D.E. Tyler (2000) The uniqueness of S and M-functionals under non-elliptical distributions.
\emph{The Annals of Statistics}, \bold{28}, 1219--1243. 
\item  S. Van Aelst and G. Willems (2005) Multivariate regression S-estimators for robust estimation and
inference. \emph{Statistica Sinica}, \bold{15}, 981--1001.
}
}
\author{ Gert Willems, Stefan Van Aelst and Ella Roelant }

%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{  \code{\link{FRBmultiregMM}}, \code{\link{MMboot_multireg}},  \code{\link{Sest_multireg}}, \code{\link{MMcontrol}} }
\examples{
data(schooldata)
school.x <- data.matrix(schooldata[,1:5])
school.y <- data.matrix(schooldata[,6:8])

# compute 95\% efficient MM-estimates
MMres <- MMest_multireg(school.x,school.y)
# or using the formula interface
MMres <- MMest_multireg(cbind(reading,mathematics,selfesteem)~., data=schooldata)


# the MM-estimate of the regression coefficient matrix:
MMres$coefficients
# or alternatively 
coef(MMres)

n <- nrow(schooldata)
par(mfrow=c(2,1))
# the estimates can be considered as weighted least squares estimates with the 
# following implicit weights
plot(1:n, MMres$weights)
# Sres$outFlag tells which points are outliers based on whether or not their 
# robust distance exceeds the .975 chi-square cut-off:
plot(1:n, MMres$outFlag)
# (see also the diagnostic plot in plotDiag())
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
