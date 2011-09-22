\name{GSest_multireg}
\alias{GSest_multireg}
\alias{GSest_multireg.default}
\alias{GSest_multireg.formula}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{GS Estimates for Multivariate Regression}
\description{
  Computes GS-Estimates of multivariate regression based on Tukey's biweight 
function. 
}
\usage{
\method{GSest_multireg}{formula}(formula, data=NULL, ...)

\method{GSest_multireg}{default}(X, Y, int = TRUE, bdp = 0.5, control=GScontrol(...),
na.action=na.omit, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{ an object of class \code{\link{formula}}; a symbolic description of the model to be fit.}
  \item{data}{ data frame from which variables specified in formula are to be taken.}
  \item{X}{ a matrix or data frame containing the explanatory variables. }
  \item{Y}{ a matrix or data frame containing the response variables. }
  \item{int}{ logical: if \code{TRUE} an intercept term is added to the model (unless it is already present in \code{X})}
  \item{bdp}{ required breakdown point. Should have \eqn{0 < } \code{bdp} \eqn{\le 0.5}, the default is 0.5. }
  \item{control}{a list with control parameters for tuning the computing algorithm, see \code{\link{GScontrol}}().}
  \item{na.action}{ a function which indicates what should happen when the data contain NAs. Defaults to \code{\link{na.omit}}.}
  \item{...}{allows for specifying control parameters directly instead of via \code{control}.}
}
\details{
%Called by \code{\link{FRBmultiregGS}} and typically not to be used on its own.
Generalized S-estimators are defined by minimizing the determinant of a robust estimator of the scatter matrix of 
the differences of the residuals. Hence, this procedure is intercept free and only gives an estimate for the slope matrix. To estimate 
the intercept, we use the M-type estimator of location of Lopuhaa (1992) on the residuals with the residual scatter matrix 
estimate of the residuals as a preliminary estimate.  We use a fast algorithm similar to the one proposed by Salibian-Barrera
and Yohai (2006) for the regression case. See \code{\link{GScontrol}} for the adjustable tuning parameters of this algorithm.

The returned object inherits from class \code{mlm} such that the standard \code{\link{coef}}, \code{\link{residuals}}, \code{\link{fitted}} and \code{\link[=predict.lm]{predict}} functions can be used.
}
\value{
  An object of class \code{FRBmultireg} which extends class \code{mlm} and contains at least the following components:
  \item{coefficients}{ GS-estimates of the regression coefficients}
%  \item{intercept}{ estimate for the intercept}
  \item{residuals}{ the residuals, that is response minus fitted values}
  \item{fitted.values}{ the fitted values.}
  \item{Sigma }{ GS-estimate of the error covariance matrix}
  \item{Gamma}{GS-estimate of the error shape matrix}
  \item{scale}{ GS-estimate of the size of the multivariate errors}
  \item{weights}{implicit weights corresponding to the GS-estimates (i.e. final weights in the RWLS procedure for the intercept estimate)}
  \item{outFlag}{outlier flags: 1 if the robust distance of the residual exceeds the .975 quantile of (the square root of)
  the chi-square distribution with degrees of freedom equal to the dimension of the responses; 0 otherwise}
  \item{b,c}{tuning parameters used in Tukey biweight loss function, as determined by \code{bdp}}
  \item{method}{ a list with following components: \code{est} = character string indicating that GS-estimates were used, and
  \code{bdp} = a copy of the \code{bdp} argument }
  \item{control}{ a copy of the \code{control} argument }
}
\references{
\itemize{ 
\item H.P. Lopuhaa (1992) Highly efficient estimators of multivariate location with high breakdown point. \emph{The Annals of Statistics}, \bold{20}, 
398-413. 
\item E. Roelant, S. Van Aelst and  C. Croux (2009) Multivariate Generalized S-estimators. \emph{Journal of Multivariate Analysis}, \bold{100}, 876--887.
\item M. Salibian-Barrera and V. Yohai (2006) A fast algorithm for S-regression estimates. 
\emph{Journal of Computational and Graphical Statistics}, \bold{15}, 414-427. 
}
}
\author{ Ella Roelant, Gert Willems and Stefan Van Aelst}
%\note{ ~~further notes~~ 

% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{diagplot.FRBmultireg}}, \code{\link{FRBmultiregGS}}, \code{\link{GSboot_multireg}},  \code{\link{Sest_multireg}}, \code{\link{GScontrol}} }
\examples{
data(schooldata)
school.x <- data.matrix(schooldata[,1:5])
school.y <- data.matrix(schooldata[,6:8])
GSest <- GSest_multireg(school.x,school.y)
# or using the formula interface
GSests <- GSest_multireg(cbind(reading,mathematics,selfesteem)~., data=schooldata)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
