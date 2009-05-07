\name{FRBmultiregGS}
\alias{FRBmultiregGS}
\alias{FRBmultiregGS.default}
\alias{FRBmultiregGS.formula}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ GS-Estimates for multivariate regression with bootstrap confidence 
intervals }
\description{Computes GS-estimates for multivariate regression together with standard errors and confidence intervals
based on the Fast and Robust Bootstrap.
}
\usage{
\method{FRBmultiregGS}{formula}(formula, data, ...)

\method{FRBmultiregGS}{default}(X, Y, R = 999, bdp = 0.5, conf = 0.95, 
                control=GScontrol(...), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{ an object of class \code{\link{formula}}; a symbolic description of the model to be fit.}
  \item{data}{ data frame from which variables specified in formula are to be taken.}
  \item{X}{ a matrix or data frame containing the explanatory variables. }
  \item{Y}{ a matrix or data frame containing the response variables. }
  \item{R}{ number of bootstrap samples. }
  \item{bdp}{ required breakdown point. Should have \eqn{0 < } \code{bdp} \eqn{\le 0.5}, the default is 0.5. }
  \item{conf}{ confidence level of the bootstrap confidence intervals. Default is \code{conf=0.95}. }
  \item{control}{a list with control parameters for tuning the computing algorithm, see \code{\link{GScontrol}}().}
  \item{...}{allows for specifying control parameters directly instead of via \code{control}.}
}
\details{  Generalized S-estimators are defined by minimizing the determinant of a robust estimator of the scatter matrix of 
the differences of the residuals. Hence, this procedure is intercept free and only gives an estimate for the slope matrix. To estimate 
the intercept, we use the M-type estimator of location of Lopuhaa (1992) on the residuals with the residual scatter matrix 
estimate of the residuals as a preliminary estimate. This computation is carried out by a call to \code{\link{GSest_multireg}}(), 
which uses a fast-S-type algorithm (its tuning parameters can be changed via the \code{control} argument).  
The result of this call is also returned as the value \code{est}.
%The result is also given in the output as the value GSest, a list containing:  beta (the matrix of coefficient estimates),
%covariance (scatter matrix estimate of the residuals),scale (GS-scale residuale estimate), b (constant that appears in constraint 
%to obtain consistent estimators) and c (constant used in Tukey's biweight function).

The fast and robust bootstrap is used to calculate so-called basic bootstrap confidence intervals and bias corrected and accelerated 
confidence intervals (Davison and Hinkley, 1997, p.194 and p.204 respectively). This computation is carried out by a call to \code{\link{GSboot_multireg}}(), the result 
of which is returned as the value \code{bootest}. Bootstrap standard errors are returned as well.
% The result 
%is given in the output as the value GSboot, a list containing: centered (fast and robust bootstrap recalculation centered by the original estimates,
%GSest (original GS-estimates stacked on top of each other),SE (bootstrap standard errors), CI.bca (bias corrected and accelerated confidence intervals), CI.basic (basic bootstap 
%intervals).

In the \code{formula}-interface, a multivariate response is produced via \code{cbind}. For example \code{cbind(x4,x5) ~ x1+x2+x3}.
All arguments from the default method can also be passed to the \code{formula} method.
}
\value{
  An object of class \code{FRBmultireg}, which is a list containing the following components:
  \item{Beta }{ GS-estimate for slope}
  \item{intercept}{ estimate for the intercept}
  \item{Sigma }{ GS-estimate for the error covariance matrix}
  \item{SE }{ bootstrap standard errors corresponding to the elements in \code{Beta}}
  \item{CI.bca.lower}{ a matrix containing the lower bound of the bias corrected and accelerated confidence intervals for each element of \code{Beta}}
  \item{CI.bca.upper}{ a matrix containing the upper bound of the bias corrected and accelerated confidence intervals for each element of \code{Beta}}
  \item{CI.basic.lower}{ a matrix containing the lower bound of basic bootstrap intervals for each element of \code{Beta}}
  \item{CI.basic.upper}{ a matrix containing the upper bound of basic bootstrap intervals for each element of \code{Beta}}
  \item{est }{GS-estimates as returned by the call to \code{\link{GSest_multireg}}()}
  \item{bootest}{bootstrap results for the GS-estimates as returned by the call to \code{\link{GSboot_multireg}}()}
  \item{conf}{ a copy of the \code{conf} argument}
  \item{method}{ a list with following components: \code{est} = character string indicating that GS-estimates were used, and
  \code{bdp} = a copy of the \code{bdp} argument }
  \item{control}{ a copy of the \code{control} argument }
  \item{X, Y}{either copies of the respective arguments or the corresponding matrices produced from \code{formula}} 
}
\references{ 
\itemize{
\item A.C. Davison and D.V. Hinkley (1997). \emph{Bootstrap Methods and their Application}. Cambridge Series in
Statistical and Probabilistic Mathematics. Cambridge: Cambridge University Press.
\item H.P. Lopuhaa (1992) Highly efficient estimators of multivariate location with high breakdown point. \emph{The Annals of Statistics}, \bold{20}, 
398-413. 
\item E. Roelant, C. Croux and S. Van Aelst (2008) Multivariate Generalized S-estimators. To appear in Journal of Multivariate Analysis.
\item M. Salibian-Barrera, S. Van Aelst and G. Willems (2008) Fast and robust 
bootstrap. \emph{Statistical Methods and Applications}, \bold{17}, 41-71. 
}
}
\author{ Ella Roelant and Gert Willems }
%\note{ ~~further notes~~ 

% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{summary.FRBmultireg}}, \code{\link{print.FRBmultireg}}, \code{\link{plot.FRBmultireg}}, \code{\link{GSboot_multireg}}, \code{\link{GSest_multireg}},  
\code{\link{FRBmultiregMM}},  \code{\link{FRBmultiregS}}, \code{\link{GScontrol}} }
\examples{
data(schooldata)
school.x <- data.matrix(schooldata[,1:5])
school.y <- data.matrix(schooldata[,6:8])

#computes 25\% breakdown point GS-estimate and 99\% confidence intervals 
#based on 999 bootstrap samples:
GSres <- FRBmultiregGS(school.x, school.y, R=999, bdp = 0.25, conf = 0.99)
#or, equivalently,
GSres <- FRBmultiregGS(cbind(reading,mathematics,selfesteem)~., data=schooldata, 
          R=999, bdp = 0.25, conf = 0.99)

#the print method displays the coefficients with their bootstrap standard errors
GSres

#the summary function additionally displays the confidence intervals 
#("BCA" method by default)
summary(GSres)

summary(GSres, confmethod="basic")

#ask explicitely for the coefficient matrix:
GSres$Beta
#or for the error covariance matrix:
GSres$Sigma
                                                              
#plot some bootstrap histograms for the coefficient estimates 
#(with "BCA" intervals by default) 
plot(GSres, which=2, expl=c("education", "occupation"), resp=c("selfesteem","reading"))

#plot bootstrap histograms for all coefficient estimates
plot(GSres, which=2)
#possibly the plot-function has made a selection of coefficients to plot here, 
#since 'all' may have been too many to fit on one page, see help(plot.FRBmultireg); 
#this is platform-dependent
%# this may give an error if the number of coefficients to plot is too large to fit on one page; 
%# use arguments expl and resp

# diagnostic plot for outlier detection:
plot(GSres, which=1)
# this may take a while, since the function needs to compute GS-estimates 
# for the X matrix
}

\keyword{ multivariate }% __ONLY ONE__ keyword per line
\keyword{ robust }
