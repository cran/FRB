\name{FRBmultiregS}
\alias{FRBmultiregS}
\alias{FRBmultiregS.default}
\alias{FRBmultiregS.formula}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{S-Estimates for Multivariate Regression with Bootstrap Inference}
\description{Computes S-estimates for multivariate regression together with standard errors, confidence intervals
and p-values based on the Fast and Robust Bootstrap.
}
\usage{
\method{FRBmultiregS}{formula}(formula, data=NULL, ...)

\method{FRBmultiregS}{default}(X, Y, int = TRUE, R = 999, bdp = 0.5, conf = 0.95, 
                control=Scontrol(...), na.action=na.omit, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{ an object of class \code{\link{formula}}; a symbolic description of the model to be fit.}
  \item{data}{ data frame from which variables specified in formula are to be taken.}
  \item{X}{ a matrix or data frame containing the explanatory variables. }
  \item{Y}{ a matrix or data frame containing the response variables. }
  \item{int}{ logical: if \code{TRUE} an intercept term is added to the model (unless it is already present in \code{X})}
  \item{R}{ number of bootstrap samples. Default is \code{R=999}. }
  \item{bdp}{ required breakdown point for the S-estimates. Should have \eqn{0 < } \code{bdp} \eqn{\le 0.5}, the default is 0.5 }
  \item{conf}{ level of the bootstrap confidence intervals. Default is \code{conf=0.95} }
  \item{control}{a list with control parameters for tuning the computing algorithm, see \code{\link{Scontrol}}().}
  \item{na.action}{ a function which indicates what should happen when the data contain NAs. Defaults to \code{\link{na.omit}}.}
  \item{...}{allows for specifying control parameters directly instead of via \code{control}}
}
\details{
Multivariate S-estimates were introduced by Davies (1987) and can be highly robust while enjoying a reasonable Gaussian efficiency.
Their use in the multivariate regression setting was discussed in Van Aelst and Willems (2005).
The loss function used here is Tukey's biweight. It is tuned in order to achieve the required breakdown point \code{bdp} 
(any value between 0 and 0.5). 
  
The computation is carried out by a call to \code{\link{Sest_multireg}}(), which performs the fast-S algorithm 
(Salibian-Barrera and Yohai 2006), see \code{\link{Scontrol}} for its tuning parameters. 
The result of this call is also returned as the value \code{est}.
%, a list containing:  beta (the matrix of coefficient estimates),
%covariance (scatter matrix estimate of the residuals),scale (GS-scale residuale estimate), b (constant that appears in constraint 
%to obtain consistent estimators) and c (constant used in Tukey's biweight function).

The Fast and Robust Bootstrap (Salibian-Barrera and Zamar 2002) is used to calculate so-called
basic bootstrap confidence intervals and bias corrected and accelerated (BCa)
confidence intervals (Davison and Hinkley 1997, p.194 and p.204 respectively). 
Apart from the intervals with the requested confidence level, the function also returns p-values for each coefficient
corresponding to the hypothesis that the actual coefficient is zero. The p-values are computed as 
1 minus the smallest level for which the confidence intervals would include zero. Both BCa and basic bootstrap p-values in this sense are given.
The bootstrap calculation is carried out by a call to \code{\link{Sboot_multireg}}(), the result 
of which is returned as the value \code{bootest}. Bootstrap standard errors are returned as well.

Note: Bootstrap samples which contain too few distinct observations with positive weights are discarded 
(a warning is given if this happens). The number of samples actually used is returned via \code{ROK}. 

%, a list containing: centered (fast and robust bootstrap recalculation centered by the original estimates,
%GSest (original GS-estimates stacked on top of each other),SE (bootstrap standard errors), CI.bca (bias corrected and accelerated confidence intervals), CI.basic (basic bootstap 
%intervals).

In the \code{formula}-interface, a multivariate response is produced via \code{cbind}. For example \code{cbind(x4,x5) ~ x1+x2+x3}.
All arguments from the default method can also be passed to the \code{formula} method except for \code{int} (passing \code{int} explicitely
will produce an error; the inclusion of an intercept term is determined by \code{formula}).

The returned object inherits from class \code{mlm} such that the standard \code{\link{coef}}, \code{\link{residuals}}, \code{\link{fitted}} and \code{\link[=predict.lm]{predict}} functions can be used.
}
\value{
  An object of class \code{FRBmultireg} which extends class \code{mlm} and contains at least the following components:
  \item{coefficients}{ MM-estimates of the regression coefficients}
%  \item{intercept}{ estimate for the intercept}
  \item{residuals}{ the residuals, that is response minus fitted values}
  \item{fitted.values}{ the fitted values.}
  \item{Sigma }{ S-estimate of the error covariance matrix}
  \item{scale}{ MM-estimate of the size of the multivariate errors}
  \item{weights}{ implicit weights corresponding to the S-estimates (i.e. final weights in the RWLS procedure at the end of the fast-S algorithm)}
  \item{outFlag}{outlier flags: 1 if the robust distance of the residual exceeds the .975 quantile of (the square root of)
  the chi-square distribution with degrees of freedom equal to the dimension of the responses; 0 otherwise}  
  \item{SE }{ bootstrap standard errors corresponding to the regression coefficients.}
  \item{cov}{ bootstrap covariance matrix corresponding to the regression coefficients (in vectorized form)}
  \item{CI.bca.lower}{ a matrix containing the lower bounds of the bias corrected and accelerated confidence intervals for the regression coefficients.}
  \item{CI.bca.upper}{ a matrix containing the upper bounds of the bias corrected and accelerated confidence intervals for the regression coefficients.}
  \item{CI.basic.lower}{ a matrix containing the lower bounds of basic bootstrap intervals for the regression coefficients.}
  \item{CI.basic.upper}{ a matrix containing the upper bounds of basic bootstrap intervals for the regression coefficients.}
  \item{p.bca}{ a matrix containing the p-values based on the BCa confidence intervals for the regression coefficients.}
  \item{p.basic}{ a matrix containing the p-values based on the basic bootstrap intervals for the regression coefficients.}
  \item{est }{S-estimates as returned by the call to \code{\link{Sest_multireg}}()}
  \item{bootest}{bootstrap results for the S-estimates as returned by the call to \code{\link{Sboot_multireg}}()}
  \item{conf}{ a copy of the \code{conf} argument}
  \item{method}{ a list with following components: \code{est} = character string indicating that S-estimates were used, and
  \code{bdp} = a copy of the \code{bdp} argument }
  \item{control}{ a copy of the \code{control} argument }
  \item{X, Y}{either copies of the respective arguments or the corresponding matrices produced from \code{formula}} 
  \item{ROK }{number of bootstrap samples actually used (i.e. not discarded due to too few distinct observations 
  with positive weight)}
}
\references{
\itemize{ 
\item P.L. Davies (1987) Asymptotic behavior of S-estimates of multivariate location parameters and dispersion
matrices. \emph{The Annals of Statistics}, \bold{15}, 1269-1292.
\item A.C. Davison and D.V. Hinkley (1997) \emph{Bootstrap Methods and their Application}. Cambridge Series in
Statistical and Probabilistic Mathematics. Cambridge: Cambridge University Press.
\item M. Salibian-Barrera, S. Van Aelst and G. Willems (2008) Fast and robust 
bootstrap. \emph{Statistical Methods and Applications}, \bold{17}, 41-71. 
\item M. Salibian-Barrera and V. Yohai (2006) A fast algorithm for S-regression estimates. 
\emph{Journal of Computational and Graphical Statistics}, \bold{15}, 414-427. 
\item M. Salibian-Barrera, R.H. Zamar (2002) Bootstrapping robust estimates of 
regression. \emph{The Annals of Statistics}, \bold{30}, 556-582.
\item  S. Van Aelst and G. Willems (2005) Multivariate regression S-estimators for robust estimation and
inference. \emph{Statistica Sinica}, \bold{15}, 981-1001.
\item S. Van Aelst and G. Willems (2013). Fast and robust bootstrap for multivariate inference: The R package FRB. \emph{Journal of Statistical Software}, \bold{53}(3), 1--32.
        URL: \url{http://www.jstatsoft.org/v53/i03/}.
}
}
\author{ Gert Willems, Stefan Van Aelst and Ella Roelant }
\seealso{ \code{\link{summary.FRBmultireg}}, \code{\link{plot.FRBmultireg}}, \code{\link{Sboot_multireg}}, 
\code{\link{Sest_multireg}},  \code{\link{FRBmultiregMM}},  \code{\link{FRBmultiregGS}}, \code{\link{Scontrol}} }
\examples{
data(schooldata)
school.x <- data.matrix(schooldata[,1:5])
school.y <- data.matrix(schooldata[,6:8])

#computes 25\% breakdown point S-estimate and 99\% confidence intervals 
#based on 999 bootstrap samples:
Sres <- FRBmultiregS(school.x, school.y, R=999, bdp = 0.25, conf = 0.99)
#or, equivalently using the formula interface
\dontrun{Sres <- FRBmultiregS(cbind(reading,mathematics,selfesteem)~., data=schooldata, 
        R=999, bdp = 0.25, conf = 0.99)}
          
#the print method displays the coefficient estimates 
Sres

#the summary function additionally displays the bootstrap standard errors and p-values
#("BCA" method by default)
summary(Sres)

summary(Sres, confmethod="basic")
                                                              
#ask explicitely for the coefficient matrix:
Sres$coefficients
# or equivalently,
coef(Sres)
#For the error covariance matrix:
Sres$Sigma
                                                              
#plot some bootstrap histograms for the coefficient estimates 
#(with "BCA" intervals by default) 
plot(Sres, expl=c("education", "occupation"), resp=c("selfesteem","reading"))

#plot bootstrap histograms for all coefficient estimates
plot(Sres)
#probably the plot-function has made a selection of coefficients to plot here, 
#since 'all' was too many to  fit on one page, see help(plot.FRBmultireg); 
#this is platform-dependent
%# this may give an error if the number of coefficients to plot is too large to fit on one page; 
%# use arguments expl and resp

%# diagnostic plot for outlier detection:
%plot(Sres, which=1)
%# this may take a while, since the function needs to compute S-estimates 
%# for the X matrix
}

\keyword{ multivariate }% __ONLY ONE__ keyword per line
\keyword{ robust }
