\name{FRBmultiregS}
\alias{FRBmultiregS}
\alias{FRBmultiregS.default}
\alias{FRBmultiregS.formula}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{S-Estimates for Multivariate Regression with Bootstrap Inference}
\description{Computes S-estimates for multivariate regression together with standard errors and confidence intervals
based on the Fast and Robust Bootstrap.
}
\usage{
\method{FRBmultiregS}{formula}(formula, data, ...)

\method{FRBmultiregS}{default}(X, Y, int = TRUE, R = 999, bdp = 0.5, conf = 0.95, 
                control=Scontrol(...), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{ an object of class \code{\link{formula}}; a symbolic description of the model to be fit.}
  \item{data}{ data frame from which variables specified in formula are to be taken.}
  \item{X}{ a matrix or data frame containing the explanatory variables. }
  \item{Y}{ a matrix or data frame containing the response variables. }
  \item{int}{ logical: if \code{TRUE} an intercept term is added to the model (unless it is already present in \code{X})}
  \item{R}{ number of bootstrap samples }
  \item{bdp}{ required breakdown point for the S-estimates. Should have \eqn{0 < } \code{bdp} \eqn{\le 0.5}, the default is 0.5 }
  \item{conf}{ level of the bootstrap confidence intervals. Default is \code{conf=0.95} }
  \item{control}{a list with control parameters for tuning the computing algorithm, see \code{\link{Scontrol}}().}
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
basic bootstrap confidence intervals and bias corrected and accelerated 
confidence intervals (Davison and Hinkley 1997, p.194 and p.204 respectively). This computation is carried out by a call to \code{\link{Sboot_multireg}}(), the result 
of which is returned as the value \code{bootest}. Bootstrap standard errors are returned as well.
%, a list containing: centered (fast and robust bootstrap recalculation centered by the original estimates,
%GSest (original GS-estimates stacked on top of each other),SE (bootstrap standard errors), CI.bca (bias corrected and accelerated confidence intervals), CI.basic (basic bootstap 
%intervals).

In the \code{formula}-interface, a multivariate response is produced via \code{cbind}. For example \code{cbind(x4,x5) ~ x1+x2+x3}.
All arguments from the default method can also be passed to the \code{formula} method except for \code{int} (passing \code{int} explicitely
will produce an error; the inclusion of an intercept term is determined by \code{formula}).
}
\value{
  An object of class \code{FRBmultireg}, which is a list containing the following components:
  \item{Beta }{ S-estimate for the regression coefficients}
  \item{Sigma }{ S-estimate for the error covariance matrix}
  \item{SE }{ bootstrap standard errors corresponding to the elements in \code{Beta}}
  \item{CI.bca.lower}{ a matrix containing the lower bounds of the bias corrected and accelerated confidence intervals for each element in \code{Beta}.}
  \item{CI.bca.upper}{ a matrix containing the upper bounds of the bias corrected and accelerated confidence intervals for each element in \code{Beta}.}
  \item{CI.basic.lower}{ a matrix containing the lower bounds of basic bootstrap intervals for each element in \code{Beta}.}
  \item{CI.basic.upper}{ a matrix containing the upper bounds of basic bootstrap intervals for each element in \code{Beta}.}
  \item{est }{S-estimates as returned by the call to \code{\link{Sest_multireg}}()}
  \item{bootest}{bootstrap results for the S-estimates as returned by the call to \code{\link{Sboot_multireg}}()}
  \item{conf}{ a copy of the \code{conf} argument}
  \item{method}{ a list with following components: \code{est} = character string indicating that S-estimates were used, and
  \code{bdp} = a copy of the \code{bdp} argument }
  \item{control}{ a copy of the \code{control} argument }
  \item{X, Y}{either copies of the respective arguments or the corresponding matrices produced from \code{formula}} 
}
\references{
\itemize{ 
\item P.L. Davies (1987). Asymptotic behavior of S-estimates of multivariate location parameters and dispersion
matrices. \emph{The Annals of Statistics}, \bold{15}, 1269-1292.
\item A.C. Davison and D.V. Hinkley (1997). \emph{Bootstrap Methods and their Application}. Cambridge Series in
Statistical and Probabilistic Mathematics. Cambridge: Cambridge University Press.
\item M. Salibian-Barrera, S. Van Aelst and G. Willems (2008) Fast and robust 
bootstrap. \emph{Statistical Methods and Applications}, \bold{17}, 41-71. 
\item M. Salibian-Barrera and V. Yohai (2006) A fast algorithm for S-regression estimates. 
\emph{Journal of Computational and Graphical Statistics}, \bold{15}, 414-427. 
\item M. Salibian-Barrera, R.H. Zamar (2002) Bootstrapping robust estimates of 
regression. \emph{The Annals of Statistics}, \bold{30}, 556-582.
\item  S. Van Aelst and G. Willems (2005). Multivariate regression S-estimators for robust estimation and
inference. \emph{Statistica Sinica}, \bold{15}, 981-1001.
}
}
\author{ Gert Willems and Ella Roelant }
\seealso{ \code{\link{summary.FRBmultireg}}, \code{\link{print.FRBmultireg}}, \code{\link{plot.FRBmultireg}}, \code{\link{Sboot_multireg}}, 
\code{\link{Sest_multireg}},  \code{\link{FRBmultiregMM}},  \code{\link{FRBmultiregGS}}, \code{\link{Scontrol}} }
\examples{
data(schooldata)
school.x <- data.matrix(schooldata[,1:5])
school.y <- data.matrix(schooldata[,6:8])

#computes 25\% breakdown point S-estimate and 99\% confidence intervals 
#based on 999 bootstrap samples:
Sres <- FRBmultiregS(school.x, school.y, R=999, bdp = 0.25, conf = 0.99)
#or, equivalently,
Sres <- FRBmultiregS(cbind(reading,mathematics,selfesteem)~., data=schooldata, 
        R=999, bdp = 0.25, conf = 0.99)
          
#the print method displays the coefficients with their bootstrap standard errors
Sres

#the summary function additionally displays the confidence intervals 
#("BCA" method by default)
summary(Sres)

summary(Sres, confmethod="basic")
                                                              
#ask explicitely for the coefficient matrix:
Sres$Beta
#or for the error covariance matrix:
Sres$Sigma
                                                              
#plot some bootstrap histograms for the coefficient estimates 
#(with "BCA" intervals by default) 
plot(Sres, which=2, expl=c("education", "occupation"), resp=c("selfesteem","reading"))

#plot bootstrap histograms for all coefficient estimates
plot(Sres, which=2)
#probably the plot-function has made a selection of coefficients to plot here, 
#since 'all' was too many to  fit on one page, see help(plot.FRBmultireg); 
#this is platform-dependent
%# this may give an error if the number of coefficients to plot is too large to fit on one page; 
%# use arguments expl and resp

# diagnostic plot for outlier detection:
plot(Sres, which=1)
# this may take a while, since the function needs to compute S-estimates 
# for the X matrix
}

\keyword{ multivariate }% __ONLY ONE__ keyword per line
\keyword{ robust }
