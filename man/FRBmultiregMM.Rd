\name{FRBmultiregMM}
\alias{FRBmultiregMM}
\alias{FRBmultiregMM.default}
\alias{FRBmultiregMM.formula}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{MM-Estimates for Multivariate Regression with Bootstrap Inference}
\description{Computes MM-estimates for multivariate regression together with standard errors, confidence intervals
and p-values based on the Fast and Robust Bootstrap.
}
\usage{
\method{FRBmultiregMM}{formula}(formula, data, ...)

\method{FRBmultiregMM}{default}(X, Y, int = TRUE, R = 999, conf = 0.95, 
                control=MMcontrol(...), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{ an object of class \code{\link{formula}}; a symbolic description of the model to be fit.}
  \item{data}{ data frame from which variables specified in formula are to be taken.}
  \item{X}{ a matrix or data frame containing the explanatory variables. }
  \item{Y}{ a matrix or data frame containing the response variables. }
  \item{int}{ logical: if \code{TRUE} an intercept term is added to the model (unless it is already present in \code{X})}
  \item{R}{ number of bootstrap samples }
  \item{conf}{ level of the bootstrap confidence intervals. Default is \code{conf=0.95} }
  \item{control}{a list with control parameters for tuning the MM-estimate and its computing algorithm, 
  see \code{\link{MMcontrol}}().}
  \item{...}{allows for specifying control parameters directly instead of via \code{control}}
}
\details{
Multivariate MM-estimates combine high breakdown point and high Gaussian efficiency. They are defined by first computing an S-estimate of regression, 
then fixing the scale component of the error covariance estimate, and finally re-estimating the regression coefficients 
and the shape part of the error covariance by a more efficient M-estimate (see Tatsuoka and Tyler (2000) for MM-estimates 
in the special case of location/scatter estimation, and Van Aelst and Willems (2005) for S-estimates of multivariate regression). 
  
Tukey's biweight is used for the loss functions. By default, the first loss function (in the S-estimate) is tuned in order to obtain 50\% breakdown point. 
The default tuning of the second loss function (M-estimate) ensures 95\% efficiency at the normal model for the coefficient estimates. 
This tuning is recommended but can be changed through argument \code{control} if desired. 

The computation is carried out by a call to \code{\link{MMest_multireg}}(), which first performs the fast-S algorithm 
(see \code{\link{Sest_multireg}}) and does the M-part by reweighted least squares (RWLS) iteration. 
See \code{\link{MMcontrol}} for some adjustable tuning parameters regarding the algorithm.
The result of this call is also returned as the value \code{est}.

The Fast and Robust Bootstrap (Salibian-Barrera and Zamar 2002) is used to calculate so-called
basic bootstrap confidence intervals and bias corrected and accelerated (BCa)
confidence intervals (Davison and Hinkley 1997, p.194 and p.204 respectively). 
Apart from the intervals with the requested confidence level, the function also returns p-values for each coefficient
corresponding to the hypothesis that the actual coefficient is zero. The p-values are computed as 
1 minus the smallest level for which the confidence intervals would include zero. Both BCa and basic bootstrap p-values in this sense are given.
The bootstrap calculation is carried out by a call to \code{\link{MMboot_multireg}}(), the result 
of which is returned as the value \code{bootest}. Bootstrap standard errors are returned as well.

Note: Bootstrap samples which contain too few distinct observations with positive weights are discarded 
(a warning is given if this happens). The number of samples actually used is returned via \code{ROK}. 

In the \code{formula}-interface, a multivariate response is produced via \code{cbind}. For example \code{cbind(x4,x5) ~ x1+x2+x3}.
All arguments from the default method can also be passed to the \code{formula} method except for \code{int} (passing \code{int} explicitely
will produce an error; the inclusion of an intercept term is determined by \code{formula}).
}
\value{
  An object of class \code{FRBmultireg}, which is a list containing the following components:
  \item{Beta }{ MM-estimate for the regression coefficients}
  \item{Sigma }{ MM-estimate for the error covariance matrix}
  \item{SE }{ bootstrap standard errors corresponding to the elements in \code{Beta}}
  \item{CI.bca.lower}{ a matrix containing the lower bounds of the bias corrected and accelerated confidence intervals for each element in \code{Beta}.}
  \item{CI.bca.upper}{ a matrix containing the upper bounds of the bias corrected and accelerated confidence intervals for each element in \code{Beta}.}
  \item{CI.basic.lower}{ a matrix containing the lower bounds of basic bootstrap intervals for each element in \code{Beta}.}
  \item{CI.basic.upper}{ a matrix containing the upper bounds of basic bootstrap intervals for each element in \code{Beta}.}
  \item{p.bca}{ a matrix containing the p-values based on the BCa confidence intervals for each element in \code{Beta}.}
  \item{p.basic}{ a matrix containing the p-values based on the basic bootstrap intervals for each element in \code{Beta}.}
  \item{est }{MM-estimates as returned by the call to \code{\link{MMest_multireg}}()}
  \item{bootest}{bootstrap results for the MM-estimates as returned by the call to \code{\link{MMboot_multireg}}()}
  \item{conf}{ a copy of the \code{conf} argument}
  \item{method}{ a list with following components: \code{est} = character string indicating that MM-estimates were used, 
  \code{bdp} = a copy of \code{bdp} from the \code{control} argument, and \code{eff} = a copy of \code{eff} from the \code{control} argument}
  \item{control}{ a copy of the \code{control} argument }
  \item{X, Y}{either copies of the respective arguments or the corresponding matrices produced from \code{formula}} 
  \item{ROK}{number of bootstrap samples actually used (i.e. not discarded due to too few distinct observations 
  with positive weight)}
  \item{w}{implicit weights corresponding to the MM-estimates (i.e. final weights in the RWLS procedure)}
  \item{outFlag}{outlier flags: 1 if the robust distance of the residual exceeds the .975 quantile of (the square root of)
  the chi-square distribution with degrees of freedom equal to the dimension of the responses; 0 otherwise}
}
\references{ 
\itemize{
\item A.C. Davison and D.V. Hinkley (1997). \emph{Bootstrap Methods and their Application}. Cambridge Series in
Statistical and Probabilistic Mathematics. Cambridge: Cambridge University Press.
\item M. Salibian-Barrera, S. Van Aelst and G. Willems (2008) Fast and robust 
bootstrap. \emph{Statistical Methods and Applications}, \bold{17}, 41-71. 
\item M. Salibian-Barrera, R.H. Zamar (2002) Bootstrapping robust estimates of 
regression. \emph{The Annals of Statistics}, \bold{30}, 556-582.
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
\seealso{ \code{\link{summary.FRBmultireg}}, \code{\link{print.FRBmultireg}}, \code{\link{plot.FRBmultireg}}, \code{\link{MMboot_multireg}}, 
\code{\link{MMest_multireg}},  \code{\link{FRBmultiregS}},  \code{\link{FRBmultiregGS}}, \code{\link{MMcontrol}} }
\examples{
data(schooldata)
school.x <- data.matrix(schooldata[,1:5])
school.y <- data.matrix(schooldata[,6:8])

#computes MM-estimate and 95\% confidence intervals 
#based on 999 bootstrap samples:
MMres <- FRBmultiregMM(school.x, school.y, R=999, conf = 0.95)
#or, equivalently,
MMres <- FRBmultiregMM(cbind(reading,mathematics,selfesteem)~., data=schooldata, 
              R=999, conf = 0.95)

#the print method displays the coefficients with their bootstrap standard errors
MMres

#the summary function additionally displays the confidence intervals and p-values
#("BCA" method by default)
summary(MMres)

summary(MMres, confmethod="basic")

#ask explicitely for the coefficient matrix:
MMres$Beta
#or for the error covariance matrix:
MMres$Sigma
                                                              
#plot some bootstrap histograms for the coefficient estimates 
#(with "BCA" intervals by default) 
plot(MMres, which=2, expl=c("education", "occupation"), resp=c("selfesteem","reading"))

#plot bootstrap histograms for all coefficient estimates
plot(MMres, which=2)
#probably the plot-function has made a selection of coefficients to plot here, 
#since 'all' was too many to  fit on one page, see help(plot.FRBmultireg); 
#this is platform-dependent
%# this may give an error if the number of coefficients to plot is too large to fit on one page; 
%# use arguments expl and resp

# diagnostic plot for outlier detection:
plot(MMres, which=1)
# this may take a while, since the function needs to compute MM-estimates 
# for the X matrix
}

\keyword{ multivariate }% __ONLY ONE__ keyword per line
\keyword{ robust }
