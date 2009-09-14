\name{Sboot_multireg}
\alias{Sboot_multireg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fast and Robust Bootstrap for S-Estimates of Multivariate Regression }
\description{
Calculates bootstrapped S-estimates of multivariate regression and corresponding bootstrap confidence intervals using the Fast and Robust Bootstrap
method. 
}
\usage{
Sboot_multireg(X, Y, R, conf=0.95, ests = Sest_multireg(X, Y))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ a matrix or data frame containing the explanatory variables (possibly including intercept). }
  \item{Y}{ a matrix or data frame containing the response variables. }
  \item{R}{ number of bootstrap samples }
  \item{conf}{ level of the bootstrap confidence intervals. Default is \code{conf=0.95} }
  \item{ests}{ S-estimates as returned by \code{\link{Sest_multireg}}() }
}
\details{
Called by \code{\link{FRBmultiregS}} and typically not to
be used on its own. It requires the result of \code{\link{Sest_multireg}} applied on \code{X} and \code{Y}, 
supplied through the argument \code{ests}. If \code{ests} is not provided, \code{\link{Sest_multireg}} 
will be called with default arguments. 

The fast and robust bootstrap was first developed by Salibian-Barrera and Zamar (2002) for univariate regression MM-estimators. 

The value \code{centered} gives a matrix with \code{R} columns and \eqn{p*q+q*q} rows (\eqn{p} is the number of explanatory variables 
and \eqn{q} the number of response variables), 
containing the recalculated S-estimates of the regression coefficients and the error covariance matrix. 
Each column represents a different bootstrap sample. 
The first \eqn{p*q} rows are the coefficient estimates, the next \eqn{q*q} rows represent the covariance estimate 
(the estimates are vectorized, i.e. columns stacked on top of each other). 
The estimates are centered by the original estimates, which are also returned through \code{vecest} in vectorized form.

The output list further contains bootstrap standard errors, as well as so-called basic bootstrap confidence intervals and bias corrected and accelerated (BCa)
confidence intervals (Davison and Hinkley, 1997, p.194 and p.204 respectively). 
Also in the output are p-values defined as 1 minus the smallest confidence level for which the confidence intervals would include the 
(hypothesised) value of zero. Both BCa and basic bootstrap p-values are given. These are only useful for the
regression coefficient estimates (not really for the covariance estimates).

Bootstrap samples which contain less than \eqn{p} distinct observations with positive weights are discarded 
(a warning is given if this happens). The number of samples actually used is returned via \code{ROK}. 
}
\value{
  A list containing the following components:
  \item{centered }{a matrix of all fast/robust bootstrap recalculations
             where the recalculations are centered by original estimates (see Details) }
  \item{vecest }{a vector containing the original estimates (see Details)}
  \item{SE }{bootstrap standard errors for the estimates in \code{vecest}}
  \item{CI.bca }{a matrix containing bias corrected and accelerated confidence intervals corresponding to the 
  estimates in \code{vecest} (first column are lower limits, second column are upper limits)}
  \item{CI.basic }{ a matrix containing basic bootstrap intervals corresponding to the 
  estimates in \code{vecest} (first column are lower limits, second column are upper limits)}
  \item{p.bca }{a vector containing p-values based on the bias corrected and accelerated confidence intervals (corresponding to the 
  estimates in \code{vecest})}
  \item{p.basic }{a vector containing p-values based on the basic bootstrap intervals (corresponding to the 
  estimates in \code{vecest})}
    \item{ROK }{number of bootstrap samples actually used (i.e. not discarded due to too few distinct observations 
  with positive weight)}
}
\references{ 
\itemize{
\item A.C. Davison, D.V. Hinkley (1997) Bootstrap methods and their application. Cambridge University Press.  
\item M. Salibian-Barrera, S. Van Aelst and G. Willems (2008) Fast and robust 
bootstrap. \emph{Statistical Methods and Applications}, \bold{17}, 41-71. 
\item M. Salibian-Barrera, R.H. Zamar (2002) Bootstrapping robust estimates of 
regression. \emph{The Annals of Statistics}, \bold{30}, 556-582.
\item  S. Van Aelst and G. Willems (2005). Multivariate regression S-estimators for robust estimation and
inference. \emph{Statistica Sinica}, \bold{15}, 981-1001
}
}
\author{ Gert Willems and Ella Roelant }
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{FRBmultiregS}}, \code{\link{Sest_multireg}}, \code{\link{MMboot_multireg}} }
\examples{
data(schooldata)
school.x <- data.matrix(schooldata[,1:5])
school.y <- data.matrix(schooldata[,6:8])

#computes 1000 bootstrap recalculations starting from the S-estimator
#obtained from Sest_multireg()
bootres <- Sboot_multireg(school.x,school.y,R=1000)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
