\name{GSboot_multireg}
\alias{GSboot_multireg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fast and Robust Bootstrap for GS-Estimates  }
\description{
  Calculates bootstrapped GS-estimates and bootstrap confidence intervals using the Fast and Robust Bootstrap
method. 
}
\usage{
GSboot_multireg(X, Y, R, conf=0.95, ests = GSest_multireg(X, Y))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ a matrix or data frame containing the explanatory variables. }
  \item{Y}{ a matrix or data frame containing the response variables. }
  \item{R}{ number of bootstrap samples }
  \item{conf}{ confidence level of the bootstrap confidence intervals. Default is \code{conf=0.95} }
  \item{ests}{ GS-estimates as returned by \code{\link{GSest_multireg}}() }
}
\details{Called by \code{\link{FRBmultiregGS}} and typically not to
be used on its own. If no original GS-estimates are provided the function uses the value of \code{\link{GSest_multireg}}. 

The fast and robust bootstrap was first developed by Salibian-Barrera and Zamar (2002) for univariate regression MM-estimators. 

The value \code{centered} gives a matrix with \code{R} columns and \eqn{p*q+q*q} rows (\eqn{p} is the number of explanatory variables without intercept 
and \eqn{q} is the number of response variables), containing the recalculated GS-estimates. Each column represents a different bootstrap sample.
The first \eqn{p*q} rows are the recalculated slope estimates and the next \eqn{q*q} rows are the covariance estimates 
(the estimates are vectorized, i.e. columns stacked on top of each other). These bootstrap estimates are centered by the original
estimates, which are also returned through \code{vecest} in vectorized form.

The output list further contains bootstrap standard errors, as well as so-called basic bootstrap confidence intervals and bias corrected and accelerated 
confidence intervals (Davison and Hinkley, 1997, p.194 and p.204 respectively). 
Also in the output are p-values defined as 1 minus the smallest confidence level for which the confidence intervals would include the 
(hypothesised) value of zero. Both BCa and basic bootstrap p-values are given. These are only useful for the
regression coefficient estimates (not really for the covariance estimates).

Bootstrap samples which contain too few distinct observations with positive weights are discarded 
(a warning is given if this happens). The number of samples actually used is returned via \code{ROK}. 
}
\value{
  A list containing the following components:
  \item{centered }{a matrix of all fast and robust bootstrap recalculations
             where the recalculations are centered by the original estimates (see Details). }
  \item{vecest }{a vector containing the orginal estimates stacked on top of each other}
  \item{SE }{bootstrap standard errors for the estimates in \code{vecest}}
  \item{CI.bca }{a matrix containing bias corrected and accelerated confidence intervals, corresponding to the 
  estimates in \code{vecest} (first column are lower limits, second column are upper limits)}
  \item{CI.basic }{ a matrix containing basic bootstrap intervals, corresponding to the 
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
\item E. Roelant, C. Croux and S. Van Aelst (2008) Multivariate Generalized S-estimators. To appear in Journal of Multivariate Analysis.
\item M. Salibian-Barrera, S. Van Aelst and G. Willems (2008) Fast and robust 
bootstrap. \emph{Statistical Methods and Applications}, \bold{17}, 41-71. 
\item M. Salibian-Barrera, R.H. Zamar (2002) Bootstrapping robust estimates of 
regression. \emph{The Annals of Statistics}, \bold{30}, 556-582.
}
}
\author{ Ella Roelant and Gert Willems }
%\note{ ~~further notes~~ 

% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{GSest_multireg}}}
\examples{

data(schooldata)
school.x <- data.matrix(schooldata[,1:5])
school.y <- data.matrix(schooldata[,6:8])

#computes 10 bootstrap recalculations starting from the GS-estimator
#obtained from GSest_multireg
bootres <- GSboot_multireg(school.x,school.y,R=10)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
