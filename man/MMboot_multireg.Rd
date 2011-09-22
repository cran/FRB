\name{MMboot_multireg}
\alias{MMboot_multireg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fast and Robust Bootstrap for MM-Estimates of Multivariate Regression }
\description{
  Calculates bootstrapped MM-estimates of multivariate regression and corresponding bootstrap confidence intervals using the Fast and Robust Bootstrap
method. 
}

\usage{
MMboot_multireg(X, Y, R = 999, conf=0.95, ests = MMest_multireg(X, Y))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ a matrix or data frame containing the explanatory variables (possibly including intercept). }
  \item{Y}{ a matrix or data frame containing the response variables. }
  \item{R}{ number of bootstrap samples. Default is \code{R=999}. }
  \item{conf}{ level of the bootstrap confidence intervals. Default is \code{conf=0.95}. }
  \item{ests}{ MM-estimates as returned by \code{\link{MMest_multireg}}(). }
}
\details{
Called by \code{\link{FRBmultiregMM}} and typically not to
be used on its own. It requires the result of \code{\link{MMest_multireg}} applied on \code{X} and \code{Y}, 
supplied through the argument \code{ests}. If \code{ests} is not provided, \code{\link{MMest_multireg}} 
will be called with default arguments. 

The fast and robust bootstrap was first developed by Salibian-Barrera and Zamar (2002) for univariate regression MM-estimators and extended to multivariate regression by Van Aelst and Willems (2005). 

The value \code{centered} gives a matrix with \code{R} columns and \eqn{2*(p*q+q*q)} rows (\eqn{p} is the number of explanatory variables 
and \eqn{q} the number of response variables), 
containing the recalculated MM-estimates and initial S-estimates. 
Each column represents a different bootstrap sample. 

The first \eqn{p*q} rows are the MM-coefficient estimates, the next \eqn{q*q} rows represent the MM-estimate of the error shape matrix 
(having determinant 1). 
Then the next \eqn{q*q} rows are the S-estimate of error covariance and the final \eqn{p*q} rows are the S-estimates of the regression
coefficients (all estimates are vectorized, i.e. columns stacked on top of each other). 
These estimates are centered by the original estimates, which are also returned through \code{vecest} in vectorized form.

The output list further contains bootstrap standard errors, as well as so-called basic bootstrap confidence intervals and bias corrected and accelerated 
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
  \item{cov}{ bootstrap covariance matrix for the estimates in \code{vecest}} 
  \item{CI.bca }{a matrix containing 95\% bias corrected and accelerated confidence intervals corresponding to the 
  estimates in \code{vecest} (first column are lower limits, second column are upper limits)}
  \item{CI.basic }{ a matrix containing 95\% basic bootstrap intervals corresponding to the 
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
bootstrap. \emph{Statistical Methods and Applications}, \bold{17}, 41--71. 
\item M. Salibian-Barrera, R.H. Zamar (2002) Bootstrapping robust estimates of 
regression. \emph{The Annals of Statistics}, \bold{30}, 556--582.
\item  S. Van Aelst and G. Willems (2005) Multivariate regression S-estimators for robust estimation and
inference. \emph{Statistica Sinica}, \bold{15}, 981--1001.
}
}
\author{ Gert Willems, Ella Roelant and Stefan Van Aelst}
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{  \code{\link{FRBmultiregMM}}, \code{\link{MMest_multireg}}, \code{\link{Sboot_multireg}} }
\examples{
data(schooldata)
school.x <- data.matrix(schooldata[,1:5])
school.y <- data.matrix(schooldata[,6:8])

#computes 1000 bootstrap recalculations starting from the MM-estimator
#obtained from MMest_multireg()
bootres <- MMboot_multireg(school.x,school.y,R=1000)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
