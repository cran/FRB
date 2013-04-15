\name{summary.FRBmultireg}
\alias{summary.FRBmultireg}
\alias{print.summary.FRBmultireg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Summary Method for Objects of Class 'FRBmultireg' }
\description{
Summary method for objects of class \code{FRBmultireg}, and print method of the summary object.
}
\usage{
\method{summary}{FRBmultireg}(object, confmethod = c("BCA", "basic", "both"), digits = 3, 
print.CI=FALSE, sep="", ...)

\method{print}{summary.FRBmultireg}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ an \R object of class \code{FRBmultireg}, typically created by \code{\link{FRBmultiregS}}, \code{\link{FRBmultiregMM}} or \code{\link{FRBmultiregGS}} }
  \item{confmethod}{ which kind of bootstrap confidence intervals to be displayed: 
'BCA'= bias corrected and accelerated method, 'basic'= basic bootstrap method, 'both'=both kinds of confidence intervals }
  \item{digits}{ number of digits for printing (default is 3) }
  \item{print.CI}{ logical: Should Confidence intervals be printed?}
  \item{sep}{ Symmbol to separate columns in output. Default is \code{""} }
  \item{x}{ an \R object of class \code{summary.FRBmultireg}, resulting for example from \code{summary(\link{FRBmultiregS}(),...)} }
  \item{\dots}{ potentially more arguments to be passed to methods }
}
\details{
The \code{print} method displays in a \dQuote{familiar way} the components of the \code{summary} object, which are listed in the Value section. 
}
\value{
  \code{summary} returns an object of class \code{summary.FRBmultireg}, which contains the following components:
  \item{responses}{the names of the response variables in the fitted model}
  \item{covariates}{the names of the covariates (predictors) in the fitted model}
  \item{Betawstd}{a data frame containing the coefficient estimates and their bootstrap standard errors}
  \item{Sigma}{estimate for the error covariance matrix}
  \item{table.bca}{a list with for each response variable a matrix containing the estimates, standard errors, lower and upper limits of the 
      BCa confidence intervals, p-values and a significance code (only present when \code{confmethod="BCA"} or \code{confmethod="both"})}
  \item{table.basic}{a list with for each response variable a matrix containing the estimates, standard errors, lower and upper limits of the 
      basic bootstrap confidence intervals, p-values and a significance code (only present when \code{confmethod="basic"} or \code{confmethod="both"})}
  \item{method }{multivariate regression method that was used}
  \item{conf}{confidence level that was used}
  \item{digits}{ number of digits for printing }
}
\references{ 
\itemize{
\item S. Van Aelst and G. Willems (2013). Fast and robust bootstrap for multivariate inference: The R package FRB. \emph{Journal of Statistical Software}, \bold{53}(3), 1--32.
        URL: \url{http://www.jstatsoft.org/v53/i03/}.
        }
}
\author{ Gert Willems, Ella Roelant and Stefan Van Aelst }
\seealso{ \code{\link{FRBmultiregS}}, \code{\link{FRBmultiregMM}}, \code{\link{FRBmultiregGS}},  \code{\link{print.FRBmultireg}}, \code{\link{plot.FRBmultireg}} }
\examples{

data(schooldata)

MMres <- FRBmultiregMM(cbind(reading,mathematics,selfesteem)~., data=schooldata,
R=199, conf = 0.99,nsamp=200)
summary(MMres)  # -> print.summary.FRBmultireg() method

GSres <- FRBmultiregGS(cbind(reading,mathematics,selfesteem)~., data=schooldata, 
bdp = 0.25,R=199,nsamp=50)
summary(GSres, confmethod="both")  # -> print.summary.FRBmultireg() method

}
