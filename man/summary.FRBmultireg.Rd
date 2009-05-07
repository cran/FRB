\name{summary.FRBmultireg}
\alias{summary.FRBmultireg}
\alias{print.summary.FRBmultireg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Summary Method for Objects of Class 'FRBmultireg' }
\description{
Summary method for objects of class \code{FRBmultireg}, and print method of the summary object.
}
\usage{
\method{summary}{FRBmultireg}(object, confmethod = c("BCA", "basic", "both"), ...)
\method{print}{summary.FRBmultireg}(x, digits = 3, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ an \R object of class \code{FRBmultireg}, typically created by \code{\link{FRBmultiregS}}, \code{\link{FRBmultiregMM}} or \code{\link{FRBmultiregGS}} }
  \item{confmethod}{ which kind of bootstrap confidence intervals to be displayed: 
'BCA'= bias corrected and accelerated method, 'basic'= basic bootstrap method, 'both'=both kinds of confidence intervals }
  \item{x}{ an \R object of class \code{summary.FRBmultireg}, resulting for example from \code{summary(\link{FRBmultiregS}(),...)} }
  \item{digits}{ number of digits for printing (defaulting to 3) }
  \item{\dots}{ potentially more arguments to be passed to methods }
}
\details{
The \code{print} method displays the components of the \code{summary} object, practically as listed in the Value section. 
}
\value{
  \code{summary} returns an object of class \code{summary.FRBmultireg}, which is a list containing:
  \item{responses}{the names of the response variables in the fitted model}
  \item{covariates}{the names of the covariates (predictors) in the fitted model}
  \item{Betawstd}{a data frame containing the coefficient estimates and their bootstrap standard errors}
  \item{Sigma}{estimate for the error covariance matrix}
  \item{limits.bca}{a list with for each response variable a matrix containing the lower and upper limits of the 
      BCa confidence intervals (only present when \code{confmethod="BCA"} or \code{confmethod="both"})}
  \item{limits.basic}{a list with for each response variable a matrix containing the lower and upper limits of the 
      basic bootstrap confidence intervals (only present when \code{confmethod="basic"} or \code{confmethod="both"})}
  \item{method }{multivariate regression method that was used}
  \item{conf}{confidence level that was used}
}
\author{ Gert Willems and Ella Roelant }
\seealso{ \code{\link{FRBmultiregS}}, \code{\link{FRBmultiregMM}}, \code{\link{FRBmultiregGS}},  \code{\link{print.FRBmultireg}}, \code{\link{plot.FRBmultireg}} }
\examples{

data(schooldata)
school.x <- data.matrix(schooldata[,1:5])
school.y <- data.matrix(schooldata[,6:8])

MMres <- FRBmultiregMM(school.x, school.y, R=999, conf = 0.99)
summary(MMres)  # -> print.summary.FRBmultireg() method

GSres <- FRBmultiregGS(school.x, school.y, bdp = 0.25)
summary(GSres, confmethod="both")  # -> print.summary.FRBmultireg() method

}
