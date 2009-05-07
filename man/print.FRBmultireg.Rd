\name{print.FRBmultireg}
\alias{print.FRBmultireg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Print Method for Objects of Class 'FRBmultireg'  }
\description{
  This is the print method for objects of class \code{FRBmultireg}.
}
\usage{
\method{print}{FRBmultireg}(x, digits = 3, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an \R object of class \code{FRBmultireg}, typically created by \code{\link{FRBmultiregS}}, \code{\link{FRBmultiregMM}} or \code{\link{FRBmultiregGS}} }
  \item{digits}{ number of digits for printing (defaulting to 3) }
  \item{\dots}{ potentially more arguments }
}
\details{
  The print method for \code{FRBmultireg} displays the estimated coefficients of the multivariate linear regression model, 
  with corresponding standard errors based on Fast and Robust Bootstrap. 
}
\author{ Gert Willems and Ella Roelant }
\seealso{\code{\link{summary.FRBmultireg}}, \code{\link{FRBmultiregS}}, \code{\link{FRBmultiregMM}}, \code{\link{FRBmultiregGS}} }
\examples{
data(schooldata)
school.x <- data.matrix(schooldata[,1:5])
school.y <- data.matrix(schooldata[,6:8])

MMres <- FRBmultiregMM(school.x, school.y, R=999, conf = 0.99)

MMres # -> print.FRBmultireg() method

}
