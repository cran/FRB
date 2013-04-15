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
  \item{x}{ an \R object of class \code{FRBmultireg}, typically created by \code{\link{FRBmultiregS}}, \code{\link{FRBmultiregMM}} or \code{\link{FRBmultiregGS}} or by \code{\link{GSest_multireg}}, \code{\link{MMest_multireg}} or \code{\link{Sest_multireg}} }
  \item{digits}{ number of digits for printing (defaulting to 3) }
  \item{\dots}{ potentially more arguments }
}
\details{
  The print method for \code{FRBmultireg} displays the estimated coefficients of the multivariate linear regression model. 
}
\references{ 
\itemize{
\item S. Van Aelst and G. Willems (2013). Fast and robust bootstrap for multivariate inference: The R package FRB. \emph{Journal of Statistical Software}, \bold{53}(3), 1--32.
        URL: \url{http://www.jstatsoft.org/v53/i03/}.
        }
}
\author{ Gert Willems, Ella Roelant and Stefan Van Aelst}
\seealso{\code{\link{summary.FRBmultireg}}, \code{\link{FRBmultiregS}}, \code{\link{FRBmultiregMM}}, \code{\link{FRBmultiregGS}}, \code{\link{GSest_multireg}}, \code{\link{MMest_multireg}}, \code{\link{Sest_multireg}} }
\examples{
data(schooldata)
MMres <- FRBmultiregMM(cbind(reading,mathematics,selfesteem)~., data=schooldata)

MMres # -> print.FRBmultireg() method

}
