\name{print.FRBhot}
\alias{print.FRBhot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Print Method for Objects of Class 'FRBhot' }
\description{
  This is the print method for objects of class \code{FRBhot} representing a robust Hotelling test.
}
\usage{
\method{print}{FRBhot}(x, digits = 3, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an \R object of class \code{FRBhot}, typically created by \code{\link{FRBhotellingS}} or \code{\link{FRBhotellingMM}} }
  \item{digits}{ number of digits for printing (defaulting to 3) }
  \item{\dots}{ potentially more arguments }
}
\details{
The \code{FRBhot} print function displays practically all components in the \code{FRBhot} object as returned by \code{\link{FRBhotellingS}}
or \code{\link{FRBhotellingMM}}. In particular, apart from the test statistic and bootstrap p-value, the output shows 
simultaneous confidence intervals for the (difference in) mean(s). The corresponding robust estimates of the (difference in) mean(s)
are also given, as well as the robust estimate of (common) covariance.     
}
\author{ Gert Willems and Ella Roelant }
\seealso{ \code{\link{FRBhotellingS}}, \code{\link{FRBhotellingMM}} }
\examples{

data(delivery)
delivery.x <- delivery[,1:2]
FRBhotellingMM(delivery.x)  # -> print.FRBhot() method

## Two sample robust Hotelling test
data(hemophilia)
grp <-as.factor(hemophilia[,3])
x <- hemophilia[which(grp==levels(grp)[1]),1:2]
y <- hemophilia[which(grp==levels(grp)[2]),1:2]
#using the pooled covariance matrix to estimate the common covariance matrix
FRBhotellingS(x,y,method="pool")
#using the estimator of He and Fung to estimate the common covariance matrix
FRBhotellingS(x,y,method="HeFung")

}
