\name{delivery}
\alias{delivery}
\docType{data}
\title{ Delivery Time Data}
\description{
  Delivery Time Data, from Montgomery and Peck (1982). The aim is to explain the time required to service
a vending machine (Y) by means of the number of products stocked (X1) and the distance walked by the route driver (X2).
}
\usage{data(delivery)}
\format{
  A data frame with 25 observations on the following 3 variables.
\describe{
\item{\code{n.prod}}{Number of Products}
\item{\code{distance}}{Distance}
\item{\code{delTime}}{Delivery time} 
}
}
%\details{
%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
 Montgomery and Peck (1982, p.116)
}
\references{
P.J. Rousseeuw and A.M. Leroy (1987) Robust Regression and Outlier Detection, Wiley, page 155, table 23.
}
\examples{
data(delivery)
%## maybe str(delivery) ; plot(delivery) ...
}
\keyword{datasets}
