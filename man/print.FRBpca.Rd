\name{print.FRBpca}
\alias{print.FRBpca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Print Method for Objects of Class 'FRBpca' }
\description{
  This is the print method for objects of class \code{FRBpca}.
}
\usage{
\method{print}{FRBpca}(x, digits=3, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an \R object of class \code{FRBpca}, typically created by \code{\link{FRBpcaS}} or \code{\link{FRBpcaMM}} }
  \item{digits}{ number of digits for printing (defaulting to 3) }
  \item{\dots}{ potentially more arguments }
}
\details{
  The print method for 'FRBpca' displays the estimated standard deviations of the principal components, with corresponding 
  bootstrap confidence intervals (BCa method). 
}
\references{ 
\itemize{
\item S. Van Aelst and G. Willems (2013). Fast and robust bootstrap for multivariate inference: The R package FRB. \emph{Journal of Statistical Software}, \bold{53}(3), 1--32.
        URL: \url{http://www.jstatsoft.org/v53/i03/}.
        }
}
\author{ Gert Willems and Ella Roelant }
\seealso{ \code{\link{FRBpcaS}}, \code{\link{FRBpcaMM}},  \code{\link{summary.FRBpca}} }
\examples{

data(ForgedBankNotes)

MMpcares <- FRBpcaMM(ForgedBankNotes, R=999, conf=0.95)
MMpcares # -> print.FRBpca() method
}