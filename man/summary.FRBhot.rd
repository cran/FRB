\name{summary.FRBhot}
\alias{summary.FRBhot}
\alias{print.summary.FRBhot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Summary Method for Objects of Class 'FRBhot'  }
\description{
Summary method for objects of class \code{FRBhot}, and print method of the summary object.
}
\usage{
\method{summary}{FRBhot}(object, digits = 5, ...)
\method{print}{summary.FRBhot}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ an \R object of class \code{FRBhot}, typically created by \code{\link{FRBhotellingS}} or \code{\link{FRBhotellingMM}} }
  \item{digits}{ number of digits for printing (default is 5) }
    \item{x}{ an \R object of class \code{summary.FRBhot}, resulting from \code{summary(\link{FRBhotellingS}(),...)} or \code{summary(\link{FRBhotellingMM}(),...)} }
    \item{\dots}{ potentially more arguments to be passed to methods }
}
\details{
The \code{print} method here displays the value of the test statistic and the corresponding bootstrap p-value. It also presents the simultaneous confidence intervals for the components of the location vector (or difference
between the two location vectors), and the robust estimates for the location vector(s) and covariance matrix.
}
\value{
  \code{summary.FRBhot} simply returns its two arguments in a list.
}
\references{ 
\itemize{
\item S. Van Aelst and G. Willems (2013). Fast and robust bootstrap for multivariate inference: The R package FRB. \emph{Journal of Statistical Software}, \bold{53}(3), 1--32.
        URL: \url{http://www.jstatsoft.org/v53/i03/}.
        }
}
\author{ Gert Willems, Ella Roelant and Stefan Van Aelst}
\seealso{ \code{\link{FRBhotellingS}}, \code{\link{FRBhotellingMM}},   \code{\link{plot.FRBhot}} }
\examples{
data(ForgedBankNotes)
samplemean <- apply(ForgedBankNotes, 2, mean)
res = FRBhotellingS(ForgedBankNotes, mu0=samplemean)

summary(res) # -> print.summary.FRBhot() method

}
