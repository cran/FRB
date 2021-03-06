\name{summary.FRBpca}
\alias{summary.FRBpca}
\alias{print.summary.FRBpca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Summary Method for Objects of Class 'FRBpca'  }
\description{
Summary method for objects of class \code{FRBpca}, and print method of the summary object.
}
\usage{
\method{summary}{FRBpca}(object, confmethod = c("BCA", "basic", "both"), digits = 3, ...)
\method{print}{summary.FRBpca}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ an \R object of class \code{FRBpca}, typically created by \code{\link{FRBpcaS}} or \code{\link{FRBpcaMM}} }
  \item{confmethod}{ which kind of bootstrap confidence intervals to be displayed: 
'BCA'= bias corrected and accelerated method, 'basic'= basic bootstrap method, 'both'= both kinds of confidence intervals }
  \item{digits}{ number of digits for printing (default is 3) }
  \item{x}{ an \R object of class \code{summary.FRBpca}, resulting from \code{summary(\link{FRBpcaS}(),...)} or \code{summary(\link{FRBpcaMM}(),...)} }
  \item{\dots}{ potentially more arguments to be passed to methods }
}
\details{
The \code{print} method displays mostly the components of the \code{summary} object as listed in the Value section. 
}
\value{
  \code{summary} returns an object of class \code{summary.FRBpca}, which contains the following components:
  \item{eigvals }{eigenvalues of the shape estimate (variances of the principal components) with confidence limits}
  \item{eigvecs }{eigenvectors of the shape estimate (loadings of the principal components)}
  \item{avgangle }{bootstrap estimates of average angles between true and estimated eigenvectors}
  \item{pvars }{cumulative percentage of variance explained by first principal components with confidence limits}
  \item{method }{PCA method that was used}
  \item{digits}{ number of digits for printing }
  }
\references{ 
\itemize{
\item S. Van Aelst and G. Willems (2013). Fast and robust bootstrap for multivariate inference: The R package FRB. \emph{Journal of Statistical Software}, \bold{53}(3), 1--32.
        URL: \url{http://www.jstatsoft.org/v53/i03/}.
        }
}
\author{ Gert Willems, Ella Roelant and Stefan Van Aelst}
\seealso{ \code{\link{FRBpcaS}}, \code{\link{FRBpcaMM}},  \code{\link{print.FRBpca}}, \code{\link{plot.FRBpca}} }
\examples{
data(ForgedBankNotes)

MMpcares <- FRBpcaMM(ForgedBankNotes, R=999, conf=0.95)
summary(MMpcares) # -> print.summary.FRBpca() method

}
