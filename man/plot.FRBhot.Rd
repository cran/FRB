\name{plot.FRBhot}
\alias{plot.FRBhot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot Method for Objects of class 'FRBhot' }
\description{
Plot function for \code{FRBhot} objects: plots the bootstrap histogram of the null distribution, and
the simultaneous confidence limits (scaled)
}
\usage{
\method{plot}{FRBhot}(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an \R object of class \code{FRBhot}, typically created by \code{\link{FRBhotellingS}} or \code{\link{FRBhotellingMM}} }
  \item{...}{ potentially more arguments } 
  }
\details{
This generic \code{plot} function presents two graphs. The first (top panel) is a histogram representing the test statistics
in the bootstrap samples, which estimate the null distribution.
A red line indicates the test statistic in the original sample (but is not shown when this value exceeds 100).

The second (bottom panel) displays the simultaneous confidence intervals based on the same bootstrap result.
The intervals are scaled such that they all have the same length. Furthermore, in case of the one-sample test
the intervals are shown relative to the hypothesized value \code{mu0}. Such visualization is meant to easily
recognize the extent to which each variable is responsible for the overall deviation from the hypothesized value.


}
\references{ 
\itemize{
\item S. Van Aelst and G. Willems (2013). Fast and robust bootstrap for multivariate inference: The R package FRB. \emph{Journal of Statistical Software}, \bold{53}(3), 1--32.
        URL: \url{http://www.jstatsoft.org/v53/i03/}.
        }
}
\author{ Gert Willems, Ella Roelant and Stefan Van Aelst }
\seealso{ \code{\link{FRBhotellingS}}, \code{\link{FRBhotellingMM}} }
\examples{

## One sample robust Hotelling test
data(ForgedBankNotes)
samplemean <- apply(ForgedBankNotes, 2, mean)
res = FRBhotellingS(ForgedBankNotes, mu0=samplemean,R=99)

plot(res)

# Note that the test rejects the hypothesis that the true mean equals the
# sample mean; this is due to outliers in the data (i.e. the robustly estimated
# center apparently significantly differs from the non-robust sample mean.

# It is clear from the scaled simultaneous confidence limits that the rejection
# of the hypothesis is due to the differences in variables Bottom and Diagonal

# For comparison, the hypothesis would be accepted if only the first three
# variables were considered:
res = FRBhotellingS(ForgedBankNotes[,1:3], mu0=samplemean[1:3],R=99)
plot(res)

## Two sample robust Hotelling test
data(hemophilia)
res <- FRBhotellingMM(cbind(AHFactivity,AHFantigen)~gr,data=hemophilia,R=99)
plot(res)

# From the confidence limits it can be seen that the significant difference
# is mainly caused by the AHFactivity variable.
# the red line on the histogram indicates the test statistic value in the original
# sample (it is omitted if the statistic exceeds 100)

}
