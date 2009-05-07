\name{FRBhotellingS}
\alias{FRBhotellingS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Robust Hotelling test using the S-estimator }
\description{
Robust one-sample and two-sample Hotelling test using the S-estimator and the Fast and Robust Bootstrap.
}
\usage{
FRBhotellingS(Xdata, Ydata=NULL, mu0 = 0, R = 999, bdp = 0.5, conf = 0.95, 
                method = c("pool", "HeFung"), control=Scontrol(...), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Xdata}{ a matrix or data-frame}
  \item{Ydata}{ an optional matrix or data-frame in case of a two-sample test}
  \item{mu0}{ an optional vector of data values (or a single number which will be repeated p times) indicating the
true value of the mean (does not apply in case of the two-sample test). Default is the null vector \code{mu0=0}}
  \item{R}{ number of bootstrap samples. Default is \code{R=999}  }
  \item{bdp}{ required breakdown point. Should have \eqn{0 < } \code{bdp} \eqn{\le 0.5}, the default is 0.5 }
  \item{conf}{ confidence level for the simultaneous confidence intervals. Default is \code{conf=0.95} }
  \item{method}{ for the two-sample Hotelling test, indicates the way the common covariance matrix is estimated:
    \code{"pool"}= pooled covariance matrix, \code{"HeFung"}= using the He and Fung method }
  \item{control}{a list with control parameters for tuning the computing algorithm, see \code{\link{Scontrol}}().}
  \item{...}{allows for specifying control parameters directly instead of via \code{control}}
 }

\details{  The classical Hotelling test for testing if the mean equals a certain center or if two means 
are equal is modified into a robust one through substitution of the empirical estimates 
by the S-estimates of location and scatter. The S-estimator uses Tukey's biweight function where the constant is chosen
to obtain the desired breakdown point as specified by \code{bdp}. 
The S-estimator is computed by a call to \code{\link{Sest_loccov}}() or \code{\link{twosampleS}}(), in case of one- and two-sample 
tests respectively. These functions implement a fast-S-type algorithm, the tuning parameters of which can be changed via \code{control}.

The fast and robust bootstrap is used to mimic the distribution of the test statistic under the null hypothesis. 
For instance, the 5\% critical value for the test is given by the 95\% quantile of the recalculated statistics. 

Robust simultaneous confidence intervals for linear combinations of the center are developed similarly to the classical case 
(Johnson and Wichern, 1988, page 239). The value \code{CI} is a (p x 2) matrix with the confidence intervals for each element 
of the center. It consists of two rows, the first being the lower bound and the second the upper bound. 

For the two-sample Hotelling test we assume that the samples have an underlying distribution with the same covariance matrix.
This covariance matrix can be estimated in two different ways using the pooled covariance matrix or the two-sample 
estimator of He and Fung (He and Fung 2000), and argument \code{method} defaults to the first option. 
For more details see Roelant et al. (2008).

In the two-sample version, the null hypothesis always states that the two means are equal. For the one-sample version, the default
null hypothesis is that the mean equals zero, but the hypothesized value can be changed and specified through argument \code{mu0}.
 
See \code{\link{print.FRBhot}} for details on the output.              
}
\value{
An object of class \code{FRBhot}, which is a list containing the following components:
\item{pvalue }{p-value of the robust one or two-sample Hotelling test, determined by the fast and robust bootstrap}
\item{teststat }{the value of the robust test statistic.}
\item{teststat.boot }{the bootstrap recalculated values of the robust test statistic.}
\item{Mu }{center of the sample in case of one-sample Hotelling test}
\item{Mu1 }{center of the first sample in case of the two-sample Hotelling test}
\item{Mu2 }{center of the second sample in case of the two-sample Hotelling test}
\item{Sigma }{covariance of one sample or common covariance matrix in the case of two samples}
\item{CI }{bootstrap simultaneous confidence intervals for each component of the center}
\item{conf }{a copy of the \code{conf} argument}
\item{data }{the names of the \code{Xdata} and possibly \code{Ydata} object}
\item{meth }{a character string giving the estimator that was used}
}
\references{ 
\itemize{
\item X. He and W.K. Fung (2000) High breakdown estimation for multiple populations with 
applications to discriminant analysis. \emph{Journal of Multivariate Analysis}, \bold{72}, 151-162.
\item R.A. Johnson, D.W. Wichern (1988). Applied Multivariate Statistical Analysis, 
2nd Edition, Prentice-Hall.
\item E. Roelant, S. Van Aelst and G. Willems, (2008) Fast Bootstrap for Robust Hotelling Tests, COMPSTAT 2008: 
Proceedings in Computational Statistics (P. Brito, Ed.) Heidelberg: Physika-Verlag, 709-719.
\item M. Salibian-Barrera, S. Van Aelst and G. Willems (2008) Fast and robust 
bootstrap. \emph{Statistical Methods and Applications}, \bold{17}, 41-71. 
}
}
\author{ Ella Roelant and Gert Willems }
%\note{ ~~further notes~~ 

% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{print.FRBhot}}, \code{\link{FRBhotellingMM}}, \code{\link{Scontrol}} }
\examples{
## One sample robust Hotelling test
data(delivery)
delivery.x <- delivery[,1:2]
FRBhotellingS(delivery.x)

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

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
\keyword{ robust }
