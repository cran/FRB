\name{FRBhotellingS}
\alias{FRBhotellingS}
\alias{FRBhotellingS.default}
\alias{FRBhotellingS.formula}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Robust Hotelling test using the S-estimator }
\description{
Robust one-sample and two-sample Hotelling test using the S-estimator and the Fast and Robust Bootstrap.
}
\usage{
\method{FRBhotellingS}{formula}(formula, data=NULL, ...)

\method{FRBhotellingS}{default}(X, Y=NULL, mu0 = 0, R = 999, bdp = 0.5, conf = 0.95,
method = c("HeFung", "pool"), control=Scontrol(...),
na.action=na.omit, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{ an object of class \code{\link{formula}}; a symbolic description of the model to be fit.}
  \item{data}{ data frame from which variables specified in formula are to be taken.}
  \item{X}{ a matrix or data-frame}
  \item{Y}{ an optional matrix or data-frame in case of a two-sample test}
  \item{mu0}{ an optional vector of data values (or a single number which will be repeated p times) indicating the
true value of the mean (does not apply in case of the two-sample test). Default is the null vector \code{mu0=0}}
  \item{R}{ number of bootstrap samples. Default is \code{R=999}.  }
  \item{bdp}{ required breakdown point. Should have \eqn{0 < } \code{bdp} \eqn{\le 0.5}, the default is 0.5 }
  \item{conf}{ confidence level for the simultaneous confidence intervals. Default is \code{conf=0.95} }
  \item{method}{ for the two-sample Hotelling test, indicates the way the common covariance matrix is estimated:
    \code{"pool"}= pooled covariance matrix, \code{"HeFung"}= using the He and Fung method }
  \item{control}{a list with control parameters for tuning the computing algorithm, see \code{\link{Scontrol}}().}
    \item{na.action}{ a function which indicates what should happen when the data contain NAs. Defaults to \code{\link{na.omit}}.}
  \item{...}{allows for specifying control parameters directly instead of via \code{control}}
 }

\details{  The classical Hotelling test for testing if the mean equals a certain center or if two means 
are equal is modified into a robust one through substitution of the empirical estimates 
by the S-estimates of location and scatter. The S-estimator uses Tukey's biweight function where the constant is chosen
to obtain the desired breakdown point as specified by \code{bdp}. 
One-sample S-estimates are computed by a call to the implementation of the fast-S algorithm in the \pkg{rrcov} package of Todorov and Filzmoser (2009). For two-sample S-estimates an adaptation of the fast-S algorithm is used. 
The tuning parameters of the algorithm can be changed via \code{control}.

The fast and robust bootstrap is used to mimic the distribution of the test statistic under the null hypothesis. 
For instance, the 5\% critical value for the test is given by the 95\% quantile of the recalculated statistics. 

Robust simultaneous confidence intervals for linear combinations of the mean (or difference in means) are developed similarly to the classical case
(Johnson and Wichern, 1988, page 239). The value \code{CI} is a matrix with the confidence intervals for each element
of the mean (or difference in means), with level \code{conf}. It consists of two rows, the first being the lower bound and the second the upper bound.
Note that these intervals are rather conservative in the sense that the simultaneous confidence level holds for all linear combinations
and here only \eqn{p} of these are considered (with \eqn{p} the dimension of the data).

For the two-sample Hotelling test we assume that the samples have an underlying distribution with the same covariance matrix.
This covariance matrix can be estimated in two different ways using the pooled covariance matrix or the two-sample 
estimator of He and Fung (He and Fung 2000), and argument \code{method} defaults to the second option. 
For more details see Roelant et al. (2008).

In the two-sample version, the null hypothesis always states that the two means are equal. For the one-sample version, the default
null hypothesis is that the mean equals zero, but the hypothesized value can be changed and specified through argument \code{mu0}.
 
Bootstrap samples are discarded if the fast and robust covariance estimate is not positive definite, such that the actual number 
of recalculations used can be lower than \code{R}. This number is returned as \code{ROK}.
 
%See \code{\link{print.FRBhot}} for details on the output.              
}
\value{
An object of class \code{FRBhot} which extends class \code{htest} and contains at least the following components:
\item{statistic }{the value of the robust test statistic.}
\item{pvalue }{p-value of the robust one or two-sample Hotelling test, determined by the fast and robust bootstrap}
\item{estimate }{the estimated mean vector or vectors depending on whether it was a one-sample test or a two-sample test.} 
\item{alternative }{a character string describing the alternative hypothesis.}
\item{method }{a character string indicating what type of Hotelling test was performed.}
\item{data.name }{a character string giving the name(s) of the data.}
\item{teststat.boot }{the bootstrap recalculated values of the robust test statistic.}
\item{CI }{bootstrap simultaneous confidence intervals for each component of the center}
\item{conf }{a copy of the \code{conf} argument}
\item{Sigma }{covariance of one-sample or common covariance matrix in the case of two samples}
\item{w}{implicit weights corresponding to the S-estimates (i.e. final weights in the RWLS procedure at the end of the fast-S algorithm)}
\item{outFlag}{outlier flags: 1 if the robust distance of the observation exceeds the .975 quantile of (the square root of)
 the chi-square distribution with degrees of freedom equal to the dimension of \code{X}; 0 otherwise}
\item{ROK}{number of bootstrap samples actually used (i.e. not discarded due to non-positive definite covariance }
%\item{Mu }{center of the sample in case of one-sample Hotelling test}
%\item{Mu1 }{center of the first sample in case of the two-sample Hotelling test}
%\item{Mu2 }{center of the second sample in case of the two-sample Hotelling test}
%\item{X, Y}{copies of the \code{Xdata} and \code{Ydata} arguments as matrices}
}
\references{ 
\itemize{
\item X. He and W.K. Fung (2000) High breakdown estimation for multiple populations with 
applications to discriminant analysis. \emph{Journal of Multivariate Analysis}, \bold{72}, 151--162.
\item R.A. Johnson, D.W. Wichern (1988) Applied Multivariate Statistical Analysis, 
2nd Edition, Prentice-Hall.
\item E. Roelant, S. Van Aelst and G. Willems, (2008) Fast Bootstrap for Robust Hotelling Tests, COMPSTAT 2008: 
Proceedings in Computational Statistics (P. Brito, Ed.) Heidelberg: Physika-Verlag, 709--719.
\item M. Salibian-Barrera, S. Van Aelst and G. Willems (2008) Fast and robust 
bootstrap. \emph{Statistical Methods and Applications}, \bold{17}, 41--71. 
\item V. Todorov and P. Filzmoser (2009) An Object-Oriented Framework for Robust Multivariate Analysis. \emph{Journal of Statistical Software}, \bold{32}, 1--47. URL \url{http://www.jstatsoft.org/v32/i03/}.
}
}
\author{ Ella Roelant, Stefan Van Aelst and Gert Willems }
%\note{ ~~further notes~~ 

% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{plot.FRBhot}},\code{\link{summary.FRBhot}}, \code{\link{FRBhotellingMM}}, \code{\link{Scontrol}} }
\examples{
## One sample robust Hotelling test
data(delivery)
delivery.x <- delivery[,1:2]
FRBhotellingS(delivery.x)

## One sample robust Hotelling test
data(ForgedBankNotes)
samplemean <- apply(ForgedBankNotes, 2, mean)
res = FRBhotellingS(ForgedBankNotes, mu0=samplemean)
res
# Note that the test rejects the hypothesis that the true mean equals the
# sample mean; this is due to outliers in the data (i.e. the robustly estimated
# mean apparently significantly differs from the non-robust sample mean.

# Graphical display of the results:
plot(res)
# It is clear from the (scaled) simultaneous confidence limits that the rejection
# of the hypothesis is due to the differences in variables Bottom and Diagonal

%# For comparison, the hypothesis would be accepted if only the first three
%# variables were considered:
%\dontrun{
%res = FRBhotellingS(ForgedBankNotes[,1:3], mu0=samplemean[1:3])
%plot(res)}

## Two sample robust Hotelling test
data(hemophilia)
grp <-as.factor(hemophilia[,3])
x <- hemophilia[which(grp==levels(grp)[1]),1:2]
y <- hemophilia[which(grp==levels(grp)[2]),1:2]

#using the pooled covariance matrix to estimate the common covariance matrix
res = FRBhotellingS(x,y,method="pool")

#using the estimator of He and Fung to estimate the common covariance matrix
res = FRBhotellingS(x,y,method="HeFung")

# or using the formula interface
res = FRBhotellingS(as.matrix(hemophilia[,-3])~hemophilia[,3],method="HeFung")

# From the confidence limits it can be seen that the significant difference
# is mainly caused by the AHFactivity variable. The graphical display helps too:
plot(res)
# the red line on the histogram indicates the test statistic value in the original
# sample (it is omitted if the statistic exceeds 100)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
\keyword{ robust }
