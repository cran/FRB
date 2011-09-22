\name{diagplot}
\alias{diagplot}
\alias{diagplot.FRBmultireg}
\alias{diagplot.FRBpca}
\alias{diagplot.FRBhot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot Method for Objects of class 'FRBmultireg' }
\description{
  Diagnostic plots for objects of class \code{FRBmultireg}, \code{FRBpca} and \code{FRBhot}. It shows robust distances
  and allows detection of multivariate outliers.
}
\usage{
\method{diagplot}{FRBmultireg}(x, Xdist = TRUE, ...)

\method{diagplot}{FRBpca}(x, EIF = TRUE, ...)

\method{diagplot}{FRBhot}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an \R object of class \code{FRBmultireg} (typically created by  \code{\link{FRBmultiregS}}, \code{\link{FRBmultiregMM}} or \code{\link{FRBmultiregGS}} or by \code{\link{Sest_multireg}}, \code{\link{MMest_multireg}} or \code{\link{GSest_multireg}})
   or an \R object of class \code{FRBpca} (typically created by  \code{\link{FRBpcaS}} or \code{\link{FRBpcaMM}})
   or an \R object of class \code{FRBhot} (typically created by  \code{\link{FRBhotellingS}} or \code{\link{FRBhotellingMM}}) }
  \item{Xdist}{ logical: if TRUE, the plot shows the robust distance versus the distance in the space of the explanatory variables;
      if FALSE, it plots the robust distance versus the index of the observation}
  \item{EIF}{ logical: if TRUE, the plot shows the robust distance versus an influence measure for each point;
        if FALSE, it plots the robust distance versus the index of the observation }
  \item{\dots}{ potentially more arguments to be passed }
}
\details{
The diagnostic plots are based on the robust distances of the observations. In a multivariate sample \eqn{X_n=\{\mathbf{x}_1,...,\mathbf{x}_n\}}{X_n=\{x_1,...,x_n\}},
the robust distance \eqn{d_i} of observation \eqn{i} is given by
\eqn{d_i^2=(\mathbf{x}_i-\hat{\mu})'\hat{\Sigma}^{-1}(\mathbf{x}_i-\hat{\mu})}{d_i^2=(x_i-\mu)'\Sigma^(-1)(x_i-\mu)}.
where \eqn{\hat{\mu}}{\mu} and \eqn{\hat{\Sigma}}{\Sigma} are robust estimates of location and covariance.
Observations with large robust distance are considered as outlying.

The default diagnostic plot in the multivariate regresssion setting (i.e. for objects of type \code{FRBmultireg} and \code{Xdist=TRUE}),
shows the residual distances (i.e. the robust distances of the multivariate residuals) based on the estimates in \code{x},
versus the distances within the space of the explanatory variables. The latter are based on robust estimates of location and scatter for the
data matrix \code{x$X} (without intercept). Computing these robust estimates may take an appreciable amount of time. The estimator used
corresponds to the one which was used in obtaining \code{Xmultireg} (with the same breakdown point, for example, and the same control parameters).
On the vertical axis a cutoff line is drawn at the square root of the .975 quantile of the chi-squared distribution with degrees of
freedom equal to the number of response variables. On the horizontal axis the same quantile is drawn but now with degrees of freedom
equal to the number of covariates (not including intercept).
Those points to the right of the cutoff can be viewed as high-leverage points. These can be classified into so-called
'bad' or 'good' leverage points depending on whether they are above or below the cutoff. Points above the cutoff but to the
left of the vertical cutoff are sometimes called vertical outliers.
See also Van Aelst and Willems (2005) for example.

To avoid the additional computation time, one can choose \code{Xdist=FALSE}, in which case the residual distances are simply plotted
versus the index of the observation.

The default plot in the context of PCA (i.e. for objects of type \code{FRBpca} and \code{EIF=FALSE})
is a plot proposed by Pison and Van Aelst (2004). It shows the robust distance versus a measure of the overall empirical influence
of the observation on the (classical) principal components. The empirical influences are obtained by using the influence function of
the eigenvectors of the empirical or classical shape estimator at the normal model, and by
substituting therein the robust estimates for the population parameters.
The overall influence value is then defined by averaging the squared influence
over all coefficients in the eigenvectors.
The vertical line on the plot is an indicative cutoff value, obtained through simulation. This last part takes
a few moments of computation time.

Again, to avoid the additional computation time, one can choose \code{EIF=FALSE}, in which case the robust distances are simply plotted
versus the index of the observation.

For the result of the robust Hotelling test (i.e. for objects of type \code{FRBhot}), the method plots the robust
distance versus the index. In case of a two-sample test, the indices are within-sample and a vertical line separates
the two groups. In the two-sample case, each group has its own location estimate \eqn{\hat{\mu}}{\mu} and a common
covariance estimate \eqn{\hat{\Sigma}}{\Sigma}.
}
%\value{
%  Produces one or more arrays of histograms.
%}
\references{
\itemize{
\item G. Pison and S. Van Aelst (2004). Diagnostic Plots for Robust Multivariate Methods. \emph{Journal of Computational and Graphical Statistics}, \bold{13}, 310--329.
%\item G. Pison and S. Van Aelst (2002). Analyzing robust multivariate methods
%with a diagnostic plot.  In \emph{Proceedings in Computational Statistics 2002 (W.
%Hardle and B. Ronz, eds.)}, 165-170.
\item  S. Van Aelst and G. Willems (2005). Multivariate regression S-estimators for robust estimation and
inference. \emph{Statistica Sinica}, \bold{15}, 981--1001.
}
}
\author{ Gert Willems and Ella Roelant }
%\note{ ~~further notes~~
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{FRBmultiregS}}, \code{\link{FRBmultiregMM}}, \code{\link{FRBmultiregGS}}, \code{\link{FRBpcaS}}
, \code{\link{FRBpcaMM}}, \code{\link{FRBhotellingS}}, \code{\link{FRBhotellingMM}} }
\examples{

# for multivariate regression:
data(schooldata)
MMres <- MMest_multireg(cbind(reading,mathematics,selfesteem)~., data=schooldata)
diagplot(MMres)
# a large 'bad leverage' outlier should be noticeable (observation59)

# for PCA:
data(ForgedBankNotes)
MMres <- FRBpcaMM(ForgedBankNotes, R=10)
diagplot(MMres)
# a group of 15 fairly strong outliers can be seen which apparently would have
# a large general influence on a classical PCA analysis

# for Hotelling tests (two-sample)
data(hemophilia)
MMres <- FRBhotellingMM(cbind(AHFactivity,AHFantigen)~gr,data=hemophilia, R=10)
diagplot(MMres)
# the data seem practically outlier-free


}
