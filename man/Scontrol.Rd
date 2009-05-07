\name{Scontrol}
\alias{Scontrol}
\alias{MMcontrol}
\alias{GScontrol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Tuning parameters for multivariate S, MM and GS estimates }
\description{
  Tuning parameters for multivariate S, MM and GS estimates as used in FRB functions for multivariate regression, 
  PCA and Hotelling tests. Mainly regarding the fast-(G)S algorithm.
}
\usage{
Scontrol(nsamp = 500, k = 3, bestr = 5, convTol = 1e-10, maxIt = 50)

MMcontrol(bdp = 0.5, eff = 0.95, shapeEff = FALSE, convTol.MM = 1e-07, 
          maxIt.MM = 50, fastScontrols = Scontrol(...), ...)

GScontrol(nsamp = 100, k = 3, bestr = 5, convTol = 1e-10, maxIt = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nsamp}{ number of random subsamples to be used in the fast-(G)S algorithm }
  \item{k}{ number of initial concentration steps performed on each subsample candidate }
  \item{bestr}{ number of best candidates to keep for full iteration (i.e. concentration steps until convergence) }
  \item{convTol}{ relative convergence tolerance for estimates used in (G)S-concentration iteration}
  \item{maxIt}{ maximal number of steps in (G)S-concentration iteration }
  \item{bdp}{ breakdown point of the MM-estimates; usually equals 0.5 }
  \item{eff}{ Gaussian efficiency of the MM-estimates; usually set at 0.95 }
  \item{shapeEff}{ logical; if \code{TRUE}, \code{eff} is with regard to shape-efficiency, 
                                           otherwise location-efficiency  }
  \item{convTol.MM}{ relative convergence tolerance for estimates used in MM-iteration }
  \item{maxIt.MM}{ maximal number of steps in MM-iteration }
  \item{fastScontrols}{ the tuning parameters of the initial S-estimate }
  \item{\dots}{ allows for any individual parameter from \code{Scontrol} to be set directly }
}
\details{
The default number of random samples is lower for GS-estimates than for S-estimates, because computations regarding the former 
are more demanding.

}
\value{
A list with the tuning parameters as set by the arguments.
}
\author{ Gert Willems and Ella Roelant }

\seealso{ \code{\link{Sest_loccov}}, \code{\link{MMest_loccov}},  \code{\link{GSest_multireg}},  \code{\link{Sest_multireg}}, 
\code{\link{GSest_multireg}},  \code{\link{twosampleS}},  \code{\link{twosampleMM}}, \code{\link{FRBpcaS}}, ... }
\examples{
## Show the default settings:
str(Scontrol())
str(MMcontrol())
str(GScontrol())

}
