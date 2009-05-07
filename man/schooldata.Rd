\name{schooldata}
\alias{schooldata}
\docType{data}
\title{ School Data}
\description{
 School Data, from Charnes et al. (1981). The aim is to explain scores on 3 
different tests from 70 school sites by means of 5 explanatory variables.
}
\usage{data(schooldata)}
\format{
  A data frame with 70 observations on the following 8 variables.

\describe{
\item{\code{education}}{ education level of mother as measured in terms of percentage of high school 
graduates among female parents}
\item{\code{occupation}}{ highest occupation of a family member according to a pre-arranged rating scale}
\item{\code{visit}}{ parental visits index representing the number of visits to the school site}
\item{\code{counseling}}{ parent counseling index calculated from data on time spent with child on 
school-related topics such as reading together, etc.}
\item{\code{teacher}}{ number of teachers at a given site}
\item{\code{reading}}{ total reading score as measured by the Metropolitan Achievement Test}
\item{\code{mathematics}}{ total mathematics score as measured by the Metropolitan Achievement Test}
\item{\code{selfesteem}}{ Coopersmith Self-Esteem Inventory, intended as a measure of self-esteem}
} 
}
%\details{
%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
 Charnes et al. (1981)
}
\references{
\itemize{
 \item A. Charnes, W.W. Cooper and E. Rhodes (1981) Evaluating Program and Managerial Efficiency: An Application of Data
Envelopment Analysis to Program Follow Through. \emph{Management Science}, \bold{27}, 668-697.
}
}
\examples{
data(schooldata)
%## maybe str(schooldata) ; plot(schooldata) ...
}
\keyword{datasets}
