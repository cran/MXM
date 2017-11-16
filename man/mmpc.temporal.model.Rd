\name{Generalised linear mixed models based on temporal SES and MMPC outputs}
\alias{mmpc.temporal.model}
\title{
Generalised linear mixed model(s) based obtained from temporal SES or MMPC
}

\description{
One or more regression models obtained from SES or MMPC, are returned. 
}

\usage{
mmpc.temporal.model(target, dataset, reps = NULL, group, slopes = FALSE, wei = NULL,
mmpctemporal.Object, test = NULL) 
}

\arguments{
\item{target}{
The class variable. Provide a vector with continuous (normal), binary (binomial) or discrete (Poisson) data.
}
\item{dataset}{
The dataset; provide either a data frame or a matrix (columns = variables, rows = observations). Currently, only continuous datasets are supported. 
}
\item{reps}{
A numeric vector containing the time points of the subjects. Its length is equal to the length of the target variable. If you have clustered data, leave this NULL. 
}
\item{group}{
A numeric vector containing the subjects or groups. It must be of the same legnth as target. 
}
\item{slopes}{
Should random slopes for the ime effect be fitted as well? Default value is FALSE. 
}
\item{wei}{
A vector of weights to be used for weighted regression. The default value is NULL. 
}
\item{mmpctemporal.Object}{
An object with the results of a wald.ses run.
}
\item{test}{
The conditional independence test to use. Default value is NULL. Currently, the only available conditional independence test is the \code{\link{testIndGLMM}}, which fits linear mixed models. 
}
}

\details{
 This command is useful if you want to see all models and check for example their fitting ability, MSE in linear models for exmaple.
}

\value{
A list including:
\item{mod}{
Depending on the number of signatures requested, one or models will be returned.
}
\item{ypografi}{
A matrix (or just one vector if one signature only) with the variables of each signature, along with the BIC of the corresponding regression model.
}
}
\references{
I. Tsamardinos, M. Tsagris and V. Lagani (2015). Feature selection for longitudinal data. Proceedings of the 10th conference of the Hellenic Society for Computational Biology & Bioinformatics (HSCBB15)

Eugene Demidenko (2013). Mixed Models: Theory and Applications with R, 2nd Edition. New Jersey: Wiley \& Sons.

Pinheiro J. and D. Bates. Mixed-effects models in S and S-PLUS. Springer Science \& Business Media, 2006.
}

\author{
 Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@csd.uoc.gr>
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{ \link{SES}, \link{MMPC}, \link{cv.ses}, \link{cv.mmpc} }
}

\examples{
\dontrun{
require(lme4)
data(sleepstudy)
attach(sleepstudy)
x <- matrix(rnorm(180 * 100),ncol = 100) ## unrelated predictor variables
m1 <- SES.temporal(Reaction, Days, Subject, x)
m2 <- MMPC.temporal(Reaction, Days, Subject, x)
mod <- mmpc.temporal.model(target, dataset = x, reps = Days, group = Subject, slopes = FALSE, 
mmpctemporal.Object = m2) 
}
}



\keyword{ Regression modelling }