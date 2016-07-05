\name{GEneric regression models}
\alias{reg.fit}
\title{
Regression modelling
}

\description{
Generic regression modelling function.
}

\usage{
reg.fit(y, dataset, event = NULL, reps = NULL, group = NULL, slopes = FALSE, 
reml = FALSE, model = NULL, robust = FALSE, xnew = NULL)
}

\arguments{
\item{y}{
The target (dependent) variable. It can be a numerical variable, factor, ordinal factor, percentages, matrix, or time to event. 
If the values are proportions or percentages, i.e. strictly within 0 and 1 they are mapped into R using log( target/(1 - target) ).
If they are compositional data the additive log-ratio (multivariate logit) transformation is aplied beforehand. 
}
\item{dataset}{
The indendent variable(s). It can be a vector, a matrix or a dataframe with continuous only variables, a data frame with mixed or only categorical variables.
}
\item{event}{
This is NULL unless you have time to event data (survival regression).
}
\item{reps}{
This is NULL unless you have time measurements (longitudinal data).
}
\item{group}{
This is NULL unless you have grouped (or clustered) data or longitudinal data (is the latter case the arugment reps is required). 
}
\item{slopes}{
This is for the longitudinal data only, TRUE or FALSE. Should random slopes be added or not?
}
\item{reml}{
This is for the longitudinal or grouped data only, TRUE or FALSE. If TRUE, REML will be used, otherwise ML will be used.
}
\item{model}{
The type of model you want to use. It can be specified by the user or left NULL, if other correct arguments have been passed. 
Poissible values (apart from NULL) are: "gaussian" (default), "binary", "multinomial", "poisson", "ordinal", "Cox", "weibull", "zip0", "zipx", "beta", "median", "negbin", "longitudinal" and "grouped".  The "zip0" means that the zero part is constant, the variables are not associated with the excessive zeros. This is done with "zipx". The value "grouped" refers to grouped data, but this does not have to be given if the argument "group" is given, but not the argument "reps.
}
\item{robust}{
A boolean variable. If TRUE robust models will be returned. Currently this is supported by models "gaussian", "poisson" and "binary".
}
\item{xnew}{
If you have new data whose target values you want to predict put it here, otherwise leave it blank.
}
}

\details{
This is a generic regression function, which offers prediction as well. It is important that you pass the arguments with their names, for example if you have time to event data, write "event = ..." and not just put your event variable. This will avoid confusion. 
}

\value{
A list including:
\item{mod}{
The fitted model.
}
\item{pred}{
If you have new data the predicted values of the target (dependent) variable.
}
}

\references{
Almost the same as in \code{\link{CondIndTests}}.
}

\author{
 Michail Tsagris

R implementation and documentation: Giorgos Athineou <athineou@csd.uoc.gr> and Michail Tsagris <mtsagris@yahoo.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{ridge.reg}, \link{ses.model}, \link{mmpc.model} }
}

\examples{
#simulate a dataset with continuous data
dataset <- matrix(runif(100 * 10, 1, 100), nrow = 100 ) 
#the target feature is the last column of the dataset as a vector
target <- dataset[, 1]
dataset <- dataset[, -1]
a <- reg.fit(target, dataset)
}

\keyword{ Regression model }