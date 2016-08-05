\name{ROC and area under the curve}
\alias{auc}
\title{
ROC and area under the curve
}

\description{
Receiver operating curve and area under the curve.
}

\usage{
auc(group, preds, roc = FALSE, cutoffs = NULL) 
}

\arguments{
\item{group}{
The true labels, either a factor or a numerical vector with two numbers only. 
}
\item{preds}{
The predicted values of each group. 
}
\item{roc}{
If you want the ROC to appear set it to TRUE.
}
\item{cutoffs}{
If you provide a vector with decreasing numbers from 1 to 0 that will be used for the ROC, otherwise, the values from 1 to 0 with a step equal to -0.01 will be used.
}
}

\details{
The ara under the curve is returned. The user has the option of getting the receiver operating curve as well.
}

\value{
If the roc is set to FALSE, the area under the curve is returned only. Otherwise a list including
\item{cutoffs}{
The cutoff values.
}
\item{sensitivity}{
The sensitivity values for each cutoff value.
}
\item{specificity}{
The specificity value for each cutoff value.
}
\item{youden}{
The Youden index which is defined as the maximum value of sensitivity - specificity + 1.
}
\item{auc}{
The area under the curve.
}
}

%\references{
%
%}

\author{
 Michail Tsagris

R implementation and documentation: Giorgos Athineou <athineou@csd.uoc.gr> and Michail Tsagris <mtsagris@yahoo.gr>
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{ \link{SES}, \link{MMPC}, \link{testIndLogistic} 
}
}

\examples{
# simulate a dataset with continuous data
g = rbinom(150, 1, 0.6)
f = rnorm(150)
auc(g, f, roc = TRUE)
}

\keyword{ Area under the curve }
\keyword{ receiver operating curve }
