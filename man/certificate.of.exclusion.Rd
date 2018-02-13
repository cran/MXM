\name{Certificate of exclusion from the selected variables set using SES or MMPC}
\alias{certificate.of.exclusion}
\title{
Certificate of exclusion from the selected variables set using SES or MMPC
}

\description{
Information on why one ore more variables were not selected.
}

\usage{
certificate.of.exclusion(xIndex, sesObject = NULL, mmpcObject = NULL) 
}

\arguments{
\item{xIndex}{
A numerical vector with the  indices of the predictor variables.  
}
\item{sesObject}{
If you ran SES, wald.ses or perm.ses, give the whole SES object here, otherwise leave it NULL.
}
\item{mmpcObject}{
If you ran MMPC, wald.mmpc or prm.mmpc, give the whole SES object here, otherwise leave it NULL.
}
}

%\details{
%}

\value{
A list with the conditioning variables (if any), the test statistic and the logarithm of the p-value. In case a variable has been selected a message appears.
}

\author{
Michail Tsagris
  
R implementation and documentation: Michail Tsagris <mtsagris@csd.uoc.gr>
}

\seealso{
\code{\link{MMPC} }
}

\examples{
set.seed(123)
require("hash", quietly = TRUE) 
#simulate a dataset with continuous data
dataset <- matrix(runif(1000 * 200, 1, 100), ncol = 200)
#define a simulated class variable 
target <- 3 * dataset[, 10] + 2 * dataset[, 100] + 3 * dataset[, 20] + rnorm(1000, 0, 5)
# define some simulated equivalences
dataset[, 15] <- dataset[, 10] + rnorm(1000, 0, 2)
dataset[, 150] <- dataset[, 100] + rnorm(1000, 0, 2) 
dataset[, 130] <- dataset[, 100] + rnorm(1000, 0, 2)
# run the SES algorithm
mod1 <- SES(target, dataset, max_k = 5, threshold = 0.05, test = "testIndFisher", 
hash = TRUE, hashObject = NULL);
mod2 <- MMPC(target, dataset, max_k = 5, threshold = 0.05, test = "testIndFisher", 
hash = TRUE, hashObject = NULL);
certificate.of.exclusion(c(10,15,30,45,20), mod1)
certificate.of.exclusion(c(10,15,30,45,20), NULL, mod2)
}

