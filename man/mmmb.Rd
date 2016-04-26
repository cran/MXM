\name{mmmb}
\alias{mmmb}
\title{
mmmb: Feature selection algorithm for identifying multiple minimal, statistically-equivalent and equally-predictive feature signatures.

}

\description{
The mmmb algorithm follows a forward-backward filter approach for feature selection in order to provide a minimal, highly-predictive, feature subset of a high dimensional dataset. See also Details. 
}

\usage{
mmmb(target , dataset , max_k = 3 , threshold = 0.05 , test = "testIndFisher", 
user_test = NULL, robust = FALSE, ncores = 1, hold = FALSE)
}

\arguments{
  \item{target}{
The class variable. Provide either a string, an integer, a numeric value, a vector, a factor, an ordered factor or a Surv object. See also Details.
}
  \item{dataset}{
The dataset; provide either a data frame or a matrix (columns = variables, rows = samples). In either case, only two cases are avaialble, either all data are continuous, or categorical. 

}
  \item{max_k}{
The maximum conditioning set to use in the conditional indepedence test (see Details). Integer, default value is 3.
}
  \item{threshold}{
Threshold (suitable values in [0,1]) for asmmmbsing p-values significance. Default value is 0.05.
}
  \item{test}{
The conditional independence test to use. Default value is "testIndFisher". See also \code{link{CondIndTests}}.
}
  \item{user_test}{
A user-defined conditional independence test (provide a closure type object). Default value is NULL. If this is defined, the "test" argument is ignored.
}
\item{robust}{
A boolean variable which indicates whether (TRUE) or not (FALSE) to use a robust version of the statistical test if it is available. 
It takes more time than a non robust version but it is suggested in case of outliers. Default value is FALSE.
}
\item{ncores}{
How many cores to use. This plays an important role if you have tens of thousands of variables or really large sample sizes and tens of thousands of variables and a regression based test which requires numerical optimisation. In other cammmb it will not make a difference in the overall time (in fact it can be slower). The parallel computation is used in the first step of the algorithm, where univariate associations are examined, those take place in parallel. We have seen a reduction in time of 50\% with 4 cores in comparison to 1 core. Note also, that the amount of reduction is not linear in the number of cores.
}
\item{hold}{
After backward (or symmetry correction) phase is implemented. This will remove any possibly falsely included variables in the parents and children set of the target variable and it will slow down the algorithm. If hold is TRUE, even if some variables are identified as falsely included, they will remain. If there are highly collinear (or statistically equivalent) variables, this phase tends to remove correctly identified variables, simply because it will identify a variable wich is highly collinear with the target variable. In this case, the hold should be TRUE. Can you know this in advnace? Well, maybe you can run the SES algorithm to get an idea, or be suspicious about it. 
}
}

\details{
The idea is to run the MMPC algorithm at first and identify the parents and children (PCt) of the target variable. As a second step, the MMPC algorithm is run on the discovered variables to return PCi. The parents of the children of the target are the spouses of the target. Every variable in PCi is checked to see if it is a spouse of the target. If yes, it is included in the Markov Blanket of the target, otherwise it is thrown. If the data are continous, the Fisher correlation test is used or the Spearman correlation (more robust). If the data are categorical, the \eqn{G^2} test is used. 
}

\value{
The output of the algorithm is S3 object including:
\item{mb}{
The Markov Blanket of the target variable. The parents and children of the target variable, along with the spouses of the target, which are the parents of the children of the target variable.
}
\item{runtime}{
The run time of the algorithm. A numeric vector. The first element is the user time, the second element is the system time and the third element is the elapsed time.
}
}

\references{
Tsamardinos I., Aliferis C. F. and Statnikov, A. (2003). Time and sample efficient discovery of Markov blankets and direct causal relations. In Proceedings of the 9th ACM SIGKDD international conference on Knowledge discovery and data mining (pp. 673-678). 
}


\author{
Ioannis Tsamardinos, Vincenzo Lagani

R implementation and documentation: Giorgos Athineou <athineou@csd.uoc.gr> Vincenzo Lagani <vlagani@csd.uoc.gr>  
}

\seealso{
\code{\link{CondIndTests}, \link{MMPC}, \link{SES}}
}

\examples{
set.seed(123)
#require(gRbase) #for faster computations in the internal functions
require(hash)

#simulate a dataset with continuous data
dataset <- matrix( runif(1000 * 50, 1, 100), ncol = 50 )

#define a simulated class variable 
target <- 3 * dataset[, 10] + 2 * dataset[, 50] + 3 * dataset[, 20] + rnorm(1000, 0, 5)

aa <- mmmb(target , dataset , max_k = 3 , threshold = 0.05, test= "testIndFisher", robust = FALSE, 
ncores = 1, hold = FALSE)
ab <- SES(target, dataset, test="testIndFisher")
}

\keyword{ Markov Blanket }
\keyword{ Variable Selection }
