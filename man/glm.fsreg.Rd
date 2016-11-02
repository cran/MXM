\name{Forward selection with generalised linear regression models}
\alias{glm.fsreg}
\title{
Variable selection in generalised linear regression models with forward selection
}

\description{
Variable selection in generalised linear regression models with forward selection
}

\usage{
glm.fsreg(target, dataset, ini = NULL, threshold = 0.05, wei = NULL, tol = 2, 
heavy = FALSE, robust = FALSE, ncores = 1) 
}

\arguments{
\item{target}{
The class variable. Provide either a string, an integer, a numeric value, a vector, a factor, an ordered factor or a Surv object. See also Details.
}
\item{dataset}{
The dataset; provide either a data frame or a matrix (columns = variables, rows = samples). In either case, only two cases are avaialble, either all data are continuous, or categorical. 
}
\item{ini}{
If you have a set of variables already start with this one. Currently this can only be a matrix with continuous variables. In such cases, the dataset must also contain continuous variables only. Otherwise leave it NULL. 
}
\item{threshold}{
Threshold (suitable values in [0,1]) for assessing the p-values significance. Default value is 0.05.
}
\item{wei}{
A vector of weights to be used for weighted regression. The default value is NULL. It is not suggested when robust is set to TRUE.
}
\item{tol}{
The difference bewtween two successive values of the stopping rule. By default this is is set to 2. If for example, the BIC difference between two succesive models is less than 2, the process stops and the last variable, even though significant does not enter the model.
}
\item{heavy}{
A boolean variable specifying whether heavy computations are required or not. If for exmaple the dataset contains tens of thousands of rows, it is advised to used memory efficient GLMs and hence set this to TRUE.
}
\item{robust}{
A boolean variable which indicates whether (TRUE) or not (FALSE) to use a robust version of the statistical test if it is available. It takes more time than a non robust version but it is suggested in case of outliers. 
Default value is FALSE and this is currently supported only for the linear regression.
}
\item{ncores}{
How many cores to use. This plays an important role if you have tens of thousands of variables or really large sample sizes and tens of thousands of variables and a regression based test which requires numerical optimisation. In other cammmb it will not make a difference in the overall time (in fact it can be slower). The parallel computation is used in the first step of the algorithm, where univariate associations are examined, those take place in parallel. We have seen a reduction in time of 50\% with 4 cores in comparison to 1 core. Note also, that the amount of reduction is not linear in the number of cores.
}
%\item{maxit}{
%The maximum number of iterations required by the MM estimators in the robust binary or poisson regrression. This argument is %taken into account only if the argument robust is set to TRUE.
%}
}

%\details{
%
%}

\value{
The output of the algorithm is S3 object including:
\item{mat}{
A matrix with the variables and their latest test statistics and p-values.
}
\item{info}{
A matrix with the selected variables, their p-values and test statistics. Each row corresponds to a model which contains the variables up to that line. The BIC in the last column is the BIC of that model.
}
\item{models}{
The regression models, one at each step.
}
\item{final}{
The final regression model.
}
\item{ci_test}{
The conditional independence test used. 
}
\item{runtime}{
The run time of the algorithm. A numeric vector. The first element is the user time, the second element is the system time and the third element is the elapsed time.
}
}

%\references{
%Tsamardinos I., Aliferis C. F. and Statnikov, A. (2003). Time and sample efficient discovery of Markov blankets and direct causal relations. In Proceedings of the 9th ACM SIGKDD international conference on Knowledge discovery and data mining (pp. 673-678). 
%}


\author{
Michail Tsagris

R implementation and documentation: Giorgos Athineou <athineou@csd.uoc.gr> Michail Tsagris <mtsagris@yahoo.gr>
}

\seealso{
\code{\link{fs.reg}, \link{lm.fsreg}, \link{bic.fsreg}, \link{bic.glm.fsreg}, \link{CondIndTests}, \link{MMPC}, \link{SES}}
}

\examples{
set.seed(123)

#simulate a dataset with continuous data
dataset <- matrix( runif(1000 * 50, 1, 100), ncol = 50 )

#define a simulated class variable 
target <- rpois(1000, 10)

a <- glm.fsreg(target, dataset, threshold = 0.05, tol = 2, robust = FALSE, ncores = 1 ) 
b <- MMPC(target, dataset, max_k = 3, threshold = 0.05, test = "testIndPois")
}

\keyword{ Markov Blanket }
\keyword{ Variable Selection }
