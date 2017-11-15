\name{Symmetric conditional independence test with clustered data}
\alias{glmm.ci.mm}
\title{
Symmetric conditional independence test with clustered data
}

\description{
Symmetric conditional independence test with clustered data.
}

\usage{
glmm.ci.mm(ind1, ind2, cs = NULL, dat, id) 
}

\arguments{
\item{ind1}{
The index of the one variable to be considered. 
}
\item{ind2}{
The index of the other variable to be considered. 
}
\item{cs}{
The index or indices of the conditioning set of variable(s). If you have no variables set this equal to 0.
}
\item{dat}{
A numerical matrix with data. 
}
\item{id}{
This is a numerical vector denoting the groups or the subjects.  
}
}

\details{
Two linear random intercept models are fitted, one for each variable and the p-value of the hypothesis test that the other variable is significant is calculated. These two p-values are combined in a meta-analytic way. 
}

\value{
A vector including the test statistic, it's associated p-value and the relevant degrees of freedom. 
}

\references{
Tsagris M., Borboudakis G., Lagani V. and Tsamardinos I. (2017). Constraint-based Causal Discovery with Mixed Data. 
The 2017 ACM SIGKDD Workshop on Causal Discovery, 14/8/2017, Halifax, Nova Scotia, Canada. 
http://nugget.unisa.edu.au/CD2017/papersonly/constraint-based-causal-r1.pdf

Eugene Demidenko (2013). Mixed Models: Theory and Applications with R, 2nd Edition. New Jersey: Wiley \& Sons.
}

\author{
Michail Tsagris
  
R implementation and documentation: Michail Tsagris <mtsagris@csd.uoc.gr>
}

%\note{
  %%  ~~further notes~~
    %}
\seealso{
\code{ \link{pc.skel}, \link{condi}, \link{testIndGLMM} }
}

\examples{
## we generate two independent vectors of clustered data
s1 <- matrix(1.5, 4, 4)
diag(s1) <- 2.5
s2 <- matrix(1.5, 4, 4)
diag(s2) <- 2
x1 <- MASS::mvrnorm(10, rnorm(4), s1)  
x1 <- as.vector( t(x1) )
x2 <- MASS::mvrnorm(10, rnorm(4), s2)  
x2 <- as.vector( t(x2) )
id <- rep(1:10, each = 4)
glmm.ci.mm(1, 2, dat = cbind(x1,x2), id = id)
}
