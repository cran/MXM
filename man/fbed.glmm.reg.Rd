\name{Forward Backward Early Dropping selection regression for GLMM}
\alias{fbed.glmm.reg}
\title{
Forward Backward Early Dropping selection regression for GLMM
}

\description{
Forward Backward Early Dropping selection regression for GLMM.
}

\usage{
fbed.glmm.reg(y, x, id, alpha = 0.05, wei = NULL, K = 0, method = "LR", gam = NULL, 
backward = TRUE, type = "gaussian")
}

\arguments{
\item{y}{
The class variable. This can be a numerical vector with continuous data, binary or discrete valued data. It can also be a factor variable with two levels only.
}
\item{x}{
The dataset; provide a numerical a matrix (columns = variables, rows = samples).
}
\item{id}{
This is a numerical vector of the same size as target denoting the groups or the subjects. 
}
\item{alpha}{
Threshold (suitable values in [0,1]) for asmmmbsing p-values significance. Default value is 0.05.
}
\item{wei}{
A vector of weights to be used for weighted regression. The default value is NULL. It is not suggested when robust is set to TRUE.
}
\item{K}{
How many times should the process be repated? The default value is 0. 
}
\item{method}{
Do you want the likelihood ratio test to be performed ("LR" is the default value) or perform the selection using the "eBIC" criterion (BIC is a special case). 
}
\item{gam}{
In case the method is chosen to be "eBIC" one can also specify the \eqn{gamma} parameter. The default value is "NULL", so that the value is 
automatically calculated.
}
\item{backward}{
After the Forward Early Dropping phase, the algorithm proceeds witha the usual Backward Selection phase. The default value is set to TRUE. It is advised to perform this step as maybe some variables are false positives, they were wrongly selected. 

The backward phase using likelihood ratio test and eBIc are two different functions and can be called directly by the user. SO, if you want for example to perform a backard regression with a different threshold value, just use these two functions separately. 
}
\item{type}{
This is for the type of regression to be used, "gaussian", "logistic" or "poisson".
}
%\item{ncores}{
%How many cores to use. This plays an important role if you have tens of thousands of variables or really large sample sizes and tens of thousands of variables and a regression based test which requires numerical optimisation. In other cases it will not make a difference in the overall time (in fact it can be slower). The parallel computation is used in the first step of the algorithm, where univariate associations are examined, those take place in parallel. We have seen a reduction in time of 50\% with 4 cores in comparison to 1 core. Note also, that the amount of reduction is not linear in the number of cores.
    %}
}

\details{
The algorithm is a variation of the usual forward selection. At every step, the most significant variable enters the selected variables set. In addition, only the significant variables stay and are further examined. The non signifcant ones are dropped. This goes until no variable can enter the set. The user has the option to redo this step 1 or more times (the argument K). In the end, a backward selection is performed to remove falsely selected variables. 

Bear in mind that for the "gaussian" case, the forward phase takes place using the F test (Wald statistic calibrated against an F distribution). The backward phase though takes place using the log-likelihood ratio test.
}

\value{
A list including:
\item{res}{
A matrix with the selected variables, their test statistic and the associated p-value.
}
\item{info}{
A matrix with the number of variables and the number of tests performed (or models fitted) at each round (value of K). This refers to the
forward phase only.
}
\item{runtime}{
The runtime required.
}
\item{back.rem}{
The variables removed in the backward phase.
}
\item{back.n.tests}{
The number of models fitted in the backward phase.
}
}

\references{
Borboudakis G. and Tsamardinos I. (2017). Forward-Backward Selection with Early Dropping. 

https://arxiv.org/pdf/1705.10770.pdf
  
Eugene Demidenko (2013). Mixed Models: Theory and Applications with R, 2nd Edition. New Jersey: Wiley \& Sons.
}


\author{
Michail Tsagris
  
R implementation and documentation: Michail Tsagris <mtsagris@csd.uoc.gr>
}

\seealso{
\code{ \link{fbed.glmm.reg}, \link{glmm.bsreg}, \link{MMPC.temporal} }
}

\examples{
## require(lme4)
## data(sleepstudy)
## attach(sleepstudy)
## x <- matrix(rnorm(180 * 20),ncol = 20) ## unrelated preidctor variables
## m1 <- fbed.glmm.reg(Reaction, x, Subject) 
## m2 <- MMPC.temporal(target = Reaction, group = Subject, dataset = x)
}

