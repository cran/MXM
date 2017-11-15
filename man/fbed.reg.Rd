\name{Forward Backward Early Dropping selection regression}
\alias{fbed.reg}
\title{
Forward Backward Early Dropping selection regression
}

\description{
Forward Backward Early Dropping selection regression.
}

\usage{
fbed.reg(y, x, test = NULL, alpha = 0.05, wei = NULL, K = 0, method = "LR", 
gam = NULL, backward = TRUE) 
}

\arguments{
\item{y}{
The class variable. Provide either a string, an integer, a numeric value, a vector, a factor, an ordered factor or a Surv object. 
}
\item{x}{
The dataset; provide either a data frame or a matrix (columns = variables, rows = samples). 
}
\item{test}{
The available tests: "testIndReg", "testIndPois", "testIndNB", "testIndLogistic", "testIndMMreg" 
(in the other algorithms this the robust version of "testIndReg"), "testIndBinom", "censIndCR", 
"censIndWR", "testIndBeta", "testIndZIP", "testIndGamma, "testIndNormLog" and "testIndTobit".
}
\item{alpha}{
Threshold (suitable values in [0,1]) for asmmmbsing p-values significance. Default value is 0.05.
}
\item{wei}{
A vector of weights to be used for weighted regression. The default value is NULL. It is not suggested when robust is set 
to TRUE. If you want to use the "testIndBinom", then supply the successes in the y and the trials here.
}
\item{K}{
How many times should the process be repated? The default value is 0. 
}
\item{method}{
Do you want the likelihood ratio test to be performed ("LR" is the default value) or perform the selection using the "eBIC" criterion (BIC is a special case). 
}
\item{gam}{
In case the method is chosen to be "eBIC" one can also specify the \eqn{gamma} parameter. The default value is "NULL", so that the value is automatically calculated.
}
\item{backward}{
After the Forward Early Dropping phase, the algorithm proceeds witha the usual Backward Selection phase. The default value is set to TRUE. It is advised to perform this step as maybe some variables are false positives, they were wrongly selected. 

The backward phase using likelihood ratio test and eBIc are two different functions and can be called directly by the user. SO, if you want for example to perform a backard regression with a different threshold value, just use these two functions separately. 
}
%\item{ncores}{
%How many cores to use. This plays an important role if you have tens of thousands of variables or really large sample sizes and tens of thousands of variables and a regression based test which requires numerical optimisation. In other cases it will not make a difference in the overall time (in fact it can be slower). The parallel computation is used in the first step of the algorithm, where univariate associations are examined, those take place in parallel. We have seen a reduction in time of 50\% with 4 cores in comparison to 1 core. Note also, that the amount of reduction is not linear in the number of cores.
%}
}

\details{
The algorithm is a variation of the usual forward selection. At every step, the most significant variable enters the selected variables set. In addition, only the significant variables stay and are further examined. The non signifcant ones are dropped. This goes until no variable can enter the set. The user has the option to redo this step 1 or more times (the argument K). In the end, a backward selection is performed to remove falsely selected variables. 
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
}


\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@csd.uoc.gr>
}

\seealso{
\code{ \link{fs.reg}, \link{ebic.bsreg}, \link{bic.fsreg}, \link{MMPC} }
}

\examples{
#simulate a dataset with continuous data
dataset <- matrix( runif(100 * 50, 1, 100), ncol = 50 )
#define a simulated class variable 
target <- rt(100, 10)

a1 <- fbed.reg(target, dataset, alpha = 0.05, test = "testIndReg") 
y <- rpois(100, 10)
a2 <- fbed.reg(y, dataset, alpha = 0.05, test = "testIndPois") 
a3 <- MMPC(y, dataset)
}

