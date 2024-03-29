\name{Constraint based feature selection algorithms for longitudinal and clustered data}
\alias{SES.glmm}
\alias{MMPC.glmm}
\alias{MMPC.gee}
\alias{SES.gee}
\title{
SES.glmm/SES.gee: Feature selection algorithm for identifying multiple minimal, statistically-equivalent and equally-predictive feature signatures with correlated data
}

\description{
SES.glmm algorithm follows a forward-backward filter approach for feature selection in order to provide minimal, highly-predictive, statistically-equivalent, multiple feature subsets of a high dimensional dataset. See also Details. MMPC.glmm algorithm follows the same approach without generating multiple feature subsets. They are both adapted to longitudinal target variables.
}

\usage{
SES.glmm(target, reps = NULL, group, dataset, max_k = 3, threshold = 0.05, 
test = NULL, ini = NULL, wei = NULL, user_test = NULL, hash = FALSE, 
hashObject = NULL, slopes = FALSE, ncores = 1)

MMPC.glmm(target, reps = NULL, group, dataset, max_k = 3, threshold = 0.05, 
test = NULL, ini = NULL, wei = NULL, user_test = NULL, hash = FALSE, 
hashObject = NULL, slopes = FALSE, ncores = 1)

MMPC.gee(target, reps = NULL, group, dataset, max_k = 3, threshold = 0.05, 
test = NULL, ini = NULL, wei = NULL, user_test = NULL, hash = FALSE, 
hashObject = NULL, correl = "exchangeable", se = "jack", ncores = 1)

SES.gee(target, reps = NULL, group, dataset, max_k = 3, threshold = 0.05, 
test = NULL, ini = NULL, wei = NULL, user_test = NULL, hash = FALSE, 
hashObject = NULL, correl = "exchangeable", se = "jack", ncores = 1)
}

\arguments{
\item{target}{
The class variable. Provide a vector with continuous (normal), binary (binomial) or discrete (Poisson) data.
}
\item{reps}{
A numeric vector containing the time points of the subjects. It's length is equal to the length of the target variable. If you have 
clustered data, leave this NULL. 
}
\item{group}{
A numeric vector containing the subjects or groups. It must be of the same legnth as target. 
}
\item{dataset}{
The dataset; provide either a data frame or a matrix (columns = variables , rows = samples). Currently, only continuous datasets are supported. 
}
\item{max_k}{
The maximum conditioning set to use in the conditional indepedence test (see Details). Integer, default value is 3.
}
\item{threshold}{
Threshold (suitable values in (0, 1)) for assessing p-values significance. Default value is 0.05.
}
\item{test}{
The conditional independence test to use. Default value is NULL. Currently, the only available conditional independence tests are the 
"testIndGLMMLogistic", "testIndGLMMPois", "testIndGLMMGamma", "testIndGLMMNormLog", "testIndGLMMOrdinal", "testIndGLMMReg", "testIndLMM" 
and "testIndGLMMCR" for generalised linear mixed models. For the GEE, the available tests are "testIndGEEReg", "testIndGEEPois", 
"testIndGEELogistic", "testIndGEEGamma" and "testIndGEENormLog". 
}
\item{ini}{
This is a supposed to be a list. After running SES or MMPC with some hyper-parameters you might want to run SES again with different 
hyper-parameters. To avoid calculating the univariate associations (first step of SES and of MPPC) again, you can extract them 
from the first run of SES and plug them here. This can speed up the 
second run (and subequent runs of course) by 50\%. See the details and the argument "univ" in the output values. 
}
\item{wei}{
A vector of weights to be used for weighted regression. The default value is NULL. 
}
\item{user_test}{
A user-defined conditional independence test (provide a closure type object). Default value is NULL. If this is defined, the 
"test" argument is ignored.
}
\item{hash}{
A boolean variable which indicates whether (TRUE) or not (FALSE) to store the statistics calculated during SES execution in a 
hash-type object. Default value is FALSE. If TRUE a hashObject is produced.
}
\item{hashObject}{
A List with the hash objects generated in a previous run of SES.glmm. 
Each time SES runs with "hash=TRUE" it produces a list of hashObjects that can be re-used in order to speed up next runs of SES.

Important: the generated hashObjects should be used only when the same dataset is re-analyzed, possibly with different values of max_k and threshold.
}
\item{slopes}{
Should random slopes for the ime effect be fitted as well? Default value is FALSE. 
}
\item{correl}{
The correlation structure. For the Gaussian, Logistic, Poisson and Gamma regression this can be either "exchangeable" (compound symmetry, 
suitable for clustered data) or "ar1" (AR(1) model, suitable for longitudinal data). For the ordinal logistic regression its only the 
"exchangeable" correlation sturcture. 
}
\item{se}{
The method for estimating standard errors. This is very important and crucial. The available options for Gaussian, Logistic, Poisson and 
Gamma regression are: a) 'san.se': the usual robust estimate. b) 'jack': approximate jackknife variance estimate. c) 'j1s': if 1-step 
jackknife variance estimate and d) 'fij': fully iterated jackknife variance estimate. If you have many clusters (sets of repeated measurements) 
"san.se" is fine as it is asympotically correct, plus jacknife estimates will take longer. If you have a few clusters, then maybe it's better 
to use jacknife estimates. 

The jackknife variance estimator was suggested by Paik (1988), which is quite suitable for cases when the number of subjects is small 
(K < 30), as in many biological studies. The simulation studies conducted by Ziegler et al. (2000) and Yan and Fine (2004) showed that 
the approximate jackknife estimates are in many cases in good agreement with the fully iterated ones. 
}
\item{ncores}{
How many cores to use. This plays an important role if you have tens of thousands of variables or really large sample sizes and tens of 
thousands of variables and a regression based test which requires numerical optimisation. In other cases it will not make a difference 
in the overall time (in fact it can be slower). The parallel computation is used in the first step of the algorithm, where univariate 
associations are examined, those take place in parallel. We have seen a reduction in time of 50\% with 4 cores in comparison to 1 core. 
Note also, that the amount of reduction is definetely not linear in the number of cores.
}
}

\details{
The SES.glmm function implements the Statistically Equivalent Signature (SES) algorithm as presented in "Tsamardinos, Lagani and Pappas, 
HSCBB 2012" adapted to longitudinal data. The citation for this is "Tsagris, Lagani and tsamardinos, (2018)". These functions presented 
here are for the \bold{temporal-lonitudinal scenario}. 

The MMPC function mplements the MMPC algorithm as presented in "Tsamardinos, Brown and Aliferis. The max-min hill-climbing Bayesian network structure learning algorithm" adapted to longitudinal data. 

The output value "univ" along with the output value "hashObject" can speed up the computations of subesequent runs of SES and MMPC. 
The first run with a specific pair of hyper-parameters (threshold and max_k) the univariate associations tests and the conditional 
independence tests (test statistic and logarithm of their corresponding p-values) are stored and returned. In the next run(s) with 
different pair(s) of hyper-parameters you can use this information to save time. With a few thousands of variables you will see the 
difference, which can be up to 50\%. 

The max_k option: the maximum size of the conditioning set to use in the conditioning independence test. Larger values provide more 
accurate results, at the cost of higher computational times. When the sample size is small (e.g., \eqn{<50} observations) the max_k 
parameter should be \eqn{\leq 5}, otherwise the conditional independence test may not be able to provide reliable results.

If the dataset contains missing (NA) values, they will automatically be replaced by the current variable (column) mean value with an 
appropriate warning to the user after the execution.

If the target is a single integer value or a string, it has to corresponds to the column number or to the name of the target feature in 
the dataset. In any other case the target is a variable that is not contained in the dataset.

If the current 'test' argument is defined as NULL or "auto" and the user_test argument is NULL then the algorithm automatically selects 
only available, which is \code{\link{testIndGLMMReg}}.

Conditional independence test functions to be pass through the user_test argument should have the same signature of the included test. 
See "?testIndFisher" for an example.

For all the available conditional independence tests that are currently included on the package, please see "?CondIndTests". 

If two or more p-values are below the machine epsilon (.Machine$double.eps which is equal to 2.220446e-16), all of them are set to 0. 
To make the comparison or the ordering feasible we use the logarithm of the p-value. The max-min heuristic though, requires comparison 
and an ordering of the p-values. Hence, all conditional independence tests calculate the logarithm of the p-value.

If there are missing values in the dataset (predictor variables) columnwise imputation takes place. The median is used for the continuous 
variables and the mode for categorical variables. It is a naive and not so clever method. For this reason the user is encouraged to make 
sure his data contain no missing values. 

If you have percentages, in the (0, 1) interval, they are automatically mapped into \eqn{R} by using the logit transformation and a 
linear mixed model is fitted. If you have binary data, logistic mixed regression is applied and if you have discrete data (counts), 
Poisson mixed regression is applied. 

If you want to use the GEE methodology, make sure you load the library geepack first. 
}

\value{
The output of the algorithm is an object of the class 'SES.glmm.output' for SES.glmm or 'MMPC.glmm.output' for MMPC.glmm including:
\item{selectedVars}{
The selected variables, i.e., the signature of the target variable.
}
\item{selectedVarsOrder}{
The order of the selected variables according to increasing pvalues.
}
\item{queues}{
A list containing a list (queue) of equivalent features for each variable included in selectedVars. An equivalent signature can 
be built by selecting a single feature from each queue. Featured only in SES.
}
\item{signatures}{
A matrix reporting all equivalent signatures (one signature for each row). Featured only in SES.
}
\item{hashObject}{
The hashObject caching the statistic calculted in the current run.
}
\item{pvalues}{
For each feature included in the dataset, this vector reports the strength of its association with the target in the context of all other variables. Particularly, this vector reports the max p-values found when the association of each variable with the target is tested against different conditional sets. Lower values indicate higher association. \bold{Note that these are the logarithm of the p-values}.
}
\item{stats}{
The statistics corresponding to "pvalues" (higher values indicates higher association).
}
\item{univ}{
This is a list with the univariate associations. The test statistics and their corresponding \bold{logged p-values}. This list is very 
important for subsequent runs of SES with different hyper-parameters. After running SES with some hyper-parameters you might want 
to run SES again with different hyper-parameters. To avoid calculating the univariate associations (first step of SES or MMPC) 
again, you can take this list from the first run of SES and plug it in the argument "ini" in the next run(s) of SES or MMPC. 
This can speed up the second run (and subequent runs of course) by 50\%. See the argument "univ" in the output values. 
}
\item{max_k}{
The max_k option used in the current run.
}
\item{threshold}{
The threshold option used in the current run.
}
\item{n.tests}{
If you have set hash = TRUE, then the number of tests performed by SES or MMPC will be returned. If you have not set this to TRUE,
the number of univariate associations will be returned. So be careful with this number.
}
\item{runtime}{
The run time of the algorithm. A numeric vector. The first element is the user time, the second element is the system time and the 
third element is the elapsed time.
}
\item{slope}{
Whether random slopes for the time effects were used or not, TRUE or FALSE.
}

Generic Functions implemented for SESoutput Object:
\item{plot(object=SES.glmm.output, mode="all")}{
Plots the generated pvalues (using barplot) of the current SESoutput object in comparison to the threshold.
Argument mode can be either "all" or "partial" for the first 500 pvalues of the object.
}
}

\references{
Tsagris, M., Lagani, V., & Tsamardinos, I. (2018). Feature selection for high-dimensional glmm data. BMC bioinformatics, 19(1), 17.

I. Tsamardinos, M. Tsagris and V. Lagani (2015). Feature selection for longitudinal data. Proceedings of the 10th conference of 
the Hellenic Society for Computational Biology & Bioinformatics (HSCBB15).

I. Tsamardinos, V. Lagani and D. Pappas (2012). Discovering multiple, equivalent biomarker signatures. In proceedings of the 7th 
conference of the Hellenic Society for Computational Biology & Bioinformatics - HSCBB12.

Tsamardinos, Brown and Aliferis (2006). The max-min hill-climbing Bayesian network structure learning algorithm. Machine learning, 65(1), 31-78.

Eugene Demidenko (2013). Mixed Models: Theory and Applications with R, 2nd Edition. New Jersey: Wiley & Sons.

J. Pinheiro and D. Bates. Mixed-effects models in S and S-PLUS. Springer Science & Business Media, 2006.

Liang K.Y. and Zeger S.L. (1986). Longitudinal data analysis using generalized linear models. Biometrika, 73(1): 13-22.

Prentice R.L. and Zhao L.P. (1991). Estimating equations for parameters in means and covariances of multivariate discrete and 
continuous responses. Biometrics, 47(3): 825-839.

Heagerty P.J. and Zeger S.L. (1996) Marginal regression models for clustered ordinal measurements. Journal of the American Statistical 
Association, 91(435): 1024-1036.

Paik M.C. (1988). Repeated measurement analysis for nonnormal data in small samples. Communications in Statistics-Simulation and 
Computation, 17(4): 1155-1171.

Ziegler A., Kastner C., Brunner D. and Blettner M. (2000). Familial associations of lipid profiles: A generalised estimating 
equations approach. Statistics in medicine, 19(24): 3345-3357

Yan J. and Fine J. (2004). Estimating equations for association structures. Statistics in medicine, 23(6): 859-874.
}

\author{
Ioannis Tsamardinos, Vincenzo Lagani

R implementation and documentation: Giorgos Athineou <athineou@csd.uoc.gr> Vincenzo Lagani <vlagani@csd.uoc.gr>  
}


\seealso{
\code{\link{CondIndTests}, \link{testIndGLMMReg} }
}

\examples{
\dontrun{
require(lme4)
data(sleepstudy)
reaction <- sleepstudy$Reaction
days <- sleepstudy$Days
subject <- sleepstudy$Subject
x <- matrix(rnorm(180 * 200),ncol = 200) ## unrelated predictor variables
m1 <- SES.glmm(target = reaction, reps = days, group = subject, dataset = x)
m2 <- MMPC.glmm(target = reaction, reps = days, group = subject, dataset = x)
}
}

\keyword{ SES }
\keyword{ Mixed models }
\keyword{ Longitudinal data }
\keyword{ Variable Selection }
