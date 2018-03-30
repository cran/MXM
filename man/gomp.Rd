\name{Generalised orthogonal matching pursuit (gOMP)}
\alias{gomp}
\title{
Generalised orthogonal matching pursuit (gOMP)
}

\description{
Generalised orthogonal matching pursuit.
}

\usage{
gomp(target, dataset, tol = qchisq(0.95, 1) + log(dim(dataset)[1]), 
test = "testIndFisher" ) 
}

\arguments{
\item{target}{
The response variable, a numeric vector, a matrix or a Surv object.
}
\item{dataset}{
A matrix with continuous data, where the rows denote the samples and the columns are the variables. 
}
\item{tol}{
The tolerance value to terminate the algorithm. This is the change in the criterion value 
between two successive steps. The default value is the 95\% quantile of the \eqn{\chi^2} distribution
with 1 degree of freedom plus the logarithm of the sample size. For test = "testIndFisher" the BIC is already
calculated, so better make it qchisq(0.95, 1). 
}
\item{test}{
This denotes the parametric model to be used each time. It depends upon the nature of the target variable. 
The possible values are "testIndFisher" (or "testIndReg" for the same purpose), "testIndLogistic", 
"testIndPois", "testIndQPois", "testIndQbinom", "testIndNormLog", "testIndMVreg", "testIndNB", "testIndBeta", 
"testIndGamma", "testIndMMReg", "testIndRQ", "testIndTobit", "censIndCR" and "censIndWR".
}
}

\value{
A list including:
\item{runtime}{
The runtime of the algorithm
}
\item{res}{
A matrix with two columns. The selected variable(s) and the criterion value at every step. 
}
}

\references{
Pati Y. C., Rezaiifar R. & Krishnaprasad P. S. (1993). Orthogonal matching pursuit: Recursive function approximation with applications to wavelet decomposition. In Signals, Systems and Computers. 1993 Conference Record of The Twenty-Seventh Asilomar Conference on. IEEE.

Mazin Abdulrasool Hameed (2012). Comparative analysis of orthogonal matching pursuit and least angle regression. MSc thesis, Michigan State University.
https://www.google.gr/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&ved=0ahUKEwik9P3Yto7XAhUiCZoKHQ8XDr8QFgglMAA&url=https%3A%2F%2Fd.lib.msu.edu%2Fetd%2F1711%2Fdatastream%2FOBJ%2Fdownload%2FComparative_analysis_of_orthogonal_matching_pursuit_and_least_angle_regression.pdf&usg=AOvVaw2fRcSemcbteyWUiAKYi-8B

Razavi S. A. Ollila E., & Koivunen V. (2012). Robust greedy algorithms for compressed sensing. In Signal Processing 
Conference (EUSIPCO), 2012 Proceedings of the 20th European. IEEE.

Lozano A., Swirszcz G., & Abe N. (2011). Group orthogonal matching pursuit for logistic regression. In Proceedings of the Fourteenth International Conference on Artificial Intelligence and Statistics.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{cor.fbed}, \link{cor.fsreg}, \link{correls}, \link{fs.reg}
}
}

\examples{
x <- matrix( rnorm(100 * 100), ncol = 100)
y <- rnorm(100)
b <- MXM::gomp(y, x, test = "testIndFisher")
}

