\name{MXM-package}
\alias{MXM-package}
\docType{package}
\title{
This is an R package that currently implements feature selection methods for identifying minimal, statistically-equivalent and equally-predictive feature subsets. In addition, two algorithms for constructing the skeleton of a Bayesian network are included.
}

\description{
'MXM' stands for Mens eX Machina, meaning 'Mind from the Machine' in Latin. The package provides source code for the SES algorithm and for some appropriate statistical conditional independence tests. (Fisher and Spearman correlation, G-square test are some examples. Currently the response variable can be univariate or multivariate Euclidean, proportions within 0 and 1, compositional data without zeros and  ones, binary, nominal or ordinal multinomial, count data (handling also overdispersed and with more zeros than expected), longitudinal, clustered data, survival and case-control. Robust versions are also available in some cases and a K-fold cross validation is offered. Bayesian network related algorithms and ridge reression are also included. Read the package's help pages for more details.

MMPC and SES can handle even thousands of variables and for some tests, even many sample sizes of tens of thousands. The user is best advised to check his variables in the beginning. For some regressions, logistic and Poisson for example, we have used C++ codes for speed reasons. Thus no check is done for a variable with zero variance for instance. Something like "colVars" (in "Rfast") could be used in the first place to remove variables with zero variance. 

For more information the reader is addressed to Lagani V., Athineou G., Farcomeni A., Tsagris M. and Tsamardinos I. (2017). Feature Selection with the R Package MXM: Discovering Statistically Equivalent Feature Subsets. Journal of Statistical Software, 80(7), doi:10.18637/jss.v080.i07.
}

\details{
\tabular{ll}{
Package: \tab MXM\cr
Type: \tab Package\cr
Version: \tab 1.3.7\cr
Date: \tab 2018-07-26\cr
License: \tab GPL-2\cr
}
}

\author{
Giorgos Athineou <athineou@csd.uoc.gr>, Vincenzo Lagani <vlagani@csd.uoc.gr>, Michail Tsagris <mtsagris@csd.uoc.gr>, Giorgos Borboudakis <borbudak@csd.uoc.gr>, Ioannis Tsamardinos <tsamard@csd.uoc.gr>, Anna Roumpelaki <anna.roumpelaki@gmail.com>
}

\section{Maintainer}{
Michail Tsagris <mtsagris@csd.uoc.gr>
}

\references{
Lagani V., Athineou G., Farcomeni A., Tsagris M. and Tsamardinos I. (2017). Feature Selection with the R Package MXM: Discovering Statistically Equivalent 
Feature Subsets. Journal of Statistical Software, 80(7), doi:10.18637/jss.v080.i07.

Vincenzo Lagani, George Kortas and Ioannis Tsamardinos (2013), Biomarker signature identification in "omics" with multiclass outcome. 
Computational and Structural Biotechnology Journal, 6(7):1-7.

I. Tsamardinos, V. Lagani and D. Pappas (2012) Discovering multiple, equivalent biomarker signatures. In proceedings of the 7th conference of the 
Hellenic Society for Computational Biology & Bioinformatics - HSCBB12.

Tsagris M., Lagani V., & Tsamardinos I. (2018). Feature selection for high-dimensional glmm data. BMC bioinformatics, 19(1), 17.

Tsagris M., Borboudakis G., Lagani V. and Tsamardinos I. (2018). Constraint-based Causal Discovery with Mixed Data. 
International Journal of Data Science and Analytics. 

Tsamardinos, Brown and Aliferis (2006). The max-min hill-climbing Bayesian network structure learning algorithm. Machine learning, 65(1), 31-78.

Tsamardinos I., Aliferis C. F. and Statnikov, A. (2003). Time and sample efficient discovery of Markov blankets and direct causal relations. 
In Proceedings of the 9th ACM SIGKDD international conference on Knowledge discovery and data mining p. 673-678. 

Pati Y. C., Rezaiifar R. & Krishnaprasad P. S. (1993). Orthogonal matching pursuit: Recursive function approximation with applications to wavelet 
decomposition. In Signals, Systems and Computers. 1993 Conference Record of The Twenty-Seventh Asilomar Conference on. IEEE.

Davis G. (1994). Adaptive Nonlinear Approximations. PhD thesis. 
http://www.geoffdavis.net/papers/dissertation.pdf

Mallat S. G. & Zhang Z. (1993). Matching pursuits with time-frequency dictionaries. IEEE Transactions on signal processing, 41(12), 3397-3415.
https://www.di.ens.fr/~mallat/papiers/MallatPursuit93.pdf

Gharavi-Alkhansari M., & Huang T. S. (1998, May). A fast orthogonal matching pursuit algorithm. In Acoustics, Speech and Signal Processing, 1998. 
Proceedings of the 1998 IEEE International Conference on (Vol. 3, pp. 1389-1392). IEEE.

Chen S., Billings S. A., & Luo W. (1989). Orthogonal least squares methods and their application to non-linear system identification. 
International Journal of control, 50(5), 1873-1896.
http://eprints.whiterose.ac.uk/78100/1/acse%20report%20343...pdf

Borboudakis G. and Tsamardinos I. (2017). Forward-Backward Selection with Early Dropping. 
https://arxiv.org/pdf/1705.10770.pdf

Tsagis M. (2018). Guide on performing feature selection with the R package MXM.
http://mensxmachina.org/site/wp-content/uploads/2018/04/Guide-on-performing-feature-selection-with-the-R-package-MXM.pdf
}

\note{
Acknowledgments:
The research leading to these results has received funding from the European Research Council under the 
European Union's Seventh Framework Programme (FP/2007-2013) / ERC Grant Agreement n. 617393.

Michail Tsagris would like to express his acknowledgments to Marios Dimitriadis and Manos Papadakis, 
undergraduate students in the department of computer science, university of Crete, for their programming tips and advice. 
Dr David Gomez Cabrero Lopez is also acknowledged. Dr Uwe Ligges and Prof Kurt Hornik from the CRAN team are greatly 
acknowledged for their assistance. Prof Achim Zeileis is grealty acknowledged for this help with the quasi Poisson and quasi 
binomial regression models. Christina Chatzipantsiou and Kleio Verrou are acknowledged for their suggestions. 

\bold{Disclaimer:} Professor Tsamardinos is the creator of this package. Dr Lagani supervised Mr Athineou built it
and Dr Tsagris is the maintainer.
}

\keyword{ SES }
\keyword{ Multiple Feature Signatures }
\keyword{ Feature Selection }

\seealso{
\code{\link{SES}, \link{MMPC}, \link{censIndCR},\link{testIndFisher}, \link{testIndLogistic}, \link{gSquare}, \link{testIndRQ}}
}

