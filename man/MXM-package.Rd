\name{MXM-package}
\alias{MXM-package}
\docType{package}
\title{
This is an R package that currently implements feature selection methods for identifying minimal, statistically-equivalent and equally-predictive feature subsets. In addition, two algorithms for constructing the skeleton of a Bayesian network are included.
}

\description{
'MXM' stands for Mens eX Machina, meaning 'Mind from the Machine' in Latin. The package provides source code for the SES algorithm and for some appropriate statistical conditional independence tests. (Fisher and Spearman corelation, G-square test are some examples. Currently the response variable can be univariate or multivariate Euclidean, proportions within 0 and 1, compositional data without zeros and  ones, binary, nominal or ordinal multinomial, count data (handling also overdispersed and with more zeros than expected), longitudinal, clustered data, survival and case-control. Robust versions are also available in some cases and a K-fold cross validation is offered. Bayesian network related algorithms and ridge reression are also included. Read the package's help pages for more details.

MMPC and SES can handle even thousands of variables and for some tests, even many sample sizes of tens of thousands. The user is best advised to check his variables in the beginning. For some regressions, logistic and Poisson for example, we have used C++ codes for speed reasons. Thus no check is done for a variable with zero variance for instance. Something like colVars could be used in the first place to remove variables with zero variance. 

For more information the reader is addressed to Lagani V., Athineou G., Farcomeni A., Tsagris M. and Tsamardinos I. (2017). Feature Selection with the R Package MXM: Discovering Statistically Equivalent Feature Subsets. Journal of Statistical Software, 80(7):1-25, doi:10.18637/jss.v080.i07.
}

\details{
\tabular{ll}{
Package: \tab MXM\cr
Type: \tab Package\cr
Version: \tab 1.1.1\cr
Date: \tab 2017-10-07\cr
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
Tsamardinos, Brown and Aliferis (2006). The max-min hill-climbing Bayesian network structure learning algorithm. Machine learning, 65(1), 31-78.

I. Tsamardinos, V. Lagani and D. Pappas (2012) Discovering multiple, equivalent biomarker signatures. In proceedings of the 7th conference of the Hellenic 
Society for Computational Biology & Bioinformatics - HSCBB12.

Tsamardinos I., Aliferis C. F. and Statnikov, A. (2003). Time and sample efficient discovery of Markov blankets and direct causal relations. In Proceedingsof the 9th ACM SIGKDD international conference on Knowledge discovery and data mining p. 673-678. 
}

\note{
Acknowledgments:
The research leading to these results has received funding from the European Research Council under the 
European Union's Seventh Framework Programme (FP/2007-2013) / ERC Grant Agreement n. 617393.

Michail Tsagris would like to express his acknowledgments to Marios Dimitriadis and Manos Papadakis, 
undergraduate students in the department of computer science, university of Crete, for their programming tips and advice. 
Dr David Gomez Cabrero Lopez is also greatly acknowledged. Dr Uwe Ligges is greatly acknowledged for his assistance. 

\bold{Disclaimer:} Professor Tsamardinos is the creator of this package. Dr Lagani supervised Mr Athineou built it
and Dr Tsagris is the maintainer.
}

\keyword{ SES }
\keyword{ Multiple Feature Signatures }
\keyword{ Feature Selection }

\seealso{
\code{\link{SES}, \link{MMPC}, \link{censIndCR},\link{testIndFisher}, \link{testIndLogistic}, \link{gSquare}, \link{testIndRQ}}
}

