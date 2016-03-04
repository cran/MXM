\name{MXM-package}
\alias{MXM-package}
\docType{package}
\title{
This is an R package that currently implements feature selection methods for identifying minimal, statistically-equivalent and equally-predictive feature subsets. In addition, two algorithms for constructing the skeleton of a Bayesian network are included.
}
\description{
'MXM' stands for Mens eX Machina, meaning 'Mind from the Machine' in Latin. The package provides source code for the SES algorithm and for some appropriate statistical conditional independence tests (Fisher and Spearman corelation, G-square test, linear, logistic, poisson and Cox regression are some of the included tests). An automated k-fold cross validation function for the SES algortighm is included. Two algorithms for obtaining the skeleton of a Bayesian network and ridge reression are also included. Read the package's help pages for more details.
}
\details{
\tabular{ll}{
Package: \tab MXM\cr
Type: \tab Package\cr
Version: \tab 0.7\cr
Date: \tab 2016-03-04\cr
License: \tab GPL-2\cr
}
}
\author{
Giorgos Athineou <athineou@csd.uoc.gr>, Vincenzo Lagani <vlagani@csd.uoc.gr>, Michail Tsagris <mtsagris@yahoo.gr>, Giorgos Borboudakis <borbudak@csd.uoc.gr>, Ioannis Tsamardinos <tsamard@csd.uoc.gr>
}
\section{Maintainers }{
Giorgos Athineou <athineou@csd.uoc.gr>, Vincenzo Lagani <vlagani@csd.uoc.gr>, Michail Tsagris <mtsagris@yahoo.gr>
}
\references{
I. Tsamardinos, V. Lagani and D. Pappas (2012) Discovering multiple, equivalent biomarker signatures. In proceedings of the 7th conference of the Hellenic Society for Computational Biology & Bioinformatics - HSCBB12.
}

\keyword{ SES }
\keyword{ Multiple Feature Signatures }
\keyword{ Feature Selection }
\keyword{ Variable Selection }
\seealso{
\code{\link{SES},\link{censIndLR},\link{testIndFisher}, \link{testIndLogistic}, \link{gSquare}, \link{censIndLR}}
}

