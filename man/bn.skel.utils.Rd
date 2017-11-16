\name{Utilities for the skeleton of a (Bayesian) Network}
\alias{bn.skel.utils}
\title{
Utilities for the skeleton of a (Bayesian) Network
}

\description{
Utilities for the skeleton of a (Bayesian) Network
}

\usage{
bn.skel.utils(mod, G = NULL, roc = TRUE, alpha = 0.01) 
}

\arguments{
\item{mod}{
An object as retured by pc.skel, glmm.pc.skel or mmhc.skel. 
}
\item{G}{
The true adjacency matrix with 1 indicating an edge and zero its absence. Symmetric or not is not important. If this is not available, 
leave it NULL.
}
\item{roc}{
Do you want a graph with the ROC curve be returned? Default value is TRUE.
}
\item{alpha}{
The significance level ( suitable values in (0, 1) ) for assessing the p-values. Default value is 0.01.
}
}

\details{
Given the true adjaceny matrix one can evaluate the estimated adjacency matrix, skeleton, of the Pc or the MMHC algorithm. 
Area under the curve, false discovery rate and sorting of the edges based on their p-values are the main utilites of this function.
}

\value{
A list including:
\item{fdr}{
The false discovery rate as estimated using the Benjamini-Hochberg correction. 
}
\item{area}{
This is a list with the elements of the \code{\link{auc}} function. The area under the curve, the sensitivy and specificity for a range of values, the Youden index, etc.
}
\item{sig.pvalues}{
A matrix with the row and column of each significant p-value sorted in asending order. As we move down the matrix, the p-values increase and hence the strength of the associations decreases.
}
}
\references{
Spirtes P.,  Glymour C. and Scheines R. (2001). Causation, Prediction, and Search. The MIT Press, Cambridge, MA, USA, 3nd edition.

Tsamardinos I. and Brown L.E. Bounding the False Discovery Rate in Local Bayesian Network Learning. AAAI, 2008.

Benjamini Y. and Hochberg Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. 
Journal of the Royal Statistical Society Series B, 57(1), 289-300.
}

\author{
 Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@csd.uoc.gr>
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{ \link{pc.skel}, \link{mmhc.skel}, \link{corfs.network}, \link{local.mmhc.skel}  }
}

\examples{
## simulate a dataset with continuous data
y <- rdag2(300, p = 25, nei = 3)
x <- y$x
a <- pc.skel( x, method = "comb.fast", alpha = 0.05 ) 
bn.skel.utils(a, y$G, roc = FALSE, alpha = 0.01) 
}
