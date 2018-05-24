\name{MMPC.temporal.output-class}
\Rdversion{1.1}
\docType{class}
\alias{MMPC.temporal.output-class}
\alias{MMPC.temporal.output}
\alias{plot,MMPC.temporal.output-method}
\alias{MMPC.temporal.output-method}
\alias{plot,MMPC.temporal.output,ANY-method}

\title{Class \code{"MMPC.temporal.output"}}

\description{
MMPC.temporal output object class.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("MMPC.temporal.output", ...)}.
}

\section{Slots}{
  \describe{
    \item{\code{selectedVars}:}{Object of class \code{"numeric"}}
    \item{\code{selectedVarsOrder}:}{Object of class \code{"numeric"}}
    \item{\code{hashObject}:}{Object of class \code{"list"}}
    \item{\code{pvalues}:}{Object of class \code{"numeric"}}
    \item{\code{stats}:}{Object of class \code{"numeric"}}
    \item{\code{univ}:}{Object of class \code{"list"}}
    \item{\code{max_k}:}{Object of class \code{"numeric"}}
    \item{\code{threshold}:}{Object of class \code{"numeric"}}
    \item{\code{n.tests}:}{Object of class \code{"numeric"}}
    \item{\code{runtime}:}{Object of class \code{"proc_time"}}
    \item{\code{test}:}{Object of class \code{"character"}}
    \item{\code{slope}:}{Object of class \code{"logical"}}  }
}

\section{Methods}{
  \describe{
    \item{plot}{\code{plot(x = "MMPC.temporal.output", mode = "all")}: Generic function for plotting the generated pvalues of the MMPC.temporal.output object. Argument mode = "all" for plotting all the pvalues or mode="partial" for partial plotting the first 500 pvalues }
	 }
}

\author{
Giorgos Athineou <athineou@csd.uoc.gr>
}

\seealso{
\code{\link{MMPC.temporal}, \link{SES.temporal}}
}

\examples{
showClass("MMPC.temporal.output")
}

\keyword{MMPC.temporal output}
