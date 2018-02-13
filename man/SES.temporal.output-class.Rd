\name{SES.temporal.output-class}
\Rdversion{1.1}
\docType{class}
\alias{SES.temporal.output-class}
\alias{SES.temporal.output}
\alias{plot,SES.temporal.output-method}
\alias{SES.temporal.output-method}
\alias{plot,SES.temporal.output,ANY-method}


\title{Class \code{"SES.temporal.output"}}
\description{
SES.temporal output object class.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SES.temporal.output", ...)}.
}

\section{Slots}{
  \describe{
    \item{\code{selectedVars}:}{Object of class \code{"numeric"}}
    \item{\code{selectedVarsOrder}:}{Object of class \code{"numeric"}}
    \item{\code{queues}:}{Object of class \code{"list"}}
    \item{\code{signatures}:}{Object of class \code{"matrix"}}
    \item{\code{hashObject}:}{Object of class \code{"list"}}
    \item{\code{pvalues}:}{Object of class \code{"numeric"}}
    \item{\code{stats}:}{Object of class \code{"numeric"}}
    \item{\code{univ}:}{Object of class \code{"list"}}
    \item{\code{max_k}:}{Object of class \code{"numeric"}}
    \item{\code{threshold}:}{Object of class \code{"numeric"}}
    \item{\code{runtime}:}{Object of class \code{"proc_time"}}
    \item{\code{test}:}{Object of class \code{"character"}}
    \item{\code{slope}:}{Object of class \code{"logical"}}  }
}

\section{Methods}{
  \describe{
    \item{plot}{\code{plot(x = "SES.temporal.output", mode = "all")}: Generic function for plotting the generated pvalues of the SES.temporal.output object. Argument mode = "all" for plotting all the pvalues or mode="partial" for partial plotting the first 500 pvalues }
	 }
}

\author{
Giorgos Athineou <athineou@csd.uoc.gr>
}

\seealso{
\code{\link{MMPC.temporal}, \link{SES.temporal}}
}

\examples{
showClass("SES.temporal.output")
}

\keyword{SES.temporal output}
