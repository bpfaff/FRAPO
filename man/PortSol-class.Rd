\name{PortSol-class}
\Rdversion{1.1}
\docType{class}
\alias{PortSol-class}
\alias{Weights}
\alias{Solution}
\alias{show,PortSol-method}
\alias{Solution,PortSol-method}
\alias{Weights,PortSol-method}
\alias{update,PortSol-method}

\title{Class \code{"PortSol"}}

\description{
  This class is intended to hold the results for the weights
  of an optimal portfolio. Currently, this class is used for
  minimum-variance and equal-risk-contributed portfolios. It can further
  be used to store the results of optimal factor weights according to
  one of the aforementioned portfolio types.
}

\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("PortSol", ...)}.
  }
       
\section{Slots}{
  \describe{
    \item{\code{weights}:}{Numeric, vector of optimal weights.}
    \item{\code{opt}:}{List, the result of the call to the optimizing function.}
    \item{\code{type}:}{Character, the type of the optimized portfolio.}
    \item{\code{call}:}{The call to the function that created the object.}
  }
}

\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "PortSol")}: Returns the
      portfolio type as text with the optimal weights from the object.}
    \item{Solution}{\code{signature(object = "PortSol")}: Returns the
      list object of the optimizer, \emph{i.e.} the slot \code{opt} from
      the object.} 
    \item{Weights}{\code{signature(object = "PortSol")}: Returns the list
      object of the optimizer, \emph{i.e.} the slot \code{weights} from
      the object.} 
    \item{update}{\code{signature(object = "PortSol")}: updates object
      by calling the issuing function with altered arguments.}  
  }
}

\author{
  Bernhard Pfaff
}

\examples{
showClass("PortSol")
}

\keyword{classes}
