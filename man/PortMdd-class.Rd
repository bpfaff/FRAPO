\name{PortMdd-class}
\Rdversion{1.1}
\docType{class}
\alias{PortMdd-class}

\title{Class \code{"PortMdd"}}

\description{
  This class is intended to hold the results from a portfolio
  optimisation with a constraint on its maximum draw down.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("PortMdd",
  ...)}. This class extends the \code{"PortSol"} class.

}
\section{Slots}{
  \describe{
    \item{\code{MaxDD}:}{Numeric, the maximum draw down.}
    \item{\code{DrawDown}:}{timeSeries, the hsitoric portfolio's draw downs.}
    \item{\code{weights}:}{Numeric, vector of optimal weights.}
    \item{\code{opt}:}{List, the result of the call to GLPK.}
    \item{\code{type}:}{Character, the type of the optimized portfolio.}
    \item{\code{call}:}{The call to the function that created the object.}
  }
}

\section{Extends}{
Class \code{"\linkS4class{PortSol}"}, directly.
}

\section{Methods}{
No methods defined with class "PortMdd" in the signature.
}

\author{
  Bernhard Pfaff
}

\seealso{
  \code{"\linkS4class{PortSol}"}, \code{"\linkS4class{PortAdd}"}, \code{"\linkS4class{PortCdd}"} 
}

\examples{
showClass("PortMdd")
}

\keyword{classes}
