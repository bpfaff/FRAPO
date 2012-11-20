\name{PortAdd-class}
\Rdversion{1.1}
\docType{class}
\alias{PortAdd-class}

\title{Class \code{"PortAdd"}}

\description{
  This class is intended to hold the results from a portfolio
  optimisation with a constraint on its average draw down.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("PortAdd",
  ...)}. This class extends the \code{"PortSol"} class.

}
\section{Slots}{
  \describe{
    \item{\code{AveDD}:}{Numeric, the average draw down.}
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
No methods defined with class "PortAdd" in the signature.
}

\author{
  Bernhard Pfaff
}

\seealso{
  \code{"\linkS4class{PortSol}"}, \code{"\linkS4class{PortMdd}"},
  \code{"\linkS4class{PortCdd}"}  
}

\examples{
showClass("PortAdd")
}

\keyword{classes}
