\name{PortCdd-class}
\Rdversion{1.1}
\docType{class}
\alias{PortCdd-class}

\title{Class \code{"PortCdd"}}

\description{
  This class is intended to hold the results from a portfolio
  optimisation with a constraint on its average draw down.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("PortCdd",
  ...)}. This class extends the \code{"PortSol"} class.

}
\section{Slots}{
  \describe{
    \item{\code{CDaR}:}{Numeric, the conditional draw down at risk.}
    \item{\code{thresh}:}{Numeric, threshold value for draw downs at the
      \eqn{\alpha} level.}
    \item{\code{DrawDown}:}{timeSeries, the hsitoric portfolios draw downs.}
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
No methods defined with class "PortCdd" in the signature.
}

\author{
  Bernhard Pfaff
}

\seealso{
  \code{"\linkS4class{PortSol}"}, \code{"\linkS4class{PortMdd}"}, \code{"\linkS4class{PortAdd}"} 
}

\examples{
showClass("PortCdd")
}

\keyword{classes}
