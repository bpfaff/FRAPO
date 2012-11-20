\name{PortDD-class}
\Rdversion{1.1}
\docType{class}
\alias{PortDD-class}
\alias{DrawDowns}
\alias{plot}
\alias{DrawDowns,PortDD-method}
\alias{plot,PortDD,missing-method}

\title{Class \code{"PortDD"}}

\description{
  Class union of \code{"PortAdd"}, \code{"PortCdd"} and \code{"PortMdd"} 
}

\section{Objects from the Class}{
  A virtual Class: No objects may be created from it.}

\section{Methods}{
  \describe{
    \item{DrawDowns}{\code{signature(object = "PortDD")}: Returns the
          portfolio draw downs.} 
    \item{plot}{\code{signature(object = "PortDD")}: Time series plot of
          draw downs.}
  }
}

\author{
  Bernhard Pfaff
}

\note{
  This virtual class is intended for specifying methods that are common
  to all type of draw down portfolios.
}

\seealso{
  \code{"\linkS4class{PortAdd}"}, \code{"\linkS4class{PortMdd}"},
  \code{"\linkS4class{PortCdd}"}, \code{\link{PMinCDaR}},
  \code{\link{PCDaR}}, \code{\link{PAveDD}}, \code{\link{PMaxDD}}   
}

\examples{
showClass("PortDD")
}

\keyword{classes}
