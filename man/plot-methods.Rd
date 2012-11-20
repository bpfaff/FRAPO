\name{plot-methods}
\docType{methods}
\alias{plot-methods}
\alias{plot,PortDD-method}

\title{Methods for Function \code{plot} in Package \pkg{graphics}}

\description{
  Additional arguments to the \code{plot}-method pertinent to the
  defined S4-classes in this package are detailed below.
}

\usage{
\S4method{plot}{PortDD}(x, main = NULL, xlab = NULL, ylab = NULL,
col = c("black", "red"), grid = TRUE, invert = TRUE, \ldots)  
}

\arguments{
  \item{x}{\code{PortDD}: an object that belongs to this virtual class.}
  \item{main}{\code{character}: The title of the plot.}
  \item{xlab}{\code{character}: The description of the x-axis.}
  \item{ylab}{\code{character}: The description of the y-axis.}
  \item{col}{\code{character}: Two-element vector of the names of the
  colors for the portfolio's draw downs and the optimal level.}
  \item{grid}{\code{Logical}: Whether to superimpose a grid on the plot.}
  \item{invert}{\code{Logical}: Whether the draw downs shall be plotted
  as negative numbers; the default is \code{TRUE}.} 
  \item{\ldots}{Ellipsis argument is passed to the generic \code{plot}
  function.} 
}

\author{Bernhard Pfaff}

\keyword{methods}
