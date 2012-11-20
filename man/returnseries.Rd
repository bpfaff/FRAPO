\name{returnseries}
\Rdversion{1.1}
\alias{returnseries}
\alias{returnseries-methods}
\alias{returnseries,data.frame-method}
\alias{returnseries,matrix-method}
\alias{returnseries,mts-method}
\alias{returnseries,numeric-method}
\alias{returnseries,timeSeries-method}
\alias{returnseries,ts-method}
\alias{returnseries,xts-method}
\alias{returnseries,zoo-method}
\title{
  Continuous and discrete returns
}
\description{
  Either continuous returns or discrete returns are computed for an
  object. The returns can be expressed as percenatges and the first \code{NA}
  value can be trimmed.
}
\usage{
returnseries(y, method = c("continuous", "discrete"), percentage = TRUE,
             trim = FALSE, compound = FALSE)
}

\arguments{
  \item{y}{Objects of classes: numeric, matrix, data.frame, ts, mts,
    timeSeries, zoo and xts are supported.}
  \item{method}{Character, the type of return to be computed.}
  \item{percentage}{Logical, if \code{TRUE} (the default) the returns
  are expressed as percenatges.}
  \item{trim}{Logical, if \code{FALSE} (the default) the first value is
  set to \code{NA} such that the length of the return series coincides
  with the length of the series in levels.}
  \item{compound}{Logical, if \code{FALSE} (the default), then simple
  returns are computed and otherwise compounded returns.}
}

\section{Methods}{
  \describe{
    \item{y = "data.frame"}{The calculation is
      applied per column of the data.frame and only if all columns are numeric.}
    \item{y = "matrix"}{The calculation is
      applied per column of the matrix.}
    \item{y = "mts"}{The calculation is
      applied per column of the mts object. The attributes are preserved
      and an object of the same class is returned.}
   \item{y = "numeric"}{Calculation of the es trend.}
   \item{y = "timeSeries"}{The calculation is applied per column of the
  timeSeries object and an object of the same class is returned.} 
   \item{y = "ts"}{Calculation of the returns. The attributes are
  preserved and an object of the same class is returned.}  
   \item{y = "xts"}{Calculation of the returns. The attributes are
  preserved and an object of the same class is returned.} 
   \item{y = "zoo"}{Calculation of the returns. The attributes are
  preserved and an object of the same class is returned.} 
 }
}

\value{
  An object of the same class as \code{y}, containing the truncated series.
}

\author{
  Bernhard Pfaff
}

\examples{
data(StockIndex)
y <- StockIndex[, "SP500"]
ret <- returnseries(y)
head(ret)
}
\keyword{ts}
\keyword{methods}
\concept{returns}
\concept{discrete returns}
\concept{continuous returns}
\concept{compound returns}
