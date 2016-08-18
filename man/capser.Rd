\name{capser}
\Rdversion{1.1}
\alias{capser}
\alias{capser-methods}
\alias{capser,data.frame-method}
\alias{capser,matrix-method}
\alias{capser,mts-method}
\alias{capser,numeric-method}
\alias{capser,timeSeries-method}
\alias{capser,ts-method}
\title{
  Capping a series to bounds
}
\description{
  The values of a series that are absolute greater than \code{min}
  and/or \code{max} are capped to these specified values.
}
\usage{
capser(y, min, max)
}

\arguments{
  \item{y}{Objects of classes: numeric, matrix, data.frame, ts, mts,
    and timeSeries are supported.}
  \item{min}{Numeric, minimum value for the series.}
  \item{max}{Numeric, maximim value for the series.}
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
   \item{y = "timeSeries"}{The calculation is
      applied per column of the timeSeries object and an object of the
      same class is returned.}
   \item{y = "ts"}{Calculation of the es
     trend. The attributes are preserved and an object of the same class is returned.}
   \item{y = "xts"}{Calculation of the es
     trend. The attributes are preserved and an object of the same class is returned.}
   \item{y = "zoo"}{Calculation of the es
     trend. The attributes are preserved and an object of the same class is returned.}
 }
}

\value{
  An object of the same class as \code{y}, containing the truncated series.
}

\author{
  Bernhard Pfaff
}
\seealso{ 
  \code{\link{trdbilson}}, \code{\link{trdbinary}}, \code{\link{trdes}},
  \code{\link{trdhp}}, \code{\link{trdsma}}, \code{\link{trdwma}}  
}

\examples{
data(StockIndex)
y <- StockIndex[, "SP500"]
cs <- capser(y, min = 100, max = 200)
head(cs)
}
\keyword{ts}
\keyword{smooth}
\keyword{methods}
\concept{capping}
\concept{bounding}
\concept{truncating}