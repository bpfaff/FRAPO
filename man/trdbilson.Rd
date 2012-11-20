\name{trdbilson}
\Rdversion{1.1}
\alias{trdbilson}
\alias{trdbilson-methods}
\alias{trdbilson,data.frame-method}
\alias{trdbilson,matrix-method}
\alias{trdbilson,mts-method}
\alias{trdbilson,numeric-method}
\alias{trdbilson,timeSeries-method}
\alias{trdbilson,ts-method}
\alias{trdbilson,xts-method}
\alias{trdbilson,zoo-method}
\title{
  Bilson Trend
}
\description{
  Calculation of the Bilson Trend as a technical trading indicator.
}
\usage{
trdbilson(y, exponent)
}

\arguments{
  \item{y}{Objects of classes: numeric, matrix, data.frame, ts, mts,
    timeSeries, zoo and xts are supported.}
  \item{exponent}{Numeric, the value for \eqn{\alpha} in the equation below.}
}

\details{
  The Bilson trend is calculated according to the formula:
  \deqn{
    z = sign(y) \times |y|^{(1 - |y|^\alpha)}    
    }
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
   \item{y = "numeric"}{Calculation of the bilson trend.}
   \item{y = "timeSeries"}{The calculation is
      applied per column of the timeSeries object and an object of the
      same class is returned.}
   \item{y = "ts"}{Calculation of the bilson
     trend. The attributes are preserved and an object of the same class is returned.}
   \item{y = "xts"}{Calculation of the bilson
     trend. The attributes are preserved and an object of the same class is returned.}
   \item{y = "zoo"}{Calculation of the bilson
     trend. The attributes are preserved and an object of the same class is returned.}
 }
}

\value{
  An object of the same class as \code{y}, containing the computed
  Bilson trend values.
}

\author{
  Bernhard Pfaff
}
\seealso{
  \code{\link{trdbinary}}, \code{\link{trdes}},
  \code{\link{trdhp}}, \code{\link{trdsma}},
  \code{\link{trdwma}}, \code{\link{capser}} 
}

\examples{
data(StockIndex)
y <- StockIndex[, "SP500"]
yret <- diff(log(y))
bilson <- trdbilson(yret, exponent = 2)
head(bilson)
}
\keyword{ts}
\keyword{smooth}
\keyword{methods}
\concept{Bilson}
\concept{Trend}
\concept{Bilson Trend}