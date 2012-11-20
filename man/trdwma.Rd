\name{trdwma}
\Rdversion{1.1}
\alias{trdwma}
\alias{trdwma-methods}
\alias{trdwma,data.frame-method}
\alias{trdwma,matrix-method}
\alias{trdwma,mts-method}
\alias{trdwma,numeric-method}
\alias{trdwma,timeSeries-method}
\alias{trdwma,ts-method}
\alias{trdwma,xts-method}
\alias{trdwma,zoo-method}
\title{
 Weighted Moving Average
}
\description{
  Calculation of a right ended weighted moving average with weights
  according to \code{weights}.
}
\usage{
trdwma(y, weights, trim = FALSE)
}

\arguments{
  \item{y}{Objects of classes: numeric, matrix, data.frame, ts, mts,
    timeSeries, zoo and xts are supported.}
  \item{weights}{Numeric, a vector containing the weights.}
  \item{trim}{Logical, if \code{FALSE} (the default) the first
  value is set to \code{NA}, otherwise the object is trimmed by the
  first obeservation.}
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

\details{
  If the sum of the weights is greater than unity, a warning is issued.
}

\value{
  An object of the same class as \code{y}, containing the computed
  weighted moving averages.
}

\author{
  Bernhard Pfaff
}
\seealso{
  \code{\link{filter}}, \code{\link{trdbilson}},
  \code{\link{trdbinary}}, \code{\link{trdhp}},
  \code{\link{trdes}}, \code{\link{trdsma}},
  \code{\link{capser}}
}

\examples{
data(StockIndex)
y <- StockIndex[, "SP500"]
wma <- trdwma(y, weights = c(0.4, 0.3, 0.2, 0.1))
head(wma, 30)
}
\keyword{ts}
\keyword{smooth}
\keyword{methods}
\concept{weighted moving averages}
\concept{Trend}
\concept{moving averages}
