\name{trdsma}
\Rdversion{1.1}
\alias{trdsma}
\alias{trdsma-methods}
\alias{trdsma,data.frame-method}
\alias{trdsma,matrix-method}
\alias{trdsma,mts-method}
\alias{trdsma,numeric-method}
\alias{trdsma,timeSeries-method}
\alias{trdsma,ts-method}
\title{
  Simple Moving Average
}
\description{
  Calculation of a right ended simple moving average with equal weights
  determined by \code{n.periods}.
}
\usage{
trdsma(y, n.periods, trim = FALSE)
}

\arguments{
  \item{y}{Objects of classes: numeric, matrix, data.frame, ts, mts,
    and timeSeries are supported.}
  \item{n.periods}{Integer, the number of periods to be included in the
  calculation of the simple moving average.}
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

\value{
  An object of the same class as \code{y}, containing the computed
  simple moving averages.
}

\author{
  Bernhard Pfaff
}
\seealso{
  \code{\link{filter}}, \code{\link{trdbilson}},
  \code{\link{trdbinary}}, \code{\link{trdhp}},
  \code{\link{trdwma}}, \code{\link{capser}},
  \code{\link{trdes}}
}

\examples{
data(StockIndex)
y <- StockIndex[, "SP500"]
sma <- trdsma(y, n.periods = 24)
head(sma, 30)
}
\keyword{ts}
\keyword{smooth}
\keyword{methods}
\concept{Simple Moving Averages}
\concept{Trend}
\concept{Moving Average}