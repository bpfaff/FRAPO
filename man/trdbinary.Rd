\name{trdbinary}
\Rdversion{1.1}
\alias{trdbinary}
\alias{trdbinary-methods}
\alias{trdbinary,data.frame-method}
\alias{trdbinary,matrix-method}
\alias{trdbinary,mts-method}
\alias{trdbinary,numeric-method}
\alias{trdbinary,timeSeries-method}
\alias{trdbinary,ts-method}
\title{
  Binary Trend
}
\description{
  Calculation of the Binary Trend as a technical trading indicator.
}
\usage{
trdbinary(y)
}

\arguments{
  \item{y}{Objects of classes: numeric, matrix, data.frame, ts, mts,
    and timeSeries are supported.}
}

\details{
  The Binary trend is calculated according to the formula:
  \deqn{
    z = sign(y) \times \min(|4 / \pi \arctan(y)|, 1)   
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
   \item{y = "numeric"}{Calculation of the binary trend.}
   \item{y = "timeSeries"}{The calculation is
      applied per column of the timeSeries object and an object of the
      same class is returned.}
   \item{y = "ts"}{Calculation of the binary
     trend. The attributes are preserved and an object of the same class is returned.}
   \item{y = "xts"}{Calculation of the binary
     trend. The attributes are preserved and an object of the same class is returned.}
   \item{y = "zoo"}{Calculation of the binary
     trend. The attributes are preserved and an object of the same class is returned.}
 }
}

\value{
  An object of the same class as \code{y}, containing the computed
  Binary trend values.
}

\author{
  Bernhard Pfaff
}
\seealso{
  \code{\link{trdbilson}}, \code{\link{trdes}},
  \code{\link{trdhp}}, \code{\link{trdsma}},
  \code{\link{trdwma}}, \code{\link{capser}} 
}

\examples{
data(StockIndex)
y <- StockIndex[, "SP500"]
yret <- diff(log(y))
binary <- trdbinary(yret)
head(binary)
}
\keyword{ts}
\keyword{smooth}
\keyword{methods}
\concept{Binary}
\concept{Trend}
\concept{Binary Trend}