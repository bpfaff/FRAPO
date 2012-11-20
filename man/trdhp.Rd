\name{trdhp}
\Rdversion{1.1}
\alias{trdhp}
\alias{trdhp-methods}
\alias{trdhp,data.frame-method}
\alias{trdhp,matrix-method}
\alias{trdhp,mts-method}
\alias{trdhp,numeric-method}
\alias{trdhp,timeSeries-method}
\alias{trdhp,ts-method}
\alias{trdhp,xts-method}
\alias{trdhp,zoo-method}
\title{
  Hodrick-Prescott Filter
}
\description{
  Calculation of the Hodrick-Prescott filter as a technical trading indicator.
}
\usage{
trdhp(y, lambda)
}

\arguments{
  \item{y}{Objects of classes: numeric, matrix, data.frame, ts, mts,
    timeSeries, zoo and xts are supported.}
  \item{lambda}{Numeric, the value for \eqn{\lambda} in the equation below.}
}

\details{
  The Hodrick-Prescott filter is calculated according to the formula:
  \deqn{
    \min (\tau_t) =  \sum_{t = 1}^T (y_t - \tau_t)^2 + \lambda
    \sum_{t = 2}^{T-1} (\Delta^2 \tau_{t+1})^2
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
  Hodrick-Prescott values.
}

\references{
  Hodrick, R. and E.C. Prescott (1997), Postwar U.S. Business Cycles: An
  Empirical Investigation, \emph{Journal of Money, Credit and Banking} 29(1).
  }

\author{
  Bernhard Pfaff
}

\seealso{
  \code{\link{trdbinary}}, \code{\link{trdes}},
  \code{\link{trdbilson}}, \code{\link{trdsma}},
  \code{\link{trdwma}}, \code{\link{capser}}  
}

\examples{
data(StockIndex)
y <- StockIndex[, "SP500"]
hp <- trdhp(y, lambda = 1600)
head(hp)
}
\keyword{ts}
\keyword{smooth}
\keyword{methods}
\concept{Hodrick-Prescott}
\concept{Trend}
\concept{Filter}
