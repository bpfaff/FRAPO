\name{trdes}
\Rdversion{1.1}
\alias{trdes}
\alias{trdes-methods}
\alias{trdes,data.frame-method}
\alias{trdes,matrix-method}
\alias{trdes,mts-method}
\alias{trdes,numeric-method}
\alias{trdes,timeSeries-method}
\alias{trdes,ts-method}
\alias{trdes,xts-method}
\alias{trdes,zoo-method}
\title{
  Exponentially Smoothed Trend
}
\description{
  Calculation of the exponentially smoothed trend as a technical trading indicator.
}
\usage{
trdes(y, lambda, init = NULL)
}

\arguments{
  \item{y}{Objects of classes: numeric, matrix, data.frame, ts, mts,
    timeSeries, zoo and xts are supported.}
  \item{lambda}{Numeric, the smoothing parameter for \eqn{\lambda} in
    the equation below. The value for the parameter must be in the
  interval \eqn{0 < \lambda < 1}.}
  \item{init}{The initial value in the recursive calculation of the
  filter. Specifies the initial values of the time series just prior to
  the start value, in reverse time order. The default, \emph{i.e.}
  \code{NULL}, is a set of zeros.}
}

\details{
  The exponetially smoothed trend is calculated according to the formula:
  \deqn{
    z_t = \lambda y_t + (1 - \lambda) * z_{t-1}    
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
  exponetially smoothed values.
}

\author{
  Bernhard Pfaff
}
\seealso{
  \code{\link{filter}}, \code{\link{trdbilson}},
  \code{\link{trdbinary}}, \code{\link{trdhp}},
  \code{\link{trdsma}}, \code{\link{trdwma}},
  \code{\link{capser}}
}

\examples{
data(StockIndex)
y <- StockIndex[, "SP500"]
yret <- diff(log(y))
es <- trdes(yret, lambda = 0.95)
head(es)
}
\keyword{ts}
\keyword{smooth}
\keyword{methods}
\concept{Exponential Smoothing}
\concept{Trend}
\concept{ES Trend}