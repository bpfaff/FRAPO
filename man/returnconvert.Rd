\name{returnconvert}
\Rdversion{1.1}
\alias{returnconvert}
\alias{returnconvert-methods}
\alias{returnconvert,data.frame-method}
\alias{returnconvert,matrix-method}
\alias{returnconvert,mts-method}
\alias{returnconvert,numeric-method}
\alias{returnconvert,timeSeries-method}
\alias{returnconvert,ts-method}
\title{
  Convert Returns from continuous to discrete and vice versa
}
\description{
  Either continuous returns or discrete returns can be converted into
  the other type. 
}
\usage{
returnconvert(y, convdir = c("cont2disc", "disc2cont"), percentage = TRUE)
}

\arguments{
  \item{y}{Objects of classes: numeric, matrix, data.frame, ts, mts,
    and timeSeries are supported.}
  \item{convdir}{Character, the type of return conversion.}
  \item{percentage}{Logical, if \code{TRUE} (the default) the returns,
  \code{y}, are expressed as percentages.}
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
   \item{y = "numeric"}{Calculation of the returns.}
   \item{y = "timeSeries"}{The calculation is
      applied per column of the timeSeries object and an object of the
      same class is returned.}
   \item{y = "ts"}{Calculation of the returns. The attributes are
  preserved and an object of the same class is returned.} 
   \item{y = "xts"}{Calculation of the returns. The attributes are
  preserved and an object of the same class is returned.} 
   \item{y = "zoo"}{Calculation of the returns. The attributes are
  preserved and an object of the same class is returned.} 
 }
}

\value{
  An object of the same class as \code{y}, containing the converted returns.
}

\author{
  Bernhard Pfaff
}

\examples{
data(StockIndex)
yc <- diff(log(StockIndex[, "SP500"])) * 100
yd <- returnseries(StockIndex[, "SP500"], method = "discrete",
                   percentage = TRUE, trim = TRUE)
yconv <- returnconvert(yd, convdir = "disc2cont",
                       percentage = TRUE)
all.equal(yc, yconv)
}
\keyword{ts}
\keyword{methods}
\concept{returns}
\concept{discrete returns}
\concept{continuous returns}
