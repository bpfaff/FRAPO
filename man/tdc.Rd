\name{tdc}
\alias{tdc}

\title{
  Tail Dependence Coefficient
}

\description{
  This function returns the pairwise tail dependence coefficients
  between \eqn{N} series. The TDCs are estimated non-parametrically by
  either the empirical tail copula or based on the stable
  tail-dependence function.
}

\usage{
tdc(x, method = c("EmpTC", "EVT"), lower = TRUE, k = NULL, ...)
}

\arguments{
  \item{x}{Matrix, or an object that can be coerced to it.}
  \item{method}{Character, the type of non-parametric estimation.}
  \item{lower}{Logical, if \code{TRUE} (default), lower TDC are computed
  and upper TDC, else.}
  \item{k}{Integer, the threshold value for the order statistic. If left
  \code{NULL}, then \code{k = sqrt(nrow(x))} is used.}
  \item{\dots}{Ellipsis, arguments are passed down to \code{rank}.}
}

\details{
  For a matrix or an object that can be coerced to it with \code{ncol(x)
  >= 2}, the pair wise tail dependencies are estimated
  non-parametrically and returned as a symmetric matrix. The threshold
  value \code{k} is the upper/lower bound for the order statistics to be 
  considered. The diagonal elements are always equal to one, because a
  series has a dependence of one with itself, of course.
}

\value{
  A matrix with the tail dependent coefficients.
}

\references{
  Schmidt, R. and Stadtm\"uller, U., Nonparametric estimation of tail
  dependence, \emph{The Scandinavian Journal of Statistics}, 33,
  307--335. 
}

\author{
  Bernhard Pfaff
}

\seealso{
  \code{\link{PMTD}}
}

\examples{
data(StockIndex)
Rets <- returnseries(StockIndex, method = "discrete", trim = TRUE,
                     percentage = TRUE)
tdc(Rets, method = "EmpTC")
tdc(Rets, method = "EVT")
}

\keyword{algebra}
\concept{TDC}
\concept{tail dependence coefficient}
\concept{lower tail dependence coefficient}
\concept{upper tail dependence coefficient}
\concept{non-parametric tail dependence coefficient}

