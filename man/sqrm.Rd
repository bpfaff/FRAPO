\name{sqrm}
\alias{sqrm}

\title{
  Square root of a quadratic matrix
}

\description{
  This function returns the square root of a quadratic and
  diagonalisable matrix.
}

\usage{
sqrm(x, ...)
}

\arguments{
  \item{x}{\code{matrix}, must be quadratic.}
  \item{...}{The ellipsis argument is passed down to \code{eigen()}.}
}
\details{
  The computation of the square root of a matrix is based upon its eigen
  values and corresponding eigen vectors. The square matrix \eqn{A} is
  diagonisable if there is a matrix \eqn{V} such that \eqn{D = V^{-1}AV},
  whereby \eqn{D} is a diagonal matrix. This is only achieved if the eigen
  vectors of the \eqn{(n \times n)} matrix \eqn{A} constitute a basis of
  dimension \eqn{n}. The square root of \eqn{A} is then \eqn{A^{1/2} = V
  D^{1/2} V'}.  
}

\value{
  A \code{matrix} object and a scalar in case a \eqn{(1 \times 1)} matrix has been
  provided. 
}

\author{
  Bernhard Pfaff
}

\seealso{
  \code{\link[base]{eigen}}
}

\examples{
data(StockIndex)
S <- cov(StockIndex)
SR <- sqrm(S)
all.equal(crossprod(SR), S)
}

\keyword{array}
\keyword{algebra}
\concept{square root of matrix}
