\name{PGMV}
\alias{PGMV}

\title{
  Global Minimum Variance Portfolio
}

\description{
  This function returns the solution of the global minimum variance portfolio
  (long-only). 
}

\usage{
PGMV(Returns, percentage = TRUE, ...)
}

\arguments{
  \item{Returns}{A rectangular array of return data.}
  \item{percentage}{Logical, whether the weights shall be returned as
  decimals or percentages (default).}
  \item{\dots}{Arguments are passed down to \code{cov}.}
}

\value{
  An object of formal class \code{"PortSol"}.
}

\author{
  Bernhard Pfaff
}

\note{
   The optimisation is conducted by calling \code{solve.QP()}. 
}

\seealso{
  \code{"\linkS4class{PortSol}"}
}

\examples{
data(MultiAsset)
Rets <- returnseries(MultiAsset, method = "discrete", trim = TRUE)
PGMV(Rets)
}

\keyword{optimize}
\concept{GMV}
\concept{global minimum variance}