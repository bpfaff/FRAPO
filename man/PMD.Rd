\name{PMD}
\alias{PMD}

\title{
  Most Diversified Portfolio
}

\description{
  This function returns the solution of the most diversified portfolio
  (long-only). 
}

\usage{
PMD(Returns, percentage = TRUE, ...)
}

\arguments{
  \item{Returns}{A rectangular array of return data.}
  \item{percentage}{Logical, whether the weights shall be returned as
  decimals or percentages (default).}
  \item{\dots}{Arguments are passed down to \code{cov()}.}
}

\details{
  The optimisation problem is akin to that of a global minimum-variance
  portfolio, but instead of using the variance-covariance matrix of the
  asset returns, the correlation matrix is utilised as dispersion
  measure. The weights are then recovered by rescaling the optimal
  solution with the assets' standard deviations and normalizing, such
  that the weights sum to one.
}

\value{
  An object of formal class \code{"PortSol"}.
}

\references{
  Choueifaty, Y. and Coignard, Y. (2008): Toward Maximum
  Diversification, \emph{Journal of Portfolio Management}, Vol. 34,
  No. 4, 40--51.

  Choueifaty, Y. and Coignard, Y. and Reynier, J. (2011): Properties of
  the Most Diversified Portfolio, Working Paper, \url{http://papers.ssrn.com}
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
PMD(Rets)
}

\keyword{optimize}
\concept{MDP}
\concept{diversification ratio}
\concept{most diversified portfolio}