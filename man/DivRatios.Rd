\name{DivMeasures}
\alias{dr}
\alias{cr}
\alias{rhow}

\title{
  Diversification Measures 
}

\description{
  These functions compute the diversification ratio, the volatility
  weighted average correlation and concentration ratio of a portfolio.
}

\usage{
dr(weights, Sigma)
cr(weights, Sigma)
rhow(weights, Sigma)
}

\arguments{
  \item{weights}{Vector: portfolio weights.}
  \item{Sigma}{Matrix: Variance-covariance matrix of portfolio assets.}
}

\details{
  The diversification ratio of a portfolio is defined as:

  \deqn{DR(\omega) = \frac{\sum_{i = 1}^N \omega_i \sigma_i}{\sqrt{\omega'
	\Sigma \omega}}
  }

  for a portfolio of \eqn{N} assets and \eqn{\omega_i} signify the
  weight of the i-th asset and \eqn{\sigma_i} its standard deviation and
  \eqn{\Sigma} the variance-covariance matrix of asset returns. The
  diversification ratio is therefore the weighted average of the assets'
  volatilities divided by the portfolio volatility.

  The concentration ration is defined as:

  \deqn{CR = \frac{\sum_{i = 1}^N (\omega_i \sigma_i)^2}{(\sum_{i = 1}^N
  \omega_i \sigma_i)^2} 
  }

  and the volatility-weighted average correlation of the assets as:

  \deqn{\rho(\omega) = \frac{\sum_{i > j}^N (\omega_i \sigma_i \omega_j
  \sigma_j)\rho_{ij}}{\sum_{i > j}^N (\omega_i \sigma_i \omega_j
  \sigma_j)}  
  }

  The following equation between these measures does exist:

  \deqn{DR(\omega) = \frac{1}{\sqrt{\rho(\omega) (1 - CR(\omega)) +
	CR(\omega)}}}
}

\value{
  \code{numeric}, the value of the diversification measure.
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

\seealso{
  \code{\link{PMD}}
}

\examples{
data(MultiAsset)
Rets <- returnseries(MultiAsset, method = "discrete", trim = TRUE)
w <- Weights(PMD(Rets))
V <- cov(Rets)
DR <- dr(w, V)
CR <- cr(w, V)
RhoW <- rhow(w, V)
test <- 1 / sqrt(RhoW * (1 - CR) + CR)
all.equal(DR, test) 
}

\keyword{algebra}
