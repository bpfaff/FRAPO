\name{mrc}

\alias{mrc}

\title{
  Marginal Contribution to Risk
}

\description{
  This function returns the marginal contributions to portfolio risk,
  whereby the latter is defined in terms of the portfolio standard
  deviation. 
}

\usage{
mrc(weights, Sigma, percentage = TRUE)
}

\arguments{
  \item{weights}{Vector: portfolio weights.}
  \item{Sigma}{Matrix: Variance-covariance matrix of portfolio assets.}
  \item{percentage}{\code{Logical}, whether the marginal risk
  contributions shall be returned as percentages that sum to 100
  (default) or as decimal numbers.}
}

\details{
  The marginal contributions to risk are computed for a given dispersion
  matrix and weight vector.
}

\value{
  \code{numeric}, the marginal risk contributions of the portfolio's
  asset. 
}

\author{
  Bernhard Pfaff
}

\keyword{algebra}
