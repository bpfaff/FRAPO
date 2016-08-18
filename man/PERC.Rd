\name{PERC}
\alias{PERC}

\title{
  Equal risk contributed portfolios
}

\description{
  This function solves for equal risk contributed portfolio weights.
}

\usage{
PERC(Sigma, par = NULL, percentage = TRUE, optctrl = ctrl(), ...)
}

\arguments{
  \item{Sigma}{Matrix, the variance-covariance matrix of asset returns}
  \item{par}{Vector, the initial values of the weights.}
  \item{percentage}{Logical, whether the weights shall be returned as
    decimals or percentages (default).}
  \item{optctrl}{Object of class \code{Rcpp_CTRL}.}
  \item{\dots}{Ellipsis argument is passed down to \code{nlminb()}.}
}

\details{
  The objective function is the standard deviation of the marginal risk
  contributions, which is minimal, \emph{i.e.} zero, if all
  contributions are equal. The weights are rescaled to sum to unity. 
}

\value{
  An object of formal class \code{"PortSol"}.
}

\references{
  Maillard, S. and Roncalli, T. and Teiletche, J.: The Properties
  of Equally Weighted Risk Contribution Portfolios, \emph{Journal of
  Portfolio Management}, Vol. 36, No. 4, Summer 2010, 60--70. 
}

\author{
  Bernhard Pfaff
}

\note{
  The optimisation is conducted by calling \code{nlminb()}. Hereby, the
  arguments \code{lower = 0} and \code{upper = 1} have been specified.  
}

\seealso{
  \code{"\linkS4class{PortSol}"}
}

\examples{
data(MultiAsset)
Rets <- returnseries(MultiAsset, method = "discrete", trim = TRUE,
                     percentage = TRUE)
V <- cov(Rets)
ERC <- PERC(V)
ERC
w <- Weights(ERC)
w * V %*% w
}

\keyword{optimize}
\concept{ERC}
\concept{equal risk contribution}
