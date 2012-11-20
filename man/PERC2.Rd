\name{PERC2}
\alias{PERC2}

\title{
  Equal risk contributed portfolios (alternative)
}

\description{
  This function solves for equal risk contributed portfolio weights by
  employing the \code{solnp()} function.
}

\usage{
PERC2(Sigma, par = NULL, percentage = TRUE, ...)
}

\arguments{
  \item{Sigma}{Matrix, the variance-covariance matrix of asset returns}
  \item{par}{Vector, the initial values of the weights.}
  \item{percentage}{Logical, whether the weights shall be returned as
  decimals or percentages (default).}
  \item{\dots}{Ellipsis argument is passed down to \code{solnp()}.}
}

\details{
  The objective function is the standard deviation of the marginal risk
  contributions, which is minimal, \emph{i.e.} zero, if all
  contributions are equal. The weights are rescaled to sum to unity. The
  lower and upper bounds are set to \code{LB = rep(0, N)} and \code{UB =
  rep(1, N)}, respectively.  
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
  The optimisation is conducted by calling \code{sonlp()} contained in
  the package \pkg{Rsolnp}. Hereby, the arguments \code{LB = 0} and
  \code{UB = 1} have been specified.   
}

\seealso{
  \code{"\linkS4class{PortSol}"}, \code{\link{PERC}}
}

\examples{
\dontrun{
library(Rsolnp)
data(MultiAsset)
Rets <- returnseries(MultiAsset, method = "discrete", trim = TRUE,
                     percentage = TRUE)
V <- cov(Rets)
## Budget Constraint
Budget <- function(x, Sigma) sum(x)
ERC <- PERC2(V, eqfun = Budget, eqB = 1)
ERC
w <- Weights(ERC) / 100
w * V %*% w
}
}

\keyword{optimize}
\concept{ERC}
\concept{equal risk contribution}
