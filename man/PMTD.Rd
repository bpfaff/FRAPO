\name{PMTD}
\alias{PMTD}

\title{
  Minimum Tail Dependent Portfolio
}

\description{
  This function computes the solution of a minimum tail dependent
  portfolio (long-only).
}

\usage{
PMTD(Returns, method = c("EmpTC", "EVT"), k = NULL, percentage = TRUE,
     optctrl = ctrl(),...)
}

\arguments{
  \item{Returns}{A rectangular array of return data.}
  \item{method}{Character, the type of non-parametric estimation.}
  \item{k}{Integer, the threshold value for the order statistic. If left
  \code{NULL}, then \code{k = sqrt(nrow(x))} is used.}
  \item{percentage}{Logical, whether the weights shall be returned as
    decimals or percentages (default).}
  \item{optctrl}{Object of class \code{Rcpp_CTRL}.}
  \item{\dots}{Arguments are passed down to \code{rank}.}
}

\details{
  Akin to the optimisation of a global minimum-variance portfolio, the
  minimum tail dependennt portfolio is determined by replacing the
  variance-covariance matrix with the matrix of the lower tail
  dependence coefficients as returned by \code{tdc}.
}

\value{
  An object of formal class \code{"PortSol"}.
}

\author{
  Bernhard Pfaff
}

\note{
  The optimisation is conducted by calling \code{cccp()}. 
}

\seealso{
  \code{\link{tdc}}, \code{"\linkS4class{PortSol}"}
}

\examples{
data(StockIndex)
Rets <- returnseries(StockIndex, method = "discrete", trim = TRUE,
                     percentage = TRUE)
PMTD(Rets)
}

\keyword{optimize}
\concept{minimum tail dependence}
