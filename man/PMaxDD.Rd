\name{PMaxDD}
\alias{PMaxDD}

\title{
  Portfolio optimisation with maximum draw down constraint
}

\description{
  This function returns the result of a long-only portfolio optimization
  whereby the portfolio's (historic) draw down is constrained to an
  upper limit.
}

\usage{
PMaxDD(PriceData, MaxDD = 0.1, softBudget = FALSE, ...)
}

\arguments{
  \item{PriceData}{A rectangular array of price data.}
  \item{MaxDD}{Numeric, the upper bound of the maximum draw down.}
  \item{softBudget}{Logical, whether the budget constraint shall be
  implemented as a soft constraint, \emph{i.e.} the sum of the weights
  can be less than one. The default is to use an equality constraint.}
  \item{\dots}{Arguments are passed down to \code{Rglpk_solve_LP}}
}

\details{
  This function implements a long-only portfolio optimisation with a
  maximum draw down constraint (see references below). The problem can
  be stated in the form of a linear program and GLPK is used as solver.
}

\value{
  An object of formal class \code{"PortMdd"}.
}

\references{
  Chekhlov, A. and Uryasev, S. and Zabarankin, M., Portfolio
  Optimization with Drawdown Constraints, Department of Industrial and
  Systems Engineering, University of Florida, \emph{Research Report
    2000-5}, 2000, Gainesville, FL.
  Chekhlov, A. and Uryasev, S. and Zabarankin, M., Drawdown Measure in
  Portfolio Optimization, \emph{International Journal of Theoretical and
  Applied Finance}, 2005, 8(1), 13--58.
}

\author{
  Bernhard Pfaff
}

\note{
  A warning is issued in case the solver had exit status not equal to
  zero. 
}

\seealso{
    \code{"\linkS4class{PortSol}"}, \code{"\linkS4class{PortMdd}"},
    \code{"\linkS4class{PortDD}"}, \code{\link{PCDaR}},
    \code{\link{PAveDD}}, \code{\link{PMinCDaR}}  
}

\examples{
\dontrun{
data(StockIndex)
popt <- PMaxDD(PriceData = StockIndex, MaxDD = 0.1, softBudget = TRUE)
}
}

\keyword{optimize}
\concept{portfolio optimisation}
\concept{draw down}
\concept{maximum draw down}