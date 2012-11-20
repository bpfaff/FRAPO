\name{ESCBFX}
\alias{ESCBFX}
\docType{data}
\title{
  ESCB FX Reference Rates
}
\description{
  Daily spot rates of major currencies against the EUR.
}

\usage{data(ESCBFX)}

\format{
  A data frame with 3,427 daily observations of the spot currency rates
  rates AUD, CAD, CHF, GBP, HKD, JPY and USD against EUR. The sample
  starts in 1999-01-04 and ends in 2012-04-04.
}

\details{
  The data has been retrieved from the Statistical Data Warehouse (SDW)
  Internet-Site of the ECB. In case of missing data entries due to
  holidays, the last observed data point has been carried forward.
}

\source{
  \url{http://sdw.ecb.europa.eu}
}

\examples{
data(ESCBFX)
}
\keyword{datasets}
