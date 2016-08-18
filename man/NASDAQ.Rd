\name{NASDAQ}
\alias{NASDAQ}
\docType{data}
\title{
  NASDAQ
}
\description{
  Weekly price data of 2,196 NASDAQ constituents.
}

\usage{data(NASDAQ)}

\format{
  A data frame with 265 weekly observations of 2196 members of the
  NASDAQ index. The sample starts at 2003-03-03 and ends in  2008-03-24.
}

\details{
  The data set was used in the reference below. The authors adjusted the
  price data for dividends and have removed stocks if two or more
  consecutive missing values were found. In the remaining cases the NA
  entries have been replaced by interpolated values.
}

\source{
  \url{http://host.uniroma3.it/docenti/cesarone/DataSets.htm}\cr
  \url{http://finance.yahoo.com/}
}

\references{
  Cesarone, F. and Scozzari, A. and Tardella, F.: Portfolio selection
  problems in practice: a comparison between linear and quadratic
  optimization models, Working Paper, Universita degli Studi Roma Tre,
  Universita Telematica delle Scienze Umane and Universita di Roma, July
  2010.\cr
  \url{http://arxiv.org/ftp/arxiv/papers/1105/1105.3594.pdf}
}

\examples{
data(NASDAQ)
}
\keyword{datasets}
