\name{FTSE100}
\alias{FTSE100}
\docType{data}
\title{
  FTSE 100
}
\description{
  Weekly price data of 79 FTSE 100 constituents.
}

\usage{data(FTSE100)}

\format{
  A data frame with 265 weekly observations of 79 members of the FTSE
  100 index. The sample starts at 2003-03-03 and ends in  2008-03-24.
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
data(FTSE100)
}
\keyword{datasets}
