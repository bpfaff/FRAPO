\name{EuroStoxx50}
\alias{EuroStoxx50}
\docType{data}
\title{
  EURO STOXX 50
}
\description{
  Weekly price data of 48 EURO STOXX 50 constituents.
}

\usage{data(EuroStoxx50)}

\format{
  A data frame with 265 weekly observations of 48 members of the EURO
  STOXX 50 index. The sample starts at 2003-03-03 and ends in
  2008-03-24.   
}

\details{
  The data set was used in the reference below. The authors adjusted the
  price data for dividends and have removed stocks if two or more
  consecutive missing values were found. In the remaining cases the NA
  entries have been replaced by interpolated values.
}

\source{
  \url{http://w3.uniroma1.it/Tardella/datasets.html}\cr
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
data(EuroStoxx50)
}
\keyword{datasets}
