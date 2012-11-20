\name{INDTRACK3}
\alias{INDTRACK3}
\docType{data}
\title{
  INDTRACK3: FTSE 100 Index and Constituents
}
\description{
  Weekly price data of of the FTSE 100 index and 89 constituents.
}

\usage{data(INDTRACK3)}

\format{
  A data frame with 291 weekly observations of of the index and 89
  members of the FTSE 100 index. The sample starts in March 1991 and
  ends in September 1997.     
}

\details{
  The data set was used in the first two references below. Stocks with 
  missing values during the sample period have been discarded. The data
  was downloaded from DATASTREAM and has been anonymized. The first
  columne refers to the index data itself. See the attached
  license file that is part of this package: \file{BeasleyLicence}.
}

\source{
  \url{http://people.brunel.ac.uk/~mastjjb/jeb/info.html}\cr
  \url{http://people.brunel.ac.uk/~mastjjb/jeb/orlib/legal.html}
}

\references{
  Canakgoz, N.A. and J.E. Beasley (2008), Mixed-integer programming
  approaches for index tracking and enhanced indexation, \emph{European
  Journal of Operational Research}, Vol. 196, 384--399.\cr
  Beasley, J.E. and N. Meade and T.-J. Chang (2003), An evolutionary
  heuristic for the index tracking problem, \emph{European Journal of
    Operational Research}, Vol. 148, 621--643.\cr
  Beasley, J. E. (1990), OR-Library: Distributing Test Problems by
  Electronic Mail, \emph{Journal of the Operational Research Society},
  Vol. 41, No. 11, 1069--1072.
}

\examples{
data(INDTRACK3)
}
\keyword{datasets}
