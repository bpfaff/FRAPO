\name{MultiAsset}
\alias{MultiAsset}
\docType{data}
\title{
  Multi Asset Index Data
}
\description{
  Month-end price data of stock and bond indices and gold.
}

\usage{data(MultiAsset)}

\format{
  A data frame with 85 month-end observations of stock and bond indices
  and gold, represented by the ETF SPDR Gold. The sample starts at
  2004-11-30 and ends in 2011-11-30.       
}

\details{
  The data set has been obtained from Yahoo Finance and hereby the
  unadjusted closing prices have been retrieved. If a month-end value
  was not reported, the value of the previous day has been used. The
  Yahoo mnemonics with the respective item description are listed below:

  \describe{
    \item{GSPC}{United States: S \& P 500 Index (Equity)}
    \item{RUA}{United States: Russell 3000 Index (Equity)}
    \item{GDAXI}{Germany: DAX (XETRA) Index (Equity)}
    \item{FTSE}{United Kingdom: FTSE 100 Index (Equity)}
    \item{N225}{Japan: Nikkei 225 Index (Equity)}
    \item{EEM}{iShares: MSCI Emerging Markets Index (Equity)}
    \item{DJCBTI}{United States: Dow Jones CBOT Treasury Index (Bonds)}
    \item{GREXP}{Germany: REX-Performance Index (Bonds)}
    \item{BG05.L}{United Kingdom: Gilt All Index (Bonds)}
    \item{GLD}{United States: SPDR Gold Shares (Commodities)}
    }  
}

\source{
  \url{http://finance.yahoo.com/}
}

\examples{
data(MultiAsset)
}
\keyword{datasets}
