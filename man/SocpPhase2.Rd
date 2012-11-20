\name{SocpPhase2}
\alias{SocpPhase2}

\title{SOCP: Initialising objective variable z in dual form}

\description{
  This function determines values for \code{z}, whence they have not
  been specified by the user.
}

\usage{
SocpPhase2(f, A, b, N, x, control)
}

\arguments{
  \item{f}{\code{vector}: the parameters of the objective function in
  its primal form.}
  \item{A}{\code{matrix}; the parameter matrix of the cone constraints.}
  \item{b}{\code{vector}: the parameter vector of the cone constraints.}
  \item{N}{\code{vector}: the count of rows pertinent to each cone constraint.}
  \item{x}{\code{vector}: initial point of SOCP in its primal form.}
  \item{control}{\code{list}: the list of control parameters for SOCP.}
}

\value{A vector with the initial point for \code{z} (dual form of
  SOCP).}

\references{
  Lobo, M. and Vandenberghe, L. and Boyd, S., \emph{SOCP: Software for
  Second-order Cone Programming, User's Guide}, Beta Version, April
  1997, Stanford University.
}

\author{
  Bernhard Pfaff
}

\note{
  This function has been ported from the \pkg{Rsocp} package contained
  in the Rmetrics-Project on R-Forge.
}

\seealso{
  \code{\link{Socp}}, \code{\link{SocpPhase1}}, \code{\link{SocpControl}}
}

\keyword{optimize}

