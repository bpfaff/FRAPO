\name{SocpPhase1}
\alias{SocpPhase1}

\title{SOCP: Initialising objective variable x in primal form}

\description{
  This function determines values for \code{x}, whence they have not
  been specified by the user. Here, a feasibility problem is solved
  first and its solution is then used as an initial point for the
  original problem.}

\usage{
SocpPhase1(f, A, b, N, control)
}

\arguments{
  \item{f}{\code{vector}: the parameters of the objective function in
  its primal form.}
  \item{A}{\code{matrix}; the parameter matrix of the cone constraints.}
  \item{b}{\code{vector}: the parameter vector of the cone constraints.}
  \item{N}{\code{vector}: the count of rows pertinent to each cone constraint.}
  \item{control}{\code{list}: the list of control parameters for SOCP.}
}

\details{
  The finding of an initial point \code{x} is described in the user's
  guide, sectionb 2.8.}

\value{A vector with the initial point for \code{x}.}

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
  \code{\link{Socp}}, \code{\link{SocpPhase2}}, \code{\link{SocpControl}}
}

\keyword{optimize}
