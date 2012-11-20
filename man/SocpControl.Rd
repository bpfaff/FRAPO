\name{SocpControl}
\alias{SocpControl}

\title{Control Variables for Socp}

\description{
  This function returns a \code{list} object of control parameters that
  are passed down to the \code{C}-function \code{SOCP}. It's default
  values are used in \code{Socp}.
}

\usage{
SocpControl(abs.tol = 1e-18, rel.tol = 1e-16, target = 0,
            max.iter = 500, Nu = 10, out.mode = 0, BigM.K = 2,
            BigM.iter = 5)
}

\arguments{
  \item{abs.tol}{Absolute tolerance.}
  \item{rel.tol}{Relative tolerance.}
  \item{target}{Target value < 0, only used if rel.tol < 0.}
  \item{max.iter}{The maximum number of iterations, socp is aborted if
  more are required for convergence.} 
  \item{Nu}{The parameter that controls the rate of convergence, Nu > 1,
  recommended range 5 to 50.} 
  \item{out.mode}{Specifies what should be output: 0 - nothing, 1 -
  duality gap for initial point and after each iteration, 2 - duality
  gap and deviation from centrality, for initial point and after each
  iteration.}
  \item{BigM.K}{Iterataion parameter. The default values is \code{BigM.K
  = 2}.}
  \item{BigM.iter}{Iterataion parameter. The default values is
  \code{BigM.iter = 5}.}
}

\details{
  For details about these control parameters, the reader is referred to
  the reference below, in particular sections 2.7, 2.8 and 4.3 to 4.5. A
  pdf-version of the user's guide is shipped in the packages \code{doc}
  subdirectory. 
}

\value{A \code{list} object with the control parameters.}

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
  \code{\link{Socp}}
}

\keyword{optimize}
