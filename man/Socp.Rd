\name{Socp}
\alias{Socp}

\title{Second-order Cone Programming}

\description{
  The function solves second-order cone problem by
  primal-dual interior point method. It is a wrapper function to the
  \code{C}-routines written by Lobo, Vandenberghe and Boyd (see
  reference below).
}

\usage{
Socp(f, A, b, C, d, N,
     x = NULL, z = NULL, w = NULL, control = list())
}

\arguments{
  \item{f}{Vector defining linear objective, \code{length(f)==length(x)} }
  \item{A}{Matrix with the \eqn{A_i} vertically stacked:  \eqn{A = [
  A_1; A_2; \ldots; A_L]}.}
  \item{b}{Vector with the \eqn{b_i} vertically stacked: \eqn{b = [ b_1;
  b_2; \ldots; b_L]}.} 
  \item{C}{ Matrix with the \eqn{c_i'} vertically stacked: \eqn{C = [ c_1';
        c_2'; \ldots; c_L']}.}
  \item{d}{Vector with the \eqn{d_i} vertically stacked: \eqn{d = [ d_1;
  d_2; \ldots;  d_L]}.} 
  \item{N}{Vector of size \code{L}, defining the size of each constraint.}
  \item{x}{Primal feasible initial point. Must satisfy: \eqn{|| A_i*x +
  b_i || < c_i' * x + d_i} for \eqn{i = 1, \ldots, L}.}
  \item{z}{Dual feasible initial point.}
  \item{w}{Dual feasible initial point.}
  \item{control}{A list of control parameters.}
}

\details{
  The primal formulation of an SOCP is given as:
  \deqn{minimise f' * x}
  subject to
  \deqn{||A_i*x + b_i|| <= c_i' * x + d_i}
  for \eqn{i = 1,\ldots, L}. Here, \eqn{x} is the \eqn{(n \times 1)}
  vector to be optimised. The dual form of an SOCP is expressed as:
  \deqn{maximise \sum_{i = 1}^L -(b' * z_i + d_i * w_i)}
  subject to
  \deqn{\sum_{i = 1}^L (A_i' * z_i + c_i * w_i) = f}
  and
  \deqn{||z_i || = w_i}
  for \eqn{i = 1,\ldots, L}, given strictly feasible primal and dual
  initial points.\cr

  The algorithm stops, if one of the following criteria is met:

  \enumerate{
    \item \code{abs.tol} -- maximum absolute error in objective
    function; guarantees that for any x:  \eqn{abs(f'*x - f'*x\_opt) <=
    abs\_tol}. 
    
    \item \code{rel.tol} -- maximum relative error in objective
    function; guarantees that for any x: \eqn{abs(f'*x -
    f'*x\_opt)/(f'*x\_opt) <= rel\_tol  (if f'*x\_opt > 0)}. Negative
    value has special meaning, see target next.
    
    \item \code{target} -- if \eqn{rel\_tol<0}, stops when
    \eqn{f'*x < target or -b'*z >= target}.
    
    \item \code{max.iter} -- limit on number of algorithm outer iterations.
    Most problems can be solved in less than 50 iterations. Called with
    \code{max_iter = 0} only checks feasibility of \code{x} and \code{z},
    (and returns gap and deviation from centrality).

    \item The target value is reached. \code{rel\_tol} is negative and
    the primal objective \eqn{p} is less than the \code{target}.  
  }
}

\value{A \code{list}-object with the following elements:
  
  \item{x}{Solution to the primal problem.}
  \item{z}{Solution to the dual problem.}
  \item{iter}{Number of iterations performed.}
  \item{hist}{see \code{out_mode} in \code{\link{SocpControl}}.}
  \item{convergence}{A logical code. \code{TRUE} indicates successful
    convergence.}
  \item{info}{A numerical code. It indicates if the convergence was
    successful.}
  \item{message}{A character string giving any additional information
  returned by the optimiser.}
}

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
  in the Rmetrics-Project on R-Forge. In contrast to the former
  implementation, allowance is made for specifying more than one cone
  constraint. 
}

\seealso{
  \code{\link{SocpPhase1}}, \code{\link{SocpPhase2}}, \code{\link{SocpControl}}
}

\keyword{optimize}
