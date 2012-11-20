\name{BookEx}
\alias{BookEx}
\alias{listEx}
\alias{showEx}
\alias{saveEx}
\alias{editEx}
\alias{runEx}

\title{Utility functions for handling book examples}

\description{
  Utility functions for returning a list of the included examples and
  displaying, executing, saving and editing the example codes are
  provided.
}

\usage{
listEx()
showEx(Example)
saveEx(Example)
editEx(Example, ...)
runEx(Example, ...)
}

\arguments{
  \item{Example}{\code{Characater}, the name of the example as contained
    in \code{listEx()}.
  }
  \item{...}{Ellipsis argument. See details.}
}

\details{
  The ellipsis arguments in the function \code{editEx()} are passed down
  to the function \code{file.edit()}. If the option \code{editor} is
  unset and/or a different editor shall be employed for openening the
  example code, then the ellipsis argument can be utilised by
  \code{editor = "foo"}, wherey \code{foo} is the name of the editor to
  be used.\cr
  The ellipsis arguments in the function \code{runEx()} are passed down
  to the function \code{source()}. 
}

\value{
  \code{listEx}\cr
  Returns a character vector of the examples' names.\cr
  \code{showEx}\cr
  Returns the example of of \code{Example} to the console.\cr
  \code{saveEx}\cr
  Returns a \code{logical} whether the saving of the R code example into
  the working directory was successful.\cr
  \code{editEx}\cr
  Opens a copy of the example code in an editor.\cr
  \code{runEx}\cr
  Executes the example code.
}

\author{
  Bernhard Pfaff
}

\seealso{
  \code{\link{file.edit}}, \code{\link{source}} 
}

\examples{  
\dontrun{
listEx()
showEx(Example = "Part1Chapter3Ex2")
saveEx(Example = "Part1Chapter3Ex2")
runEx(Example = "Part1Chapter3Ex2", echo = TRUE)
editEx(Example = listEx()[1], editor = "emacs")
}
}

\keyword{IO}

