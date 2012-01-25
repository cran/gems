\name{generateHazardMatrix}
\alias{generateHazardMatrix}

\title{
  generate template for transition functions
}
\description{
  This function simplifies generating the matrix of transition functions.
}

\usage{
generateHazardMatrix(statesNumber)
}

\arguments{
  \item{statesNumber}{
  the number of states to be considered. 
}
}

\value{
  a \code{transition.structure} of dimension \eqn{N \times N}{N x N}, where \eqn{N} is the number of states and with value "impossible" for all potential transitions. 
}

\author{
  Luisa Salazar Vizcaya, Nello Blaser, Thomas Gsponer
}

\seealso{
  \code{\link{transition.structure-class}},
   \code{\link{simulateCohort}}
}

\keyword{utilities}