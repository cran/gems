\name{generateParameterMatrix}
\alias{generateParameterMatrix}

\title{
  generate a template for mean parameters
}
\description{
  This function simplifies generating the matrix of mean parameters from a matrix of transition functions. 
}
\usage{
  generateParameterMatrix(hf)
}

\arguments{
  \item{hf}{
  a \code{transition.structure} of dimension \eqn{N \times N}{N x N}, where \eqn{N} is the number of states.
}
}

\value{
  a \code{transition.structure} of dimension \eqn{N \times N}{N x N}, whose components are \code{lists} of the right length for the parameters in the corresponding hazard function \code{hf}. 
}

\author{
  Luisa Salazar Vizcaya, Nello Blaser, Thomas Gsponer
}

\seealso{
  \code{\link{transition.structure-class}},
  \code{\link{simulateCohort}}
}

\keyword{utilities}