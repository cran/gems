\name{generateParameterCovarianceMatrix}
\alias{generateParameterCovarianceMatrix}

\title{
  generate a template for parameter covariances
}
\description{
  This function simplifies generating the matrix of parameter covariances from a matrix of mean parameters.
}
\usage{
  generateParameterCovarianceMatrix(mu)
}

\arguments{
  \item{mu}{
  a \code{transition.structure} of dimension \eqn{N \times N}{N x N}, whose components \code{list} the mean values for the parameters in the \code{transitionFunction}.
}
}

\value{
  a \code{transition.structure} of dimension \eqn{N \times N}{N x N} of covariance matrices for the \code{parameters}.
}

\author{
  Luisa Salazar Vizcaya, Nello Blaser, Thomas Gsponer
}

\seealso{
 
  \code{\link{transition.structure-class}},
  \code{\link{generateParameterMatrix}},
  \code{\link{simulateCohort}}
}

\keyword{utilities}