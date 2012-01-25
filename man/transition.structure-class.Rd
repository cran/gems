\name{transition.structure-class}
\Rdversion{1.1}
\docType{class}
\alias{transition.structure-class}
\alias{[[,transition.structure-method}
\alias{[[<-,transition.structure-method}
\alias{print,transition.structure-method}

\title{Class \code{"transition.structure"}}
\description{
This S4 class provides a structure to specify different characteristics
of transitions, such as transition functions functions, parameters or parameter covariances. 
}
\section{Objects from the Class}{
Objects are created by calls to the functions \code{generateHazardMatrix}, \code{generateParameterMatrix}, \code{generateParameterCovarianceMatrix}.
}
\section{Slots}{
  \describe{
    \item{\code{states.number}:}{Object of class \code{"numeric"}: number of states}
    \item{\code{list.matrix}:}{Object of class \code{"matrix"}: a \code{list} with two dimensions, where list element [i,j] correspond to transitions from i to j}
  }
}
\section{Methods}{
  \describe{
    \item{[[}{\code{signature(x = "transition.structure")}: extract subsets of slot \code{list.matrix} of \code{transition.structure} }
    \item{[[<-}{\code{signature(x = "transition.structure")}: replace subsets of slot \code{list.matrix} of \code{transition.structure} }
    \item{print}{\code{signature(x = "transition.structure")}}: prints the slot \code{list.matrix} of the \code{transition.structure}
   }
}
\author{
  Luisa Salazar Vizcaya, Nello Blaser, Thomas Gsponer
}

\seealso{
  \code{\link{generateHazardMatrix}}, \code{\link{generateParameterMatrix}}, \code{\link{generateParameterCovarianceMatrix}}
}
\examples{
  showClass("transition.structure")
}
\keyword{classes}
