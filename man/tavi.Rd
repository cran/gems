\name{tavi}
\alias{tavi}
\docType{data}
\title{
  tavi data set
}
\description{
The simulated data set for each patient contains data for kidney injuries, bleeding complications and the combined endpoint of stroke or death. The data was simulated from the original data following the steps described in the package vignette. 
}
\usage{data(tavi)}
\format{
  A data frame with 194 observations on the following 7 variables.
  \describe{
    \item{\code{id}}{a character vector that contains the patient id's} 
    \item{\code{kidney}}{a numeric vector; indicator variable that show if an event has occurred}
    \item{\code{kidney.dur}}{a numeric vector; times at which the events occurred or the patients were censored}
    \item{\code{bleeding}}{a numeric vector; indicator variable that show if an event has occurred}
    \item{\code{bleeding.dur}}{a numeric vector; times at which the events occurred or the patients were censored}
    \item{\code{death}}{a numeric vector; indicator variable that show if an event has occurred}
    \item{\code{death.dur}}{a numeric vector; times at which the events occurred or the patients were censored}
  }
}
\examples{
head(data(tavi))
}
\keyword{datasets}
