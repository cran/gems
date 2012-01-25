\name{simulateCohort}
\alias{simulateCohort}

\title{
  Simulate cohort
}
\description{
  Simulates a cohort of patients from a set of functions associated to each possible transition in a multistate model. The multistate model is not required to be a Markov model and may take the history of previous events into account. In the basic version, it allows to simulate from transition-specific hazard function, whose parameters are multivariable normally distributed.
  For each state, all transition-specific hazard functions and their parameters need to be specified. For simulating one transition, all possible event times are simulated and the minimum is chosen. Then simulation continues from the corresponding state until an absorbing state of time \code{to} is reached. 
}
\usage{
simulateCohort(transitionFunctions, 
  parameters, 
  cohortSize = 1000, 
  parameterCovariances = generateParameterCovarianceMatrix(parameters), 
  timeToTransition = array(FALSE, dim = dim(transitionFunctions@list.matrix)),
  baseline = matrix(NA,nrow = cohortSize), 
  initialState=rep(1, cohortSize),
  absorbing=transitionFunctions@states.number,
  to = 100,
  report.every = 100, 
  sampler.steps = 1000)
}

\arguments{
  \item{transitionFunctions}{
  a \code{transition.structure} of dimension \eqn{N \times N}{N x N} that contains the hazard functions
  %history: vector, each element corresponds to a transition number, 0 if transition does not occur}
}
  \item{parameters}{
  a \code{transition.structure} of dimension \eqn{N \times N}{N x N} that contains the parameters
}
  \item{cohortSize}{
  a \code{numeric} indicating the number of patients to be simulated.
}
  \item{parameterCovariances}{
  a \code{transition.structure} of dimension \eqn{N \times N}{N x N} of covariance matrices for the \code{parameters}. 
}
  \item{timeToTransition}{
  a \code{logical} \code{matrix}; TRUE for all transitions whose \code{transitionFunction} is specified as the time until transition instead of as a hazard function or as a \code{character}.
}
  \item{baseline}{
  a \code{matrix} or \code{data.frame} of dimension \eqn{cohortSize \times M}{cohortSize x M} with \eqn{M} baseline characteristics of subjects to be simulated.
}
  \item{initialState}{
  a \code{numeric} of length \code{cohortSize} with the initial state for each subject simulated.
}
  \item{absorbing}{
  a \code{numeric} containing all absorbing states. 
}
  \item{to}{
  final time of the simulation.
}
  \item{report.every}{
  a \code{numeric} to check progress of simulation.
  }
  \item{sampler.steps}{
  a \code{numeric} indicating number of steps for discretization of hazard functions
  }
}

\details{
  The \code{transitionFunctions} contains hazard functions or time to event function associated to each possible transition. The elements of this \code{list} can be either expressed as an explicit R \code{function} or as a \code{character} ("impossible", "Weibull", "multWeibull", "exponential") in order to express impossible transitions or parametric forms for the distributions of time to event. If the functions should depend on time, baseline characteristics or be \emph{history-dependent}, the function arguments \emph{t}, \emph{bl} or \emph{history} can be used. Time \emph{t} refers to the time since entry into the current state. For the time since the initial state, use \code{t+sum(history)}.

  The components of the \code{parameters} argument \code{list} the mean values for the parameters in the \code{transitionFunction}. If the corresponding \code{transitionFunction} is a \code{function}, the parameters should appear in the same order as in the \code{function}, leaving out \emph{t}, \emph{bl} and \emph{history}. If the corresponding \code{transitionFunction} is the \code{character} "Weibull", the first argument is the shape and the second one the scale. If the corresponding \code{transitionFunction} is the \code{character} "multWeibull", specify weights, shapes, scales in this order.
  
  Note that when using the \code{parameterCovariances} argument it is the users responsibility to ensure that the functions are parametrized such that \code{parameters} for each transition are multivariate normally distributed and mutually independent.
}

\value{
  an object of class \code{"ArtCohort"} with \code{time.to.state} slot of dimension \eqn{cohortSize \times N}{cohortSize x N} with entry times for each patient into each of the states.
}

\author{
  Luisa Salazar Vizcaya, Nello Blaser, Thomas Gsponer
}

\examples{
# Here is an example model with 3 states and 2 possible transitions. 
 
# number of states in the model
statesNumber <- 3 

# cohort size
cohortSize <- 100

# specification of hazard functions
hazardf <- generateHazardMatrix(statesNumber)
hazardf[[1,2]] <- function(t, r1, r2) 
{
  ifelse(t<=2, r1 , r2)
}
hazardf[[2,3]] <- "Weibull" 

# list of parameters for the hazard functions
mu <- generateParameterMatrix(hazardf) 
mu[[1,2]] <- list(0.33,  0.03) # r1, r2 
mu[[2,3]] <- list(1,0.84) # shape, scale

# time
maxTime <- 10

# simulate the cohort
cohort <- simulateCohort(
  transitionFunctions = hazardf,
  parameters = mu,
  cohortSize = cohortSize,
  to=maxTime)

# output
head(cohort)

# transition probability
tr <- transitionProbabilities(cohort, times=seq(0,4,.1))
plot(tr, ci=FALSE)

# cumulative incidence
inc <- cumulativeIncidence(cohort, times=seq(0,4,.1))
plot(inc, ci=FALSE, states=c(2,3))
}

\seealso{
  \code{\link{generateHazardMatrix}}, \code{\link{generateParameterMatrix}}, \code{\link{generateParameterCovarianceMatrix}}, \code{\link{ArtCohort-class}}, \code{\link{transitionProbabilities}}, \code{\link{cumulativeIncidence}}
}

\keyword{main function}