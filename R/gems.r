########################################################################
## This program is Open Source Software: you can redistribute it      ##
## and/or modify it under the terms of the GNU General Public License ##
## as published by the Free Software Foundation, either version 3 of  ##
## the License, or (at your option) any later version.                ##
##                                                                    ##
## This program is distributed in the hope that it will be useful,    ##
## but WITHOUT ANY WARRANTY; without even the implied warranty of     ##
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  ##
## General Public License for more details.                           ##
##                                                                    ##
## You should have received a copy of the GNU General Public License  ##
## along with this program. If not, see http://www.gnu.org/licenses/. ##
########################################################################

# classes

setClass("transition.structure", representation(states.number =  "numeric", list.matrix = "matrix") )

setClass("FuncInteractions")
setClass("FuncParametersDist", representation(func =  "function", mu = "list"), contains = "FuncInteractions")

setClass("ArtCohort",
         representation(
           states.number =  "numeric",
           size =  "numeric",
           baseline = "matrix",
           follow.up = "numeric",
           parameters= "transition.structure",
           parametersCovariances = "transition.structure",
           timeToTransition = "matrix",
           transitionFunctions = "transition.structure",
           time.to.state =  "data.frame"
           ))

setClass("PosteriorProbabilities", representation(states =  "character", times =  "numeric", probabilities = "matrix",
                                                  lower= "matrix", upper = "matrix", type="character"))


# methods
setMethod("[[", "transition.structure", function(x, i, j){
  x@list.matrix[[i,j]]
})
setMethod("[[<-", "transition.structure", function(x, i, j, value){
  x@list.matrix[[i,j]] <- value
  return(x)
})
setMethod("[", "PosteriorProbabilities", function(x, i, j){
  x@probabilities[i,j]
})
setMethod("[", "ArtCohort", function(x, i, j){
  x@time.to.state[i,j]
})
possibleTransitions = function(object) 0
implicitGeneric("possibleTransitions")
setGeneric("possibleTransitions")
setMethod( "possibleTransitions", "transition.structure", function(object){
  aux1 <- auxcounter(object@states.number)
  aux1[aux1 >0] <- 1
  aux2 <- matrix(0,object@states.number,object@states.number)
  aux2[object@list.matrix == "impossible"] = 1
  aux2[is.null(object@list.matrix)] = 1
  aux3 <-  aux1-aux2
  aux4 <- which(aux3 ==1, arr.ind = TRUE)
  colnames(aux4) <-  c("From", "To")
  rownames(aux4) <-  character()
  aux4
})
setMethod("plot", "PosteriorProbabilities", function(x, ci=FALSE, main = paste(x@type, "after starting in State", x@states[1], "at time 0"), states=1:dim(x@probabilities)[2],
                                                     lwd=c(2,2), col=c('blue','green3'), lty=c(1,2), xlab="Time", ylab="Probability", ...){
  if (ci){
    plotPrevalence(x@times, x@probabilities, x@states, lower=x@lower, upper=x@upper, main=main, states=states,
                   lwd=lwd, col=col, lty=lty, xlab=xlab, ylab=ylab, ...)
    if (sum(complete.cases(x@lower))==0) {
      warning("Too few simulations for prediction intervals")
    }
    par(mfrow=c(1,1))
  }
  else {
    plotPrevalence(x@times, x@probabilities, x@states, lower=NA, upper=NA, main=main, states=states,
                   lwd=lwd, col=col, lty=lty, xlab=xlab, ylab=ylab, ...)
    par(mfrow=c(1,1))
  }
})

setMethod( "update", "ArtCohort", function(object, newsize, addbaseline=matrix(NA, nrow = newsize - object@size), newInitialStates=rep(1, newsize - object@size)){
  aux1 =
    simulateCohort(
      transitionFunction = object@transitionFunctions,
      parameters = object@parameters,
      cohortSize = newsize - object@size ,
      parameterCovariance = object@parametersCovariances,
      timeToTransition= object@timeToTransition,
      baseline = addbaseline,
      initialState=newInitialStates,
      to=5
      )

  newcohort <- object
  newcohort@baseline <- rbind(object@baseline, addbaseline)
  newcohort@size <- newsize
  aux2 <- rbind(object@time.to.state, aux1@time.to.state)
  rownames(aux2) <- paste("Patient", 1:newsize)
  newcohort@time.to.state <-  aux2
  newcohort
})

setMethod("head", "ArtCohort", function(x, ...){
  head(x@time.to.state, ...)
})

setMethod("tail", "ArtCohort", function(x, ...){
  tail(x@time.to.state, ...)
})

setMethod("head", "PosteriorProbabilities", function(x, ...){
  head(x@probabilities, ...)
})

setMethod("tail", "PosteriorProbabilities", function(x, ...){
  tail(x@probabilities, ...)
})

setMethod("print", "transition.structure", function(x){
  mat <- x@list.matrix
  mchar <- lapply(mat, function(x) {
    y <- x
    if (is.function(x)) {
      y <- paste("function(", paste(names(formals(x)), collapse=", "), ")", sep="")
      }
    return(y)
    })
  dim(mchar) <- dim(mat)
  dimnames(mchar) <- dimnames(mat)
  print(mchar)
})

setMethod( "summary", "ArtCohort", function(object){
  summ <- matrix(list(
    (object@states.number),
    (object@size),
    (object@baseline),
    (object@follow.up),
    (object@parameters),
    (object@parametersCovariances),
    (object@timeToTransition),
    (object@transitionFunctions),
    (object@time.to.state)),9,2)

  aux2  <-  getSlots("ArtCohort")
  summ[,2] <-  noquote(aux2)
  colnames(summ) <- c("Mode", "Class")
  rownames(summ) <- names(aux2)

  print(noquote(matrix( c(class(object), mode(object)),ncol =  2, nrow =  1,dimnames =  list("",c("Class", "Mode")))))

  summ.all <-  list()
  summ.all$general <-  summ
  summ.all$events  <- noquote(summary(object@time.to.state)[,2:object@states.number])

  aux3 <- unlist(apply(object@time.to.state, 2, function(x) length(which(is.na(x) == 0))))

  message("
Number of patients entering each state: ")
  print(aux3)
  message("
object and events: ")
  summ.all
})

# functions
auxcounter <-
  function (statesNumber)
{
### "auxcounter" associates an index to each position in the transition matrix:
#   Input: number of states
#   Output: matrix M with transitions numbers, where element M_ij is the
#     transition number of the transition from State i to state j.
  m = matrix(0, ncol = statesNumber, nrow = statesNumber)
  for (i in 1:(statesNumber - 1)) {
    m[i, (i + 1):statesNumber] <- 1:(statesNumber - i) +
      (i - 1) * (statesNumber - i/2)
  }
  dimnames(m) <- list(from = paste("State", 1:statesNumber),
                      to = paste("State", 1:statesNumber))
  return(m)
}

auxposition <-
  function (matrix, value)
{
    ### "auxposition" finds the position of a value in a matrix:
    #   Input: a matrix and a value
    #   Output: position of (1st occurence of) value in matrix
  which(matrix == value, arr.ind = TRUE)[1, ]
}

### The function "createCohorts" is the main internal function.
createCohorts <-
function (hazardf, statesNumber, cohortSize, mu, sigma = matrix(0,
                                                   nrow = length(unlist(mu)), ncol = length(unlist(mu))), historyl = FALSE,
          startingStates = rep(1, cohortSize), absorbing=cohortSize, impossible = NULL, fixpar = NULL,
          direct = NULL, bl0 = matrix(0, nrow = cohortSize), to = 100)
{
  hazardf <- t(hazardf)[t(auxcounter(statesNumber)) > 0]

  if (length(hazardf) < (statesNumber * (statesNumber - 1)/2)) {
    length(hazardf) <- statesNumber * (statesNumber - 1)/2
  }
  for (tr in 1:length(hazardf)) {
    if (length(hazardf[[tr]]) > 0) {
      if (!is.function(hazardf[[tr]])) {
        if (!hazardf[[tr]] %in% c("Weibull", "multWeibull", "Exponential", "impossible") &&
            !is.na(hazardf[[tr]]))
          try(hazardf[[tr]] <- as.function(hazardf[[tr]]))
      }
    }
  }

  mu <- as.list(unlist(t(mu), recursive = FALSE))

  if (length(mu) > 0) {
    mu <- lapply(mu, as.numeric)
  }

  if (dim(sigma)[1] != length(unlist(mu)))
    stop("size of parameters and parameterCovariances are inconsistent")

  if (length(startingStates) == 1) rep(startingStates, cohortSize)

  if (cohortSize>1) {try(bl0 <- as.matrix(bl0))}
  else if (is.null(dim(bl0))) {try(bl0 <- t(as.matrix(bl0)))}
  if (nrow(bl0) != cohortSize)
    warning("baseline is not of the right dimension")

  for (i in 1:length(hazardf)) {
    if (is.element(i, impossible)) {
            if (historyl == FALSE)
              hazardf[[i]] = function(t, bl) 99 * to
            else hazardf[[i]] = function(t, history, bl) 99 *
              to
          }
  }
  parametric = numeric()
  k = 0
  for (i in 1:max(auxcounter(statesNumber))) {
    if (is.character(hazardf[[i]])) {
      k = k + 1
      parametric[k] = i
      if (hazardf[[i]] %in% c("Weibull", "multWeibull", "Exponential")) {
        hazardf[[i]] <- switch(hazardf[[i]],
                               Weibull=function(t, shape, scale, history, bl) rweibull(t, shape, scale),
                               multWeibull=function(t, w, shapes, scales, history, bl) multPar(t, w, shapes, scales),
                               Exponential=function(t, rate, history, bl) rexp(t, rate))
        if (!historyl) formals(hazardf[[i]])$history <- NULL
      }
      else {
        stop(paste("Transition function for transition", i, "not recognized", sep=" "))
      }
    }
  }
  allFunctions = mainFunctions(statesNumber = statesNumber,
    Mu = mu, sigma = sigma, cohortSize = cohortSize, history = historyl,
    hazardf, impossible = c(impossible, fixpar, direct))
  parametric = sort(c(parametric, impossible, direct))
  cohorts <- sapply(1:cohortSize, function(i) historical(gf = allFunctions[[i]],
                                                         statesNumber = statesNumber, parametric = parametric,
                                                         historyl = historyl, startingState = startingStates[i],
                                                         absorbing=absorbing, bl = bl0[i, ], to = to)[[1]])
  dimnames(cohorts) <- list(paste("State", 1:statesNumber),
                            paste("Patient", 1:cohortSize))
  return(cohorts)
}
data_prep <-
  function (progressionObject, maxtime)
{
  states_number <- dim(progressionObject)[1]
  cohortSize <- dim(progressionObject)[2]
  times <- t(progressionObject)
  status <- array(TRUE, dim = dim(times))
  status[is.na(times)] <- FALSE
  times[times>=maxtime]<-maxtime
  status[times >= maxtime] <- FALSE
  times[status == FALSE] <- maxtime + 1
  tmat <- auxcounter(states_number)
  tmat[tmat==0] <- NA
  prepData <- msprep(times, status, trans = tmat,
                     start = list(state = 1, time = 0))
  return(prepData)
}
fold <-
  function (vector, object)
{
  kk <- 0
  lc <- 0
  newobject <- list()
  for (i in 1:length(object)) {
    k <- 0
    newobject[[i]] <- list()
    for (j in 1:(length(object[[i]]@mu))) {
      l <- length((object[[i]]@mu)[[j]])
      newobject[[i]][[j]] = c(vector[(1 + k + kk):(k +
                      kk + l)])
      k <- k + l
      lc <- lc + l
    }
    kk <- lc
  }
  return(newobject)
}
generateHazardMatrix <-
  function (statesNumber)
{
  hf <- list()
  length(hf) <- statesNumber * (statesNumber - 1)/2
  for (i in 1:length(hf)) {
    hf[[i]] <- "impossible"
  }
  tmat <- auxcounter(statesNumber)
  for (trans in 1:(statesNumber * (statesNumber - 1)/2)) {
    attributes(hf)$names[trans] <- paste("Transition", trans,
                                         "- from state", auxposition(tmat, trans)[1], "to state",
                                         auxposition(tmat, trans)[2])
  }
  hfmat <- matrix(list(), nrow = statesNumber, ncol = statesNumber)
  for (i in 1:(statesNumber * (statesNumber - 1)/2)) {
    hfmat[tmat == i][[1]] <- hf[[i]]
  }
  dimnames(hfmat) <- list(from = paste("State", 1:statesNumber),
                          to = paste("State", 1:statesNumber))

  hfclass <- new("transition.structure")
  hfclass@states.number <- dim(hfmat)[1]
  hfclass@list.matrix <- hfmat
  return(hfclass)
}

generateParameterCovarianceMatrix <-
  function (mu)
{
  mu  =  mu@list.matrix
  SigmaMat <- matrix(list(), nrow = dim(mu)[1], ncol = dim(mu)[2])
  dimnames(SigmaMat) <- list(from = paste("State", 1:dim(mu)[1]),
                             to = paste("State", 1:dim(mu)[1]))
  lengthMu <- numeric(length(mu))
  for (trans in 1:length(mu)) {
    lengthMu[trans] <- length(unlist(mu[[trans]]))
    if (lengthMu[trans] > 0) {
      SigmaMat[[trans]] <- matrix(0, nrow = lengthMu[trans],
                                  ncol = lengthMu[trans])
    }
  }
  Sigmaclass <- new("transition.structure")
  Sigmaclass@states.number <- dim(SigmaMat)[1]
  Sigmaclass@list.matrix <- SigmaMat
  return(Sigmaclass)
}
generateParameterMatrix <-
  function (hf)
{
  hfMat <- hf@list.matrix
  tmat <- auxcounter(dim(hfMat)[1])
  hf <- t(hfMat)[t(tmat > 0)]
  MuMat <- matrix(list(), nrow = dim(hfMat)[1], ncol = dim(hfMat)[2])
  dimnames(MuMat) <- list(from = paste("State", 1:dim(hfMat)[1]),
                          to = paste("State", 1:dim(hfMat)[1]))
  lengthMu <- numeric(length(hf))
  for (trans in 1:length(hf)) {
    if (class(hf[[trans]]) == "character") {
      lengthMu[trans] <- switch(hf[[trans]],
                                Weibull=2,
                                multWeibull=3,
                                Exponential=1,
                                impossible=0)
    }
    else if (class(hf[[trans]]) == "function") {
      if (length(formals(hf[[trans]])) > 0) {
        n1 <- attributes(formals(hf[[trans]]))$names
        n2 <- n1[!(n1 %in% c("t", "bl", "history", "process"))]
        lengthMu[trans] <- length(n2)
      }
    }
  }
  for (trans in 1:length(hf)) {
    if (lengthMu[trans] > 0) {
      MuMat[tmat == trans][[1]] <- as.list(numeric(lengthMu[trans]))
    }
  }
  Muclass <- new("transition.structure")
  Muclass@states.number <- dim(MuMat)[1]
  Muclass@list.matrix <- MuMat
  return(Muclass)
}

historical <- function(gf, statesNumber, parametric, historyl, startingState, absorbing, bl, to){
  if (historyl) histH(gf, statesNumber, parametric, startingState, absorbing, bl, to)
  else histNoH(gf, statesNumber, parametric, startingState, absorbing, bl, to)
}

histNoH <- function(gf, statesNumber, parametric, startingState, absorbing, bl, to){
  possible <- auxcounter(statesNumber)
  # all possible transition times
  npar <- (1:length(gf))[!1:length(gf)%in%parametric]
  time0 <- rep(NA,length(gf))
  time0[parametric] <- unlist(lapply(gf[parametric], function(ff) samplerP(function(t) ff(t, bl=bl), n=1)))
  time0[npar] <- unlist(lapply(gf[npar], function(ff) sampler(n=1, function(t) ff(t, bl=bl))))
  #times in the transition matrix
  time = matrix(ncol = statesNumber, nrow =  statesNumber)
  time[unlist(lapply(1:max(possible), function(i) which(possible == i)))] <- time0
  time[possible==0] <- 99*to

  minTimeState <- apply(time, 1, function(x) c(min(x), which.min(x)))
  kk = startingState
  path = rep(NA,statesNumber)
  path[startingState]<-0
  aux =  startingState
  while(aux %in% setdiff(startingState:statesNumber, absorbing) && sum(path, na.rm=TRUE)<=to){
    aux = minTimeState[2,aux]
    path[aux] = minTimeState[1,kk] + sum(path, na.rm=TRUE)
    kk =  aux
  }
  return(list(path,NULL))
}

histH <- function(gf, statesNumber, parametric, startingState, absorbing, bl, to){
  possible <- auxcounter(statesNumber)
  formerhistory <- array(0,dim =  c(max(possible), 4))
  eventTimes   <-  matrix(NA, ncol = statesNumber ,nrow= 1)
  aux = possible[startingState,]; aux = aux[aux !=0]
  par = match(intersect(aux,parametric),aux)
  f =  gf[aux]; rm(aux)
  npar <- (1:length(f))[!1:length(f)%in%par]
  time = numeric(length(f))
  time[par] <- unlist(lapply(f[par], function(ff) samplerP(function(t) ff(t, history=0, bl=bl), n=1)))
  time[npar] <- unlist(lapply(f[npar], function(ff) sampler(n=1, function(t) ff(t, history=0, bl=bl),to=to)))
  #########################################
  mt = min(time)
  wmin = which.min(time)
  first <- possible[startingState, startingState+wmin]
  k = auxposition(possible, max(first))[2]  #current state of the patient
  formerhistory[possible[startingState,startingState+1]+wmin-1,1:4] = c(first, mt, mt, k)
  lim = mt
  ######## loop over states...
  while(k %in% setdiff(startingState:statesNumber, absorbing) && lim <= to){
    aux = possible[auxposition(possible, max(na.rm = TRUE,formerhistory[,1]))[2],]; aux =  aux[aux !=0]
    f =  gf[aux]
    par = match(intersect(aux,parametric),aux)  #those index of f (possible transition hazards corespond to parametric)
    npar <- (1:length(f))[!1:length(f)%in%par]

    timesp  =  numeric(length(f))
    timesp[par] <- unlist(lapply(f[par], function(ff) samplerP(function(t) ff(t, history=formerhistory[,2], bl=bl), n=1)))
    timesp[npar] <- unlist(lapply(f[npar], function(ff) sampler(n=1, function(t) ff(t, history=formerhistory[,2], bl=bl),to=to)))

    mt = min(timesp)
    wmin = which.min(timesp)
    auxmin <- aux[wmin]

    lim = mt + formerhistory[max(formerhistory[,1]) ,3]
    k = auxposition(possible, max(auxmin))[2]

    formerhistory[auxmin,1:4] = c(auxmin, mt, lim, k)
  }
  for(st in startingState:statesNumber){
    if(length(which(formerhistory[,4]==st)) == 0)
      eventTimes[1,st] = ifelse(st==startingState, 0, NA)
    else  eventTimes[1,st] = formerhistory[which(formerhistory[,4]==st)[1],3]
  }
  return(list(eventTimes, NULL))
}

mainFunctions <-
  function (statesNumber, Mu, sigma, cohortSize, history, functions,
            impossible)
{
  objects = setFunctions(statesNumber = statesNumber, Mu = Mu,
    history = history, functions, impossible = impossible)
  return(simFunctions(objects, covariances = sigma, history = history,
                      statesNumber = statesNumber, impossible = impossible,
                      cohortSize = cohortSize))
}
msmDataPrep <-
  function (mstateDataPrep)
{
  a <- mstateDataPrep
  cohortSize <- max(a$id)
  initState <- sapply(1:cohortSize, function(ID) min(a$from[a$id==ID]))
  unsortedId <- c(a$id[a$status == 1], seq(1, cohortSize))
  unsortedTime <- c(a$Tstop[a$status == 1], rep(0, cohortSize))
  unsortedState <- c(a$to[a$status == 1], initState)
  unsortedData <- data.frame(unsortedId, unsortedTime, unsortedState)
  sortedData <- unsortedData[order(unsortedId, unsortedTime),
                             ]
  msmData <- data.frame(id = sortedData$unsortedId, time = sortedData$unsortedTime,
                        state = sortedData$unsortedState)
}
multPar <-
  function (t, w, shape, scale)
{
  ii <- apply(rmultinom(t, 1, c(w, 1 - sum(w))), 2, which.max)
  mW <- rweibull(t, shape[ii], scale[ii])
  attr(mW, "weights") <- w
  attr(mW, "index") <- ii
  return(mW)
}
plotPrevalence <-
  function (times, prevalence, stateNames, lower, upper, typ = "separate", main = "State-wise probability over time", states,
            lwd=c(2,2), col=c('blue','green3'), lty=c(1,2), xlab="Time", ylab="Probability" , ...)
{
  if (typ == "separate") {
    if (length(col)==1) col <- rep(col, 2)
    if (length(lty)==1) lty <- rep(col, 2)
    if (length(lwd)==1) lwd <- rep(lwd, 2)
    par(mfrow = c(2, ceiling(length(states)/2)))
    if (length(states)==1) par(mfrow = c(1,1))
    par(oma = c(0, 0, 2, 0))
    for (i in states) {
      plot(times, prevalence[, i], type = "l", col = col[1],
           lwd = lwd, ylim = c(0, 1), ylab = ylab,
           xlab = xlab, main = stateNames[[i]], ...)
      try(lines(times, lower[,i], col=col[2], lty=lty[2]), silent=TRUE)
      try(lines(times, upper[,i], col=col[2], lty=lty[2]), silent=TRUE)
    }
    mtext(main, outer = TRUE, cex = 1.5)
  }
  else if (typ == "stacked") {
    par(mfrow = c(1, 1))
    totalPrev <- prevalence[, 1]
    plot(times, totalPrev, col = 1, type = "l", ylim = c(0,
                                                  1), ylab = "Probability", xlab = "Time")
    for (i in 2:dim(prevalence)[2]) {
      totalPrev <- totalPrev + prevalence[, i]
      lines(times, totalPrev, col = i)
    }
  }
}

posteriorProbabilities = function(object, times, M=100, stateNames = paste("State", as.list(1:dim(cohorts)[1]))) {
  if (class(object)=="ArtCohort") cohorts <- t(object@time.to.state)
  else cohorts <- t(object)
  statesNumber <- dim(cohorts)[1]

  # Prediction interval:
  if (dim(cohorts)[2]>=1000) {
    dd <- data.frame(t(cohorts))
    dd$cc <- c(rep(1:M, times=floor(dim(cohorts)[2])/M) , 1:(dim(cohorts)[2]%%M+1))[1:dim(cohorts)[2]]
    dd <- split(dd,dd$cc)
    prev <- list()
    ppp <- NULL
    for (ciind in 1:M){
      cohorts <- t(dd[[ciind]][,1:statesNumber])
      if (ncol(cohorts)==0) {
        prev[[ciind]] <- matrix(NA, ncol = nrow(cohorts), nrow = length(times))
      }
      if (ncol(cohorts)==1) {
        prev[[ciind]] <- matrix(0, ncol = nrow(cohorts), nrow = length(times))
        for (tInd in 1:(nrow(cohorts)-1)) {
          prev[[ciind]][times < min(c(max(times)+1,cohorts[(tInd+1):nrow(cohorts)]), na.rm=TRUE) & times >= cohorts[tInd] , tInd] <- 1
        }
        prev[[ciind]][times >= cohorts[nrow(cohorts)] , nrow(cohorts)] <- 1
      }
      if (ncol(cohorts)>=2) {
        prep1 <- data_prep(cohorts, max(times))
        prep1 <- prep1[!(prep1$Tstart == prep1$Tstop), ]
        prepData <- msmDataPrep(prep1)
        prev[[ciind]] <- prevalence(prepData, times, dim(cohorts)[1])
      }
      ppp <- rbind(ppp,prev[[ciind]])
    }
    ddd <- data.frame(ppp)
    ddd$times <- times

    prev <- as.matrix(ddply(ddd, .(times), function(x) colMeans(x)))[, 1:statesNumber]
    #prev <- as.matrix(ddply(ddd, .(times), function(x) apply(x, 2, function(y){quantile(y, .5)})))[, 1:statesNumber]
    lower <- as.matrix(ddply(ddd, .(times), function(x) apply(x, 2, function(y){quantile(y, .025)})))[, 1:statesNumber]
    upper <- as.matrix(ddply(ddd, .(times), function(x) apply(x, 2, function(y){quantile(y, .975)})))[, 1:statesNumber]
  }
  else {
    if (ncol(cohorts)==0) {
      prev <- matrix(NA, ncol = nrow(cohorts), nrow = length(times))
    }
    if (ncol(cohorts)==1) {
      prev <- matrix(0, ncol = nrow(cohorts), nrow = length(times))
      for (tInd in 1:(nrow(cohorts)-1)) {
        prev[times < min(c(max(times)+1,cohorts[(tInd+1):nrow(cohorts)]), na.rm=TRUE) & times >= cohorts[tInd] , tInd] <- 1
      }
      prev[times >= cohorts[nrow(cohorts)] , nrow(cohorts)] <- 1
    }
    if (ncol(cohorts)>=2) {
      prep1 <- data_prep(cohorts, max(times))
      prep1 <- prep1[!(prep1$Tstart == prep1$Tstop), ]
      prepData <- msmDataPrep(prep1)
      prev<- prevalence(prepData, times, dim(cohorts)[1])
    }
    lower <- matrix(NA, ncol = nrow(cohorts), nrow = length(times))
    upper <- matrix(NA, ncol = nrow(cohorts), nrow = length(times))
  }

  dimnames(prev) <- list(paste("Time", times), stateNames)
  dimnames(lower) <- list(paste("Time", times), stateNames)
  dimnames(upper) <- list(paste("Time", times), stateNames)

  pp = new("PosteriorProbabilities")
  pp@states = stateNames
  pp@times = times
  pp@probabilities = prev
  pp@lower = lower
  pp@upper = upper
  pp@type = "Transition probabilities"
  return(pp)
}

cumulativeIncidence = function(object, times, M=100, stateNames = paste("State", as.list(1:dim(cohorts)[1]))) {
  if (class(object)=="ArtCohort") cohorts <- t(object@time.to.state)
  else cohorts <- t(object)
  statesNumber <- dim(cohorts)[1]

  # Prediction interval:
  if (dim(cohorts)[2]>=1000) {
    dd <- data.frame(t(cohorts))
    dd$cc <- c(rep(1:M, times=floor(dim(cohorts)[2])/M) , 1:(dim(cohorts)[2]%%M+1))[1:dim(cohorts)[2]]
    dd <- split(dd,dd$cc)
    inc <- list()
    ppp <- NULL
    for (ciind in 1:M){
      cohorts <- t(dd[[ciind]][,1:statesNumber])

      if (ncol(cohorts)==0) {
        inc[[ciind]] <- matrix(NA, ncol = nrow(cohorts), nrow = length(times))
      }
      if (ncol(cohorts)==1) {
        inc[[ciind]] <- matrix(0, ncol = nrow(cohorts), nrow = length(times))
        for (tInd in 1:(nrow(cohorts))) {
          inc[[ciind]][times >= cohorts[tInd] , tInd] <- 1
        }
      }
      if (ncol(cohorts)>=2) {
        prep1 <- gems:::data_prep(cohorts, max(times))
        prep1 <- prep1[!(prep1$Tstart == prep1$Tstop), ]
        prepData <- gems:::msmDataPrep(prep1)
        inc[[ciind]] <- incidence(prepData, times, dim(cohorts)[1])
      }
      ppp <- rbind(ppp,inc[[ciind]])
    }
    ddd <- data.frame(ppp)
    ddd$times <- times

    inc <- as.matrix(ddply(ddd, .(times), function(x) colMeans(x)))[, 1:statesNumber]
    lower <- as.matrix(ddply(ddd, .(times), function(x) apply(x, 2, function(y){quantile(y, .025)})))[, 1:statesNumber]
    upper <- as.matrix(ddply(ddd, .(times), function(x) apply(x, 2, function(y){quantile(y, .975)})))[, 1:statesNumber]
  }

  else {
    if (ncol(cohorts)==0) {
      inc <- matrix(NA, ncol = nrow(cohorts), nrow = length(times))
    }
    if (ncol(cohorts)==1) {
      inc <- matrix(0, ncol = nrow(cohorts), nrow = length(times))
      for (tInd in 1:(nrow(cohorts))) {
        inc[times >= cohorts[tInd] , tInd] <- 1
      }
    }
    if (ncol(cohorts)>=2) {
      prep1 <- data_prep(cohorts, max(times))
      prep1 <- prep1[!(prep1$Tstart == prep1$Tstop), ]
      prepData <- msmDataPrep(prep1)
      inc <- incidence(prepData, times, dim(cohorts)[1])
    }
    lower <- matrix(NA, ncol = nrow(cohorts), nrow = length(times))
    upper <- matrix(NA, ncol = nrow(cohorts), nrow = length(times))
  }

  dimnames(inc) <- list(paste("Time", times), stateNames)
  dimnames(lower) <- list(paste("Time", times), stateNames)
  dimnames(upper) <- list(paste("Time", times), stateNames)

  pp = new("PosteriorProbabilities")
  pp@states = stateNames
  pp@times = times
  pp@probabilities = inc
  pp@lower = lower
  pp@upper = upper
  pp@type = "Cumulative incidence"
  return(pp)
}
transitionProbabilities <- posteriorProbabilities
prepareF <-
  function (func, known, historyl, historyp = numeric(2), blp = numeric(2))
{
  lsim = length(known)
  invisible(func)
  if (historyl == TRUE) {
    ofnh = function(t, sim = known, history = historyp, bl = blp) {
      hisaux = matrix(0, ncol = 1, nrow = length(history))
      hisaux[1:length(history), 1] = history[1:length(history)]
      l = list()
      l[[1]] = t
      l[2:(length(sim) + 1)] = sim[1:length(sim)]
      l[[length(sim) + 2]] = hisaux
      l[[length(sim) + 3]] = bl
      return(do.call(func, l))
    }
    return(ofnh)
  }
  else {
    ofnl = function(t, sim = known, bl = blp) {
      l = list(lsim + 2)
      l[[1]] = t
      l[2:(length(sim) + 1)] = sim[1:length(sim)]
      l[[length(sim) + 2]] = bl
      return(do.call(func, l))
    }
    return(ofnl)
  }
}
prevalence <-
  function (data, times, states_number)
{
  states <- array(0, dim = c(length(times), max(data$id)))
  for (i in 1:max(data$id)) {
    ti <- data$time[data$id == i]
    st <- data$state[data$id == i]
    for (j in 1:length(ti)) {
      states[times >= ti[j], i] <- st[j]
    }
  }
  prev <- matrix(0, ncol = states_number, nrow = length(times))
  for (i in 1:states_number) {
    prev[, i] <- apply((states == i), 1, sum)
  }
  prev <- prev/length(unique(data$id))
  return(prev)
}
incidence <-
  function (data, times, states_number)
{
  spl <- split(data, factor(data$state, levels=1:states_number))
  inc <- matrix(0, ncol = states_number, nrow = length(times))
  for (i in 1:states_number){
    if (dim(spl[[i]])[1]>0){
      sti <- sapply(times, function(x) x>=spl[[i]]$time)
      if (dim(spl[[i]])[1]>1) inc[,i] <- colSums(sti)
      else inc[,i] <- sti
    }
    else {
      inc[,i] <- 0
    }
  }
  inc <- inc/length(unique(data$id))
  return(inc)
}
sampler <-
  function (n, f, to = 100, length = 1000)
{
  time = seq(0, to, length = length)
  rate = f(time)
  rate[rate == Inf] <- 0
  rate[is.nan(rate)] <- 0
  rate[is.na(rate)] <- 0
  rate[rate < .Machine$double.eps] <- .Machine$double.eps
  samples = rpexp(n, rate, time)
  if (is.nan(samples)) {
    print(rate)
  }
  samples[is.nan(samples)] = 0
  return(samples)
}
samplerP <-
  function (f, n)
{
  f = match.fun(f)
  return(f(t = n))
}
setFunctions <-
  function (statesNumber, Mu, history, functions, impossible)
{
  trans = max(auxcounter(statesNumber))
  sample = list(trans)
  narg = numeric()
  k = 0
  for (i in 1:trans) {
    sample[[i]] = new("FuncParametersDist")
    sample[[i]]@func = functions[[i]]
    if (!is.element(i, impossible)) {
      if (history == TRUE)
        narg[i] = length(formals(functions[[i]])) - 2 -
          1
      else narg[i] = length(formals(functions[[i]])) -
        1 - 1
      if (i == 1)
        sample[[i]]@mu = Mu[1:narg[i]]
      else sample[[i]]@mu = Mu[(k + 1):(k + narg[i])]
      k = k + narg[i]
    }
  }
  return(sample)
}
simFunctions <-
  function (so, covariances, history, statesNumber, impossible,
            cohortSize)
{

  possible = max(auxcounter(statesNumber))
  lso = length(so)
  prov = setdiff(1:possible, impossible)
  for (i in prov) {
    if (i == prov[1]) {
      aux1 = unfold(so[[i]]@mu)
      aux2 = rep(1, length(aux1))
    }
    else {
      uf = unfold(so[[i]]@mu)
      aux1 = c(aux1, uf)
      aux2 = c(aux2, rep(i, length(uf)))
    }
  }

  if (sum(covariances !=0)){
    simpar = mvrnorm(cohortSize, aux1, covariances)
  }
  else{
    simpar <- matrix(aux1, nrow=cohortSize, ncol=length(aux1), byrow=TRUE)
  }

  simpar = rbind(aux2, simpar)
  auxso = so[prov]
  simpar1 <- lapply(1:cohortSize, function(k) fold(simpar[k+1,], auxso))
  lp = list(lso)
  ltp = list(cohortSize)
  for (j in 1:cohortSize) {
    k = 0
    for (i in 1:lso) {
      if (!is.element(i, impossible)) {
        k = k + 1
        lp[[i]] = prepareF(func = so[[i]]@func, known = simpar1[[j]][[k]],
            historyl = history)
      }
      else {
        lp[[i]] = so[[i]]@func
      }
    }
    ltp[[j]] = lp
  }
  return(ltp)
}
simulateCohort <-
  function(transitionFunctions,
           parameters,
           cohortSize=1000,
           parameterCovariances=generateParameterCovarianceMatrix(parameters),
           timeToTransition=array(FALSE, dim = dim(transitionFunctions@list.matrix)),
           baseline=matrix(NA, nrow = cohortSize),
           initialState=rep(1, cohortSize),
           absorbing=transitionFunctions@states.number,
           to=100){

    try(cohortSize <- as.integer(cohortSize))
    if (length(cohortSize) > 1)
      warning("cohortSize has length > 1 and only the first element will be used")
    cohortSize <- cohortSize[1]

    try(to <- as.numeric(to))
    if (length(to) > 1)
      warning("to has length > 1 and only the first element will be used")
    to <- to[1]

    statesNumber <- dim(transitionFunctions@list.matrix)[1]

    stopifnot(cohortSize>0,
              is.list(transitionFunctions@list.matrix),
              is.list(parameters@list.matrix),
              is.list(parameterCovariances@list.matrix),
              dim(transitionFunctions@list.matrix)[1]==dim(transitionFunctions@list.matrix)[2],
              identical(dim(transitionFunctions@list.matrix), dim(parameters@list.matrix)),
              identical(dim(parameters@list.matrix), dim(parameterCovariances@list.matrix)),
              identical(dim(transitionFunctions@list.matrix), dim(timeToTransition)),
              initialState %in% 1:statesNumber,
              absorbing %in% 1:statesNumber,
              length(initialState) %in% c(1, cohortSize)
              )

    transitionFunction  =  transitionFunctions@list.matrix
    parameterCovariance =  parameterCovariances@list.matrix
    parameters.in =  parameters@list.matrix

    historyDep <- length(grep(c("history"), transitionFunction))>0

    impossible=matrix(FALSE, nrow=dim(transitionFunction)[1], ncol=dim(transitionFunction)[1])
    fixpar=matrix(FALSE, nrow=dim(transitionFunction)[1], ncol=dim(transitionFunction)[1])

    hf <- transitionFunction
    for (trans in 1:length(transitionFunction)){
      if (length(hf[[trans]]) > 0 & is.character(hf[[trans]]))  {
        if (hf[[trans]]=="impossible"){
          impossible[[trans]] <- TRUE
        }
      }
      if (is.function(hf[[trans]])) {
        params <- names(formals(hf[[trans]]))
        historyIndex <- which(params=="history")
        blIndex <- which(params=="bl")
        timeIndex <- which(params=="t")
        if(length(params[-c(historyIndex, blIndex, timeIndex)])==0) fixpar[[trans]] <- TRUE

        if (historyDep){
          newParams <- c("t", params[-c(historyIndex, blIndex, timeIndex)], "history", "bl")
        }
        if (length(c(historyIndex, blIndex, timeIndex))>0 & !historyDep){
          newParams <- c("t", params[-c(historyIndex, blIndex, timeIndex)], "bl")
        }
        if (length(c(historyIndex, blIndex, timeIndex))==0 & !historyDep){
          newParams <- c("t", params, "bl")
        }

        fct <- function(x) {}
        forms <- vector("pairlist", length=length(newParams))
        for (par in 1:length(forms)){
          forms[[par]] <- formals(fct)$x
        }
        names(forms) <- newParams
        formals(hf[[trans]]) <- forms
        rm(fct, forms)
      }
    }

    sigma <- matrix(0, ncol=length(unlist(parameters.in)), nrow=length(unlist(parameters.in)))

    start <- 0
    for (matIndex in 1:dim(parameterCovariance)[1]) {
      leng <- dim(parameterCovariance[[matIndex]])[1]
      if (is.null(leng)) leng <- 0
      if (leng>0) {
        sigma[(start+1):(start+leng), (start+1):(start+leng)] <- parameterCovariance[[matIndex]]
      }
      start <- start+leng
    }

    impossible <- auxcounter(statesNumber)[impossible]
    fixpar <- auxcounter(statesNumber)[fixpar]
    direct <- auxcounter(statesNumber)[timeToTransition]

    cohort <- createCohorts(hazardf=hf,
                            statesNumber=statesNumber,
                            cohortSize=cohortSize,
                            mu=parameters.in,
                            sigma = sigma,
                            historyl = historyDep,
                            impossible = impossible,
                            fixpar = fixpar,
                            direct = direct,
                            bl0=baseline,
                            startingStates=initialState,
                            absorbing=absorbing,
                            to=to)
    cohort[cohort>to] <- NA

    simulatedcohort <- new("ArtCohort")
    simulatedcohort@baseline <- as.matrix(baseline)
    simulatedcohort@follow.up <- to
    simulatedcohort@parameters <- parameters
    simulatedcohort@parametersCovariances <- parameterCovariances
    simulatedcohort@timeToTransition <- timeToTransition
    simulatedcohort@transitionFunctions <- transitionFunctions
    simulatedcohort@states.number <- dim(cohort)[1]
    simulatedcohort@size <- dim(cohort)[2]
    simulatedcohort@time.to.state <- as.data.frame(t(cohort))

    return(simulatedcohort)
  }

unfold <-
  function (object)
{
  vector <- numeric()
  k <- 0
  for (i in 1:length(object)) {
    vector[(1 + k):(length(object[[i]]) + k)] <- c(object[[i]])
    k <- length(object[[i]]) + k
  }
  return(vector)
}

msprep <- function (time, status, data, trans, start, id, keep)
{
  if (!(is.matrix(time) | (is.data.frame(time)))) {
    if (!is.character(time))
      stop("argument \"time\" should be a character vector")
    if (missing(data))
      stop("missing \"data\" argument not allowed when time argument is character vector")
    startings <- which(apply(!is.na(trans), 2, sum) == 0)
    wh <- which(is.na(time))
    if (!all(wh %in% startings))
      stop("no NA's allowed in the \"time\" argument for non-starting states")
    tcols <- match(time[!is.na(time)], names(data))
    if (any(is.na(tcols)))
      stop("at least one of elements of \"time\" not in data")
    time <- matrix(NA, nrow(data), length(time))
    whcols <- (1:ncol(time))[!(1:ncol(time)) %in% wh]
    time[, whcols] <- as.matrix(data[, tcols])
  }
  if (!(is.matrix(status) | (is.data.frame(status)))) {
    if (!is.character(status))
      stop("argument \"status\" should be a character vector")
    if (missing(data))
      stop("missing \"data\" argument not allowed when status argument is character vector")
    startings <- which(apply(!is.na(trans), 2, sum) == 0)
    wh <- which(is.na(status))
    if (!all(wh %in% startings))
      stop("no NA's allowed in the \"status\" argument for non-starting states")
    dcols <- match(status[!is.na(status)], names(data))
    if (any(is.na(dcols)))
      stop("at least one of elements of \"status\" not in data")
    status <- matrix(NA, nrow(data), length(status))
    whcols <- (1:ncol(status))[!(1:ncol(status)) %in% wh]
    status[, whcols] <- as.matrix(data[, dcols])
  }
  time <- as.matrix(time)
  status <- as.matrix(status)
  if (!all(dim(time) == dim(status)))
    stop("unequal dimensions of \"time\" and \"status\" data")
  n <- nrow(time)
  K <- dim(trans)[1]
  if ((dim(trans)[2] != K) | (dim(time)[2] != K))
    stop("dimension of \"trans\" does not match with length of \"time\" and \"status\"")
  idname <- "id"
  if (missing(id))
    id <- 1:n
  else {
    if (!is.vector(id))
      stop("argument \"id\" is not a vector")
    else {
      if (!is.character(id)) {
        if (length(id) != n)
          stop("argument \"id\" not of correct length")
      }
      else {
        if (length(id) == 1) {
          if (n == 1)
            stop("cannot determine whether \"id\" argument indicates ")
          else {
            idname <- id
            id <- data[[id]]
          }
        }
        else {
          if (length(id) != n)
            stop("argument \"id\" not of correct length")
          id <- factor(id)
        }
      }
    }
  }
  idlevels <- NULL
  if (is.factor(id))
    idlevels <- levels(id)
  if (!missing(start)) {
    startstate <- start$state
    starttime <- start$time
    if (length(startstate) != length(starttime))
      stop("starting states and times not of equal length")
    if (length(startstate) > 1) {
      if (length(startstate) != n)
        stop("length of starting states and times different from no of subjects in data")
    }
    else {
      startstate <- rep(startstate, n)
      starttime <- rep(starttime, n)
    }
  }
  else {
    startstate <- rep(1, n)
    starttime <- rep(0, n)
  }
  msres <- msprepEngine(time = time, status = status, id = id,
                        starttime = starttime, startstate = startstate, trans = trans,
                        originalStates = (1:nrow(trans)), longmat = NULL)
  msres <- as.data.frame(msres)
  names(msres) <- c(idname, "from", "to", "trans", "Tstart",
                    "Tstop", "status")
  msres$time <- msres$Tstop - msres$Tstart
  msres <- msres[, c(1:6, 8, 7)]
  msres <- msres[order(msres[, 1], msres[, 5], msres[, 2],
                       msres[, 3]), ]
  row.names(msres) <- 1:nrow(msres)
  if (!is.null(idlevels))
    msres[, 1] <- factor(msres[, 1], 1:length(idlevels),
                         labels = idlevels)
  if (!missing(keep)) {
    if (!(is.matrix(keep) | (is.data.frame(keep)))) {
      if (is.character(keep)) {
        if (missing(data))
          stop("argument \"data\" is missing, with no default")
        nkeep <- length(keep)
        kcols <- match(keep, names(data))
        if (any(is.na(kcols)))
          stop("at least one of elements of \"keep\" not in data")
        keepname <- keep
        keep <- data[, kcols]
      }
      else {
        nkeep <- 1
        keepname <- names(keep)
        if (length(keep) != n)
          stop("argument \"keep\" has incorrect dimension")
      }
    }
    else {
      nkeep <- ncol(keep)
      keepname <- names(keep)
      if (nrow(keep) != n)
        stop("argument \"keep\" has incorrect dimension")
      if (nkeep == 1)
        keep <- keep[, 1]
    }
    if (is.null(keepname))
      keepname <- paste("keep", as.character(1:nkeep),
                        sep = "")
    if (nkeep > 0) {
      if (is.factor(msres[, 1]))
        msres[, 1] <- factor(msres[, 1])
      tbl <- table(msres[, 1])
      if (nkeep == 1) {
        ddcovs <- rep(keep, tbl)
        ddcovs <- as.data.frame(ddcovs)
        names(ddcovs) <- keepname
      }
      else {
        ddcovs <- lapply(1:nkeep, function(i) rep(keep[,
                                                       i], tbl))
        ddcovs <- as.data.frame(ddcovs)
        names(ddcovs) <- keepname
      }
      msres <- cbind(msres, ddcovs)
    }
  }
  attr(msres, "trans") <- trans
  class(msres) <- c("msdata", "data.frame")
  return(msres)
}

msprepEngine <- function (time, status, id, starttime, startstate, trans, originalStates, longmat)
{
  #if (is.null(nrow(time)))
  #  return(longmat)
  #if (nrow(time) == 0)
  #  return(longmat)
  if (length(time) == 0)
    return(longmat)
  states.to <- apply(!is.na(trans), 1, sum)
  absorbing <- which(states.to == 0)
  states.from <- apply(!is.na(trans), 2, sum)
  startings <- which(states.from == 0)
  newstate <- startstate
  newtime <- starttime
  to.remove <- NULL
  for (starting in startings) {
    subjs <- which(startstate == starting)
    nstart <- length(subjs)
    tostates <- which(!is.na(trans[starting, ]))
    transs <- trans[starting, tostates]
    nreach <- length(tostates)
    if ((nstart > 0) & (nreach > 0)) {
      Tstart <- starttime[subjs]
      ###
      if (!is.null(dim(time))) {
        Tstop <- time[subjs, tostates, drop = FALSE]
        Tstop[Tstop < Tstart] <- Inf
        stat <- status[subjs, tostates, drop = FALSE]
        smallesttime <- apply(Tstop, 1, min)
        hlp <- Tstop * 1/stat
        hlp[Tstop == 0 & stat == 0] <- Inf
        nexttime <- apply(hlp, 1, min)
        censored <- which(is.infinite(apply(hlp, 1, min)))
      }
      else {
        Tstop <- time[tostates, drop = FALSE]
        Tstop[Tstop < Tstart] <- Inf
        stat <- status[tostates, drop = FALSE]
        smallesttime <- min(Tstop)
        hlp <- Tstop * 1/stat
        hlp[Tstop == 0 & stat == 0] <- Inf
        nexttime <- min(hlp)
        censored <- which(is.infinite(min(hlp)))
      }
      wh <- which(smallesttime < nexttime)
      whminc <- setdiff(wh, censored)
      if (length(whminc) > 0) {
        whsubjs <- id[subjs[whminc]]
        whsubjs <- paste(whsubjs, collapse = " ")
        warning("From starting state ", originalStates[starting],
                ", subject ", whsubjs, " has smallest transition time with status=0, larger transition time with status=1")
      }
      nexttime[censored] <- smallesttime[censored]
      if (!is.null(dim(time))){
      if (ncol(hlp) > 1) {
        hlpsrt <- t(apply(hlp, 1, sort))
        warn1 <- which(hlpsrt[, 1] - hlpsrt[, 2] == 0)
        if (length(warn1) > 0) {
          isw <- id[subjs[warn1]]
          isw <- paste(isw, collapse = " ")
          hsw <- hlpsrt[warn1, 1]
          hsw <- paste(hsw, collapse = " ")
          warning("Starting from state ", originalStates[starting],
                  ", simultaneous transitions possible for subjects ",
                  isw, " at times ", hsw, "; smallest receiving state chosen")
        }
      }
      }
      if (length(censored) > 0) {
        if (!is.null(dim(hlp))) nextstate <- apply(hlp[-censored, , drop = FALSE], 1, which.min)
        else nextstate <- which.min(hlp[-censored , drop = FALSE])
        reachAbsorb <- (1:nstart)[-censored][which(tostates[nextstate] %in%
          absorbing)]
      }
      else {
        if (!is.null(dim(time))) nextstate <- apply(hlp, 1, which.min)
        else nextstate <- which.min(hlp)
        reachAbsorb <- (1:nstart)[which(tostates[nextstate] %in%
          absorbing)]
      }
      statmat <- matrix(0, nstart, nreach)
      if (length(censored) > 0)
        statmatmin <- statmat[-censored, , drop = FALSE]
      else statmatmin <- statmat
      if (nrow(statmatmin) > 0)
        statmatmin <- t(sapply(1:nrow(statmatmin), function(i) {
          x <- statmatmin[i, ]
          x[nextstate[i]] <- 1
          return(x)
        }))
      if (length(censored) > 0)
        statmat[-censored, ] <- statmatmin
      else statmat <- statmatmin
      mm <- matrix(c(rep(id[subjs], rep(nreach, nstart)),
                     rep(originalStates[starting], nreach * nstart),
                     rep(originalStates[tostates], nstart), rep(transs,
                                                                nstart), rep(Tstart, rep(nreach, nstart)),
                     rep(nexttime, rep(nreach, nstart)), as.vector(t(statmat))),
                   nreach * nstart, 7)
      longmat <- rbind(longmat, mm)
      to.remove <- c(to.remove, subjs[c(censored, reachAbsorb)])
      if (length(censored) > 0)
        newstate[subjs[-censored]] <- tostates[nextstate]
      else newstate[subjs] <- tostates[nextstate]
      if (length(censored) > 0)
        newtime[subjs[-censored]] <- nexttime[-censored]
      else newtime[subjs] <- nexttime
    }
  }
  if (length(to.remove) > 0) {
    if (!is.null(dim(time))){
      time <- time[-to.remove, ]
      status <- status[-to.remove, ]
    }
    else {
      time <- time[-to.remove]
      status <- status[-to.remove]
    }
    newtime <- newtime[-to.remove]
    newstate <- newstate[-to.remove]
    id <- id[-to.remove]
  }
  K <- nrow(trans)
  idx <- rep(1, K)
  idx[startings] <- 0
  idx <- cumsum(idx)
  newstate <- idx[newstate]
  if (!is.null(dim(time))) {
    Recall(time = time[, -startings], status = status[, -startings],
           id = id, starttime = newtime, startstate = newstate,
           trans = trans[-startings, -startings], originalStates = originalStates[-startings],
           longmat = longmat)
  }
  else {
    Recall(time = time[-startings], status = status[ -startings],
           id = id, starttime = newtime, startstate = newstate,
           trans = trans[-startings, -startings], originalStates = originalStates[-startings],
           longmat = longmat)
  }
}
