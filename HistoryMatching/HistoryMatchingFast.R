#' History matching
#'
#' Calculates the implausibility (for the \ell-dimensional field) for a sample of expectations and variances from basis emulators
#'
#' @param DataBasis object containing the basis used in emulation ($tBasis)
#' @param Obs observation vector (length \ell), must be centred
#' @param Expectation a matrix containing emulator expectations, where a given row contains the expectations for the q emulated basis vectors, for some x
#' @param Variance a matrix containing emulatorvariances, where a given row contains the variances for the q emulated basis vectors, for some x
#' @param Error observation error variance matrix
#' @param Disc discrepancy variance matrix
#' @param weightinv if not NULL, the inverse of W = var_err + var_disc, used for projection
#' #' 
#' @return \item{impl}{Vector of implausibilities corresponding to the rows of Expectation and Variance}
#' \item{bound}{The chi-squared bound for an \ell-dimensional field}
#' \item{nroy}{Proportion of parameter settings that are not ruled out, using bound}
#' \item{inNROY}{Vector indicating whether a parameter setting is ruled out}
#'
#' @export
HistoryMatch <- function(DataBasis, Obs, Expectation, Variance, Error, Disc, weightinv = NULL, BasisUncertainty = FALSE){
  q <- dim(Expectation)[2]
  Basis <- DataBasis$tBasis[,1:q]
  l <- dim(Basis)[1]
  stopifnot(q == dim(Variance)[2])
  W <- Error + Disc
  if (is.null(weightinv)){
    weightinv <- GetInverse(W)
  }
  R_W <- ReconError(Obs, Basis, weightinv = weightinv, scale = FALSE)
  # Project observations onto basis if required
  if (length(Obs) == l){
    ObsProj <- CalcScores(Obs, Basis, weightinv = weightinv)
  }
  # Add uncertainty from discarded basis vectors?
  if (BasisUncertainty == TRUE){
    BasisVar <- DiscardedBasisVariance(DataBasis, q, weightinv)
    W <- W + BasisVar
  }
  # Project variance matrices onto basis if required
  if (dim(Disc)[1] == l){
    WProj <- VarProj(W, Basis, weightinv = weightinv)
  }
  nn <- dim(Expectation)[1]
  impl <- as.numeric(mclapply(1:nn, function(i) ImplCoeff(Expectation[i,], Variance[i,], ObsProj, WProj, 0*WProj)))
  impl <- impl + rep(R_W, nn) 
  bound <- qchisq(0.995, l)
  nroy <- sum(impl < bound)/nn
  inNROY <- impl < bound
  return(list(impl = impl, bound = bound, nroy = nroy, inNROY = inNROY))
}

#' Coefficient implausibility
#'
#' Calculates the coefficient implausibility for a single x, given projected quantities
#'
#' @param Expectation length q vector with emulator expectations
#' @param Variance length q vector with emulator variances
#' @param Obs projected observations
#' @param Error projected observation error variance matrix
#' @param Disc projected discrepancy variance matrix
#'  
#' @return The coefficient implausibility (given the matrix used in projection)
#'
#' @export
ImplCoeff <- function(Expectation, Variance, Obs, Error, Disc){
  V <- Error + Disc + diag(Variance)
  Q <- chol(V)
  proj.output <- Expectation
  y <- backsolve(Q, as.vector(Obs - proj.output), transpose = TRUE)
  impl <- crossprod(y,y)
  return(impl)
}

#' Setting discrepancy multiple
#'
#' Scaling the discrepancy to ensure that the observations won't be ruled out
#'
#' @param Basis full basis
#' @param q where the basis is truncated
#' @param obs observations
#' @param level quantile of the chi-squared distribution to use (< 0.995)
#' @param weightinv inverse of W, to use in projection
#'  
#' @return scalar to be used as discrepancy multiplier, to ensure observations not ruled out
#'
#' @export
SetDiscrepancy <- function(tBasis, q, obs, level = 0.95, weightinv = NULL){
  TruncatedError <- ReconError(obs, tBasis[,1:q], weightinv = weightinv, scale = FALSE)
  l <- dim(tBasis)[1]
  b <- qchisq(level, l)
  DiscMultiplier <- c(TruncatedError / b)
  return(DiscMultiplier)
}


#' Prediction and History Matching
#'
#' Takes emulators, evaluates expectations and varianes for space-filling design, and history matches
#'
#' @param DataBasis object containing the basis used in emulation ($tBasis)
#' @param Obs observation vector (length \ell), must be centred
#' @param Ems a Stan emulator object
#' @param tData matrix containing parameter values
#' @param ns number of parameter settings to evaluate emulators at
#' @param Error observation error variance matrix
#' @param Disc discrepancy variance matrix
#' @param weightinv if not NULL, the inverse of W = var_err + var_disc, used for projection
#' @param Design if not NULL, passes a design at which to evaluate emulators and implausibility
#' @param PreviousWave if not NULL, provides the output of a previous PredictAndHM object, and evaluates the current NROY points from the previous design
#' #' 
#' @return \impl{Design}{Space-filling design of ns points at which the emulators were evaluated}
#' \item{Expectation}{Emulator expectations}
#' \item{Variance}{Emulator variances}
#' \item{impl}{Vector of implausibilities corresponding to the rows of Expectation and Variance}
#' \item{bound}{The chi-squared bound for an \ell-dimensional field}
#' \item{nroy}{Percentage of parameter settings that are not ruled out, using bound}
#' \item{inNROY}{Vector indicating whether a parameter setting is ruled out}
#'
#' @export
PredictAndHM <- function(DataBasis, Obs, Ems, tData, ns = 1000, Error, Disc, weightinv = NULL, BasisUncertainty = FALSE,
                         Design = NULL, PreviousWave = NULL){
  T_f <- qchisq(0.995, dim(DataBasis$tBasis)[1]) # bound for ruling out on field
  npar <- which(colnames(tData) == "Noise") - 1
  
  design_flag <- !(is.null(Design))
  
  if (is.null(Design) & is.null(PreviousWave)){
    Design <- 2*as.data.frame(randomLHS(ns, npar)) - 1
    colnames(Design) <- colnames(tData)[1:npar]
  }
  
  if (!(is.null(PreviousWave))){
    inNROY_inds <- which(PreviousWave$inNROY == TRUE)
    Design <- PreviousWave$Design[inNROY_inds,] # only evaluate at not ruled out points
  }
  
  EmOutput <- lapply(1:length(Ems), function(e) EMULATOR.gpstan(Design,Ems[[e]], FastVersion = TRUE))
  Expectation <- Variance <- matrix(0, nrow = dim(Design)[1], ncol = length(Ems))
  for (j in 1:length(Ems)){
    Expectation[,j] <- EmOutput[[j]]$Expectation
    Variance[,j] <- EmOutput[[j]]$Variance
  }
  
  FieldHM <- HistoryMatch(DataBasis, Obs, Expectation, Variance, Error, Disc, weightinv = weightinv)

  if (!(is.null(PreviousWave))){
    FieldHM$nroy <- PreviousWave$nroy * FieldHM$nroy
  }
  
  if (design_flag == TRUE){
    print("Proportion of given design not ruled out:")
  }
  else {
    print("Proportion of original space not ruled out:")
  }
  print(FieldHM$nroy)
  return(list(Design = Design, Expectation = Expectation, Variance = Variance, impl = FieldHM$impl, bound = FieldHM$bound, nroy = FieldHM$nroy, inNROY = FieldHM$inNROY))
}

