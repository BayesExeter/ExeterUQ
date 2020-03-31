#Need the directory below
#"/Users/danielwilliamson/Dropbox/BayesExeter/mogp_emulator"
twd <- getwd()
library(reticulate)
mogp_emulator <- import("mogp_emulator")
mogp_priors <- import("mogp_emulator.Priors")

#Get some high tune data
setwd("~/Dropbox/BayesExeter/ExeterUQ")
source("Demonstrations/UserExperimentDefinition.R")


load(paste("Demonstrations/", experiment_name, WAVEN, ".RData", sep=""))
load(paste("Demonstrations/", observation_name, ".RData", sep = ""))

dir.create(paste("WAVE", WAVEN, sep=""))
nmetric = length(metric)
cands = names(tData)[-which(names(tData) == metric)]
param = cands[-which(cands == "Noise")]
nparam = length(param)

source('BuildEmulator/BuildEmulator.R')

# GaussianProcess requires data as a matrix, but often you may want to do some
# regression using a data frame in R. To do this, we can split this data frame
# into inputs, targets, and a dictionary mapping column names to integer indices
# using the function below
setwd(twd)

extract_targets <- function(df, target_cols = list("y")) {
  "separate a data frame into inputs, targets, and inputdict for use with GP class"
  
  for (t in target_cols) {
    stopifnot(t %in% names(df))
  }
  
  n_targets <- length(target_cols)
  
  inputs <- matrix(NA, ncol=ncol(df) - n_targets, nrow=nrow(df))
  targets <- matrix(NA, ncol=n_targets, nrow=nrow(df))
  inputdict <- dict()
  
  input_count <- 1
  target_count <- 1
  
  for (n in names(df)) {
    if (n %in% target_cols) {
      targets[,target_count] <- as.matrix(df[n])
    } else {
      inputs[,input_count] <- as.matrix(df[n])
      inputdict[n] <- as.integer(input_count - 1)
      input_count <- input_count + 1
    }
  }
  
  if (n_targets == 1) {
    targets <- c(targets)
  }
  
  return(list(inputs, targets, inputdict))
}

target_list <- extract_targets(tData, list("SCMdata"))
inputs <- target_list[[1]]
targets <- target_list[[2]]
inputdict <- target_list[[3]]

#Fitting Emulators
EMULATOR.LIST <- list()
tau_vec = c()
cutoff_vec = c()

if ( length(EMULATOR.LIST) == WAVEN ) {
  # We shall overwrite the cutoff or tau for the last wave
  tau_vec[WAVEN]=tau
  cutoff_vec[WAVEN]=cutoff
} else {
  listnew <-lapply(1:nmetric,function(k) names(tData[1:nparam]))
  ##############################
  #ExeterUQ version
  #myem.lm <-InitialBasisEmulators(tData=tData,HowManyEmulators=nmetric,additionalVariables = listnew)
  ##############################
  #mo_gp testing code
  tEmulator <- NewEmulators(tData=tData,HowManyEmulators=nmetric,additionalVariables = listnew)
  
  #Functions to call the above
  NewEmulators <- function(tData, HowManyEmulators, additionalVariables=NULL, sigmaPrior = FALSE, 
                           nuggetPrior = FALSE, activePrior = FALSE, activeVariables = NULL, 
                           prior.params = gpstan.default.params, ...){
    lastCand <- which(names(tData)=="Noise")
    tfirst <- lastCand + 1
    if(is.null(HowManyEmulators))
      HowManyEmulators <- length(names(tData)) - lastCand
    lapply(1:HowManyEmulators, function(k) try(BuildMOGP_Emulator(Response=names(tData)[lastCand+k], tData=tData, cands=names(tData)[1:lastCand], 
                                                                  additionalVariables=additionalVariables[[k]], maxdf=ceiling(length(tData[,1])/10)+1, 
                                                                  sigmaPrior = sigmaPrior, nuggetPrior = nuggetPrior, activePrior = activePrior, 
                                                                  activeVariables = activeVariables, prior.params = prior.params, ...), silent = TRUE))
  }
 
  BuildMOGP_Emulator <- function(Response, tData, cands, additionalVariables=NULL, 
                                 canfacs = NULL, TryFouriers = TRUE, maxOrder = 2, maxdf = NULL, 
                                 CompiledModelFit = model_fit, 
                                 sigmaPrior = FALSE, nuggetPrior=FALSE, activePrior = FALSE, 
                                 activeVariables = NULL, prior.params = gpstan.default.params, 
                                 FastVersion = TRUE){
    #' Function to construct Gaussian Proceess emulator in one go
    #' 
    #' @param Response a character corresponding to the response of interest
    #' @param tData a data frame of inputs and a vector of output responses
    #' @param cands a vector of input parameters
    #' @param additionalVariables a vector of parameter names that we want to fit a GP to
    #' that didn't make it into the model
    #' @param canfacs a vector of input parameters that are factors
    #' @param TryFouriers a logical argument with TRUE (default) allowing the Fourier
    #' transformation of input parameters
    #' @param maxOrder maximum order of Fourier terms (Fourier series)
    #' @param maxdf maximum degrees of freedom
    #' @param CompiledModelFit an instance of S4 class stanmodel used for fitting a GP model
    #' @param sigmaPrior a logical argument with FALSE (default) specifying sigma prior parameters
    #' at the values found from the linear model fit
    #' @param nuggetPrior a logical argument with FALSE (default) specifying nugget parameter
    #' at fixed value and TRUE specifying a prior distribution for the nugget parameter.
    #' @param activePrior a logical argument with FALSE (default) specifying the same prior 
    #' for correlation length parameters. TRUE specifying different priors for correlation length
    #' parameters.
    #' @param activeVariables a vector of parameter names considered to be active, with default NULL
    #' @param prior.params default parameters to the prior specification. The default is
    #' gpstan.default.params.
    #' @param FastVersion TRUE value results at saving parameter values at posterior mean, 
    #' FALSE (default) saves the posterior samples for parameters.
    #' 
    #' @return A GP Emulator object
    #' 
    myem.lm <- EMULATE.lm(Response=Response, tData=tData, tcands=cands,
                          tcanfacs=canfacs,TryFouriers=TryFouriers,maxOrder=maxOrder,
                          maxdf = maxdf)
    myem.gp = EMULATE.MOGP(meanResponse=myem.lm, sigmaPrior = sigmaPrior, 
                           nuggetPrior = nuggetPrior, activePrior = activePrior, 
                           activeVariables = activeVariables, tData = tData,
                           additionalVariables=additionalVariables, FastVersion=FastVersion,
                           prior.params = prior.params)
    myem.gp
  }
  
  EMULATE.MOGP <- function(mogp=TRUE,meanResponse, CompiledModelFit = model_fit, sigmaPrior = FALSE, nuggetPrior = FALSE,
                           activePrior = FALSE, activeVariables = NULL, tData, additionalVariables = NULL, FastVersion = FALSE, 
                           prior.params = gpstan.default.params, ...){
    #' Function implement Gaussian process (GP) emulator
    #' 
    #'  @param meanResponse the output of EMULATE.lm containing a linaer model emulator
    #'  @param CompiledModelFit an instance of S4 class stanmodel used for fitting a GP model
    #'  @param sigmaPrior a logical argument with FALSE (default) specifying sigma prior parameters
    #'  at the values found from the linear model fit (NEED MORE CLARIFICATION)
    #'  @param nuggetPrior a logical argument with FALSE (default) specifying nugget parameter at
    #'  fixed value. TRUE specifying a prior distribution for the nugget parameter
    #'  @param activePrior a logical argument with FALSE (default) specifying the same prior 
    #'  for correlation length parameters. TRUE specifying different priors for correlation 
    #'  length parameters.
    #'  @param activeVariables a vector of parametere names considered to be active, with default NULL
    #'  @param tData a data frame containing the inputs and outputs and must be the same as that to create meanResponse 
    #'  @param additionalVariables is a vector of parameter names that we want to fit a GP to that didn't make it into the model
    #'  the default model, with additionalVariables=NULL, is to only treat those terms that make it into the linear model as active
    #'  @param FastVersion TRUE value results at saving parameter values at the posterior mean
    #'  FALSE (default) saves the posterior samples for parameters. 
    #'  @param prior.params a list of parameters to the prior specification. The default is 
    #'  gpstan.default.params. 
    
    #'  @return A GP emulator object
    
    #Design <- as.matrix(tData[,which((names(tData)%in%meanResponse$Names) | (names(tData)%in%meanResponse$Factors) | (names(tData)%in%additionalVariables) | (names(tData) %in% names(meanResponse$Fouriers)))])
    Design <- as.matrix(tData[,which((names(tData)%in%meanResponse$Names) |  (names(tData)%in%additionalVariables) | (names(tData) %in% names(meanResponse$Fouriers)))])
    #Perhaps can take factors out of Design? H still includes them, only effects corr lengths in Stan
    #Slightly careful handling in EMULATOR.gpstan would be required
    tF <- tData[,which(names(tData)==meanResponse$ResponseString)]
    H <- model.matrix(meanResponse$linModel)
    N1 <- dim(Design)[1]
    Np <- dim(H)[2]
    if(mogp){
      tF <- tF - meanResponse$linModel$fitted.values
      theGP <- MultiOutputGP(Design,tF)
      theGP$learn_hyperparameters()
      return(theGP)
    }
    else{
      if(sigmaPrior) prior.params$SwitchSigma = 2
      if(nuggetPrior) {
        prior.params$SwitchNugget = 2
        #prior.params$UpperLimitNugget = prior.params$BetaNugget/(prior.params$AlphaNugget-1) + 10*sqrt(prior.params$BetaNugget^2/((prior.params$AlphaNugget-1)^2*(prior.params$AlphaNugget-2)))
      } 
      #  else {
      #    prior.params$UpperLimitNugget = 2*prior.params$nugget
      #     }
      if(activePrior) prior.params$SwitchDelta = 2
      if(!sigmaPrior) {
        consEm <- EMULATE.lm(Response=meanResponse$ResponseString, tData=tData, tcands="Noise",tcanfacs=NULL,TryFouriers=TRUE,maxOrder=2,maxdf = 0)
        sd2 <- summary(consEm$linModel)$sig - summary(meanResponse$linModel)$sig
        sigsq <- summary(meanResponse$linModel)$sigma
        sigsqvar <- sd2
        # consider the constant Gaussian Process mean
        if(Np == 1) sigsqvar <- summary(meanResponse$linModel)$sig
        prior.params$SigSq <- sigsq
        prior.params$SigSqV <- sigsqvar
      }
      if(activePrior) {
        #if(is.null(activeVariables)) active.inputs <- names(tData)[which((names(tData)%in%meanResponse$Names) | (names(tData)%in%additionalVariables) |(names(tData)%in%names(meanResponse$Fouries)))]
        #else active.inputs <- activeVariables
        #inactive.inputs <- names(tData)[-which((names(tData)%in%active.inputs)|(names(tData)%in%meanResponse$ResponseString))]
        active.inputs <- colnames(Design)[colnames(Design)%in%activeVariables]
        inactive.inputs <- colnames(Design)[-which(colnames(Design)%in%active.inputs)]
        
        Design.active <- as.matrix(tData[, active.inputs])
        Design.inactive <- as.matrix(tData[, inactive.inputs])
        Design <- cbind(Design.active, Design.inactive)
        colnames(Design) <- c(active.inputs, inactive.inputs)
        p.active <- length(active.inputs)
        p.inactive <- length(inactive.inputs)
        p <- dim(Design)[2]
        init.list <- list(list(beta=array(meanResponse$linModel$coefficients, dim = Np), sigma=prior.params$SigSq, nugget = prior.params$nugget, 
                               delta_par=array(c(rep(0.05, p.active), rep(0.7, p.inactive)), dim=p)),
                          list(beta=array(meanResponse$linModel$coefficients, dim = Np), sigma=prior.params$SigSq, nugget = prior.params$nugget, 
                               delta_par=array(c(rep(0.1, p.active), rep(1, p.inactive)), dim=p)))
        StanEmulator <- sampling(CompiledModelFit, data = list(N1 = N1, pact = p.active, pinact = p.inactive, 
                                                               p = p, Np = Np, SwitchDelta = prior.params$SwitchDelta, 
                                                               SwitchNugget = prior.params$SwitchNugget, SwitchSigma = prior.params$SwitchSigma, 
                                                               SigSq = prior.params$SigSq, SigSqV = prior.params$SigSqV, 
                                                               AlphaAct = prior.params$AlphaAct, BetaAct = prior.params$BetaAct, 
                                                               AlphaInact = prior.params$AlphaInact, BetaInact = prior.params$BetaInact, 
                                                               AlphaNugget = prior.params$AlphaNugget, BetaNugget = prior.params$BetaNugget, 
                                                               AlphaRegress = prior.params$AlphaRegress, BetaRegress = prior.params$BetaRegress, 
                                                               nuggetfix = prior.params$nugget, 
                                                               X1 = Design, y1 = tF, H1 = H),
                                 #UpperLimitNugget = prior.params$UpperLimitNugget), 
                                 iter = 2000, warmup = 1000, chains = 2, cores = 2, init = init.list,
                                 pars = c('nugget', 'sigma', 'delta_par', 'beta', 'log_lik'), ...)
      } else {
        # consider the same prior specification for delta_par (correlation length parameter)
        Design <- as.matrix(tData[,which((names(tData)%in%meanResponse$Names) |  (names(tData)%in%additionalVariables) | (names(tData) %in% names(meanResponse$Fouriers)))])
        p <- dim(Design)[2]
        p.active = p.inactive = 1
        init.list <- list(list(beta=array(meanResponse$linModel$coefficients, dim = Np), sigma=prior.params$SigSq, nugget = prior.params$nugget, delta_par=array(rep(0.7, p), dim=p)),
                          list(beta=array(meanResponse$linModel$coefficients, dim = Np), sigma=prior.params$SigSq, nugget = prior.params$nugget, delta_par=array(rep(1, p), dim=p)))
        StanEmulator <- sampling(CompiledModelFit, data = list(N1 = N1, pact = p.active, pinact = p.inactive, 
                                                               p = p, Np = Np, SwitchDelta = prior.params$SwitchDelta, 
                                                               SwitchNugget = prior.params$SwitchNugget, SwitchSigma = prior.params$SwitchSigma, 
                                                               SigSq = prior.params$SigSq, SigSqV = prior.params$SigSqV, 
                                                               AlphaAct = prior.params$AlphaAct, BetaAct = prior.params$BetaAct, 
                                                               AlphaInact = prior.params$AlphaInact, BetaInact = prior.params$BetaInact, 
                                                               AlphaNugget = prior.params$AlphaNugget, BetaNugget = prior.params$BetaNugget,  
                                                               AlphaRegress = prior.params$AlphaRegress, BetaRegress = prior.params$BetaRegress, 
                                                               nuggetfix = prior.params$nugget, #UpperLimitNugget = prior.params$UpperLimitNugget,
                                                               X1 = Design, y1 = tF, H1 = H), 
                                 iter = 2000, warmup = 1000, chains = 2, cores = 2, init = init.list,
                                 pars = c('nugget', 'sigma', 'delta_par', 'beta', 'log_lik'), ...)
      }
      ParameterSamples <- rstan::extract(StanEmulator, pars = c('sigma', 'delta_par', 'beta', 'nugget'))
      if(FastVersion){
        lps <- extract_log_lik(StanEmulator)
        tMAP <- which.max(rowSums(lps))
        A <- ParameterSamples$sigma[tMAP]^2*CovMatrix(Design, ParameterSamples$delta_par[tMAP, ]) #+ diag(ParameterSamples$nugget[tMAP], nrow = dim(Design)[1], ncol = dim(Design)[1])
        QA <- chol(A)
        diff <- tF - H%*%ParameterSamples$beta[tMAP, ]
        Ldiff <- backsolve(QA, diff, transpose=TRUE) #part of the mean update that can be done offline
        FastParts <- list(tMAP=tMAP, A=A, QA=QA, Ldiff=Ldiff)
      }
      else{
        FastParts <- NULL
      }
      gp.list <- list(Design=Design, tF=tF, H=H, ParameterSamples=ParameterSamples, FastParts=FastParts, StanModel = StanEmulator, 
                      prior.params = prior.params, init.list = init.list)
      return(c(meanResponse,gp.list))
    }
  }     
  


