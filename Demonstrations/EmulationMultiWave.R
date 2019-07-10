# Wrapper to perform History Matching - so this is the inputs that the wrapper should ask?
##############################
twd <- getwd()
source("Demonstrations/UserExperimentDefinition.R")

load(paste("Demonstrations/", experiment_name, WAVEN, ".RData", sep=""))
load(paste("Demonstrations/", observation_name, ".RData", sep = ""))

dir.create(paste("WAVE", WAVEN, sep=""))
nmetric = length(metric)
cands = names(tData)[-which(names(tData) == metric)]
param = cands[-which(cands == "Noise")]
nparam = length(param)


# Produce scatter plots of model output -----------------------------------

file = paste("WAVE",WAVEN,"/Plots_Metrics.pdf",sep="")
pdf(file=file)

if(nparam*nmetric<2){ par(mfrow = c(1, 1), mar=c(4, 4, 1, 1)) 
} else if(nparam*nmetric <= 4){ par(mfrow = c(2, 2), mar=c(4, 4, 1, 1)) 
} else if(nparam*nmetric <= 9){ par(mfrow = c(3, 3), mar=c(4, 4, 1, 1)) 
} else if(nparam*nmetric <= 16){ par(mfrow = c(4, 4), mar=c(4, 4, 1, 1))
} else if(nparam*nmetric <= 25){ par(mfrow = c(5,5), mar=c(4, 4, 1, 1)) }
for ( iparam in 1:nparam ) {
  for (j in 1:nmetric) {
    ymin=tObs[j]-sqrt(tObsErr[j])
    ymax=tObs[j]+sqrt(tObsErr[j])
    plot(tData[,iparam],tData[,nparam+j+1],col=2,xlab=cands[iparam],ylab=metric[j],
         ylim=range(tData[,nparam+j+1],ymin,ymax))
    abline(h=c(tObs[j],ymin,ymax ), col = 'blue', lty =2, lwd=c(1,3,3))
  }
}
dev.off()


# Start constructing GP emulators or load emulator for previous wave --------

source('BuildEmulator/BuildEmulator.R')
EMULATOR.LIST <- list()
tau_vec = c()
cutoff_vec = c()

if(file.exists("EMULATOR_LIST_MULT_METRIC.RData")) { load("EMULATOR_LIST_MULT_METRIC.RData") }

if ( length(EMULATOR.LIST) == WAVEN ) {
  # We shall overwrite the cutoff or tau for the last wave
  tau_vec[WAVEN]=tau
  cutoff_vec[WAVEN]=cutoff
} else {
  listnew <-lapply(1:nmetric,function(k) names(tData[1:nparam]))
  myem.lm <-InitialBasisEmulators(tData=tData,HowManyEmulators=nmetric,additionalVariables = listnew)
  for(i in 1:nmetric) myem.lm[[i]]$StanModel = NULL
  file = paste("WAVE",WAVEN,"/Plots_LOO.pdf",sep="")
  pdf(file=file)
  #DIAGNOSTICS : verifying that the original SCMdata are well reproduced when eliminating the point from the emulator
  for (i in 1:nmetric) {
    tLOOs1 <- LOO.plot(StanEmulator = myem.lm[[i]], ParamNames=names(tData)[1:nparam])
  }
  dev.off()
  
  EMULATOR.LIST[[WAVEN]] = myem.lm
  tau_vec[WAVEN] = tau
  #SAVE EMULATOR FOR HISTORY MATCHING
  cutoff_vec[WAVEN] = cutoff
  save(EMULATOR.LIST, tau_vec, cutoff_vec, file = "EMULATOR_LIST_MULT_METRIC.RData")
  print(paste("A list of emulators has been saved under: ","EMULATOR_LIST_MULT_METRIC.RData",sep=""))
}


# Proceed to performing History Matching or iterative refocussing ---------

source('HistoryMatching/HistoryMatching.R')
Xp <- as.data.frame(2*randomLHS(sample_size, nparam)-1)
names(Xp) <- names(tData)[1:nparam]

if(WAVEN == 1) {
  Timps <- ManyImplausibilitiesStan(NewData=Xp, Emulator=EMULATOR.LIST[[1]], Discrepancy=Disc,
                                    Obs=tObs, ObsErr=tObsErr, is.GP=NULL,FastVersion = TRUE)
                                    #multicore=(ceiling(dim(Xp)[1]/nbatches)>1), batches=nbatches) # Multicore only if Xp is large enough
  ImpData_wave1 = cbind(Xp, Timps)
  VarNames <- names(Xp)
  valmax = tau_vec[1] + 1
  ImpListM1 = CreateImpList(whichVars = 1:nparam, VarNames=VarNames, ImpData=ImpData_wave1, nEms=length(EMULATOR.LIST[[1]]), whichMax=valmax)
  file = paste("WAVE", WAVEN,"/InputSpace_wave",WAVEN,".pdf",sep="")
  pdf(file=file)
  imp.layoutm11(ImpListM1,VarNames,VariableDensity=FALSE,newPDF=FALSE,the.title=paste("InputSpace_wave",WAVEN,".pdf",sep=""),newPNG=FALSE,newJPEG=FALSE,newEPS=FALSE)
  dev.off()
  NROY1 <- which(rowSums(Timps <= cutoff_vec[1]) >=length(EMULATOR.LIST[[1]]) - tau_vec[1])
  TMimpls_wave1 <- apply(Timps, 1,MaxImp,whichMax=valmax)
  XpNext <- Xp[NROY1, ]
  # number of plausible members divided by total size -> fraction of space retained after wave N
  print(paste("Remaining space after wave 1: ",length(NROY1)/dim(Xp)[1],sep=""))
} else {
  XpNext = Xp
  NROY.list = list()
  Impl.list = list()
  for(i in 1:(length(EMULATOR.LIST))) {
    Timps = ManyImplausibilitiesStan(NewData=XpNext, Emulator=EMULATOR.LIST[[i]], Discrepancy=Disc,
                                     Obs=tObs, ObsErr=tObsErr, is.GP=NULL,FastVersion = TRUE)
                                     #multicore=(ceiling(dim(XpNext)[1]/nbatches)>1), batches=nbatches) # Multicore only if XpNext is large enough
    valmax = tau_vec[i] + 1
    print(c("cutoff",i,cutoff_vec[i]))
    Impl.list[[i]] = matrix(apply(Timps, 1,MaxImp,whichMax=valmax), ncol = 1)
    NROY.list[[i]] = which(rowSums(Timps <= cutoff_vec[i]) >= length(EMULATOR.LIST[[i]]) - tau_vec[i])
    XpNext = XpNext[NROY.list[[i]], ]
    print(paste("Remaining space after wave",i,": ",length(NROY.list[[i]])/dim(Xp)[1],sep=""))
  }
  Impl.list[[length(EMULATOR.LIST)]] = Timps
  ImpData <- ImpDataWaveM(Xp, NROY.list, Impl.list)
  VarNames <- names(Xp)
  ImpList <- CreateImpListWaveM(whichVars = 1:nparam, VarNames=VarNames, ImpData = ImpData,
                                nEms=length(EMULATOR.LIST[[length(EMULATOR.LIST)]]), Resolution=c(15,15), whichMax=valmax)
  file = paste("WAVE", WAVEN,"/InputSpace_wave",WAVEN,".pdf",sep="")
  pdf(file=file)
  imp.layoutm11(ImpList,VarNames,VariableDensity=FALSE,newPDF=FALSE,the.title=paste("InputSpace_wave",WAVEN,".pdf",sep=""),newPNG=FALSE,newJPEG=FALSE,newEPS=FALSE,Points=NULL)
  dev.off()
}

print(paste("Created figure file: ", "InputSpace_wave",WAVEN,".pdf", sep=""))









