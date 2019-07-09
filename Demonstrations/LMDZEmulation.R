
# Load all the necessary source files -------------------------------------
twd <- getwd()
source('BuildEmulator/BuildEmulator.R')
source('HistoryMatching/HistoryMatching.R')

# Load the Ensemble 1 LMDZ ------------------------------------------------
load(file="Demonstrations/lmdzAWave1.RData")
head(tData)

# Construct a GP emulator (mean function and full GP) ---------------------
cands <- names(tData)[-which(names(tData)=="SCMdata")]
myem.lm.lmdza1 <- EMULATE.lm(Response="SCMdata", tData=tData, tcands=cands,tcanfacs=NULL,TryFouriers=TRUE,maxOrder=2,maxdf = 20)
myem.gp.lmdza1 = EMULATE.gpstan(meanResponse=myem.lm.lmdza1, tData=tData, additionalVariables=NULL,FastVersion = TRUE)

# Leave-One-Out Diagnostics for GP emulator -------------------------------
tLOOs <- LOO.plot(StanEmulator = myem.gp.lmdza1, ParamNames = myem.gp.lmdza1$Names)

# Define data necessary for History Matching ------------------------------
load(file="Demonstrations/LESens.RData")
tObs <- LESens[1,48]
tObsErr <- sd(LESens[,48])^2
Disc <- 0
Xpred <- 2*randomLHS(10000,6) - 1
Xp <- as.data.frame(Xpred)


# Perform and visualize a single wave of History Matching -----------------
names(Xp) <- names(tData)[1:6]
Timps1 <- UniImplausibilityStan(NewData=Xp, Emulator=myem.gp.lmdza1, Discrepancy=Disc, Obs=tObs, ObsErr=tObsErr, is.GP=NULL,FastVersion = TRUE)
ImpData1 <- cbind(Xp,Timps1)
ImpList <- CreateImpList(whichVars = 1:6, VarNames = VarNames, ImpData = ImpData1,Resolution = c(15,15), whichMax=1)
imp.layoutm11(ImpList,cands[1:6],VariableDensity=FALSE,newPDF=FALSE,the.title="InputSpace.png",newPNG=FALSE,newJPEG=FALSE,newEPS=FALSE,Points=NULL)
print(paste("The percentage of Wave 1 NROY space", length(which(Timps1 < 3))/dim(Xp)[1]*100, "%", sep = " "))
NROY1 <- which(Timps1 < 3)
XpNext <- Xp[NROY1, ]

# Load the Ensemble 2 LMDZ ------------------------------------------------
load(file="Demonstrations/lmdzAWave2.RData")
head(tData)


# Construct a GP emulator -------------------------------------------------
cands <- names(tData)[-which(names(tData)=="SCMdata")]
myem.lm.lmdza2 <- EMULATE.lm(Response="SCMdata", tData=tData, tcands=cands,tcanfacs=NULL,TryFouriers=TRUE,maxOrder=2,maxdf = 20)
myem.gp.lmdza2 = EMULATE.gpstan(meanResponse=myem.lm.lmdza2, tData=tData, additionalVariables=NULL, FastVersion = TRUE)


# Leave-One-Out Diagnostics for GP emulator -------------------------------
tLOOs <- LOO.plot(StanEmulator = myem.gp.lmdza2, ParamNames = myem.gp.lmdza2$Names)

# Compute implausibility at the NROY1 input space -----------------
Timps2 <- UniImplausibilityStan(NewData=XpNext, Emulator=myem.gp.lmdza2, Discrepancy=Disc, Obs=tObs, ObsErr=tObsErr, is.GP=NULL,FastVersion = TRUE)
NROY2 <- which(Timps2 < 3)

# Perform and visualize a second wave of History Matching -----------------
NROY.list = list()
Impl.list = list()
NROY.list[[1]] = matrix(NROY1, ncol = 1)
NROY.list[[2]] = matrix(NROY2, ncol = 1)

Impl.list[[1]] = matrix(Timps1, ncol = 1)
Impl.list[[2]] = matrix(Timps2, ncol = 1)

ImpData <- ImpDataWaveM(Xp, NROY.list, Impl.list)

VarNames <- names(Xp)
ImpList <- CreateImpListWaveM(whichVars = 1:6, VarNames=VarNames, 
                              ImpData = ImpData,
                              nEms= 1, Resolution=c(15,15), whichMax=1)
imp.layoutm11(ImpList,VarNames,VariableDensity=FALSE,newPDF=FALSE,
              the.title=paste("InputSpace_wave",WAVEN,".pdf",sep=""),newPNG=FALSE,newJPEG=FALSE,newEPS=FALSE,
              Points= NULL)