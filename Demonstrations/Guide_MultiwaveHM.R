source('BuildEmulator/BuildEmulator.R')
source("BuildEmulator/DannyDevelopment.R")
source("BuildEmulator/Rotation.R")
source("HistoryMatching/HistoryMatchingFast.R")
library(ncdf4)

# From CONV7.tar
load('WAVE1/EOF_ARMCU_REF_theta/ObsDat.RData') # observations (centred?)
load('WAVE1/EOF_ARMCU_REF_theta/SCMsvd.RData') # wave 1 DataBasis
load('WAVE1/EOF_ARMCU_REF_theta/rotSCM.RData') # wave 1 rotated basis
load('WAVE1/EOF_ARMCU_REF_theta/StanEmsSCM.RData') # wave 1 emulators
load('WAVE1/EOF_ARMCU_REF_theta/tDataSCM.RData') # design + projected wave 1 ensemble
load('WAVE1/EOF_ARMCU_REF_theta/FieldHM.RData') # wave 1 HM results
load('WAVE1/EOF_ARMCU_REF_theta/Disc.RData') # discrepancy

load('WAVE1/EOF_ARMCU_REF_theta/dataLES.RData') # raw observations
# 2nd column of this is the observations we're using here
summary((rotSCM$EnsembleMean + ObsDat) - dataLES[,2]) # zeros

# Assess which wave 1 runs are in the wave 1 NROY space
DiscInv <- GetInverse(Disc)
Wave1Design <- tDataSCM[,1:22]
Wave1DesignHM <- PredictAndHM(rotSCM, ObsDat, StanEmsSCM, tDataSCM, Error = 0*Disc, Disc = Disc, weightinv = DiscInv,
                              Design = Wave1Design)
Wave1DesignNROYinds <- which(Wave1DesignHM$inNROY == TRUE) # selecting those not ruled out

# Load wave 2 runs
######## load from .nc here
#setwd("WAVE2/ARMCU/REF")
#tt <- nc_open("SCM_2-001.nc")
#levels <- ncvar_get(tt, varid = "zf") # changes through time
#times <- ncvar_get(tt, varid = "time")
#ell <- dim(levels)[1] * length(times)
#file_names <- list.files()
#extracted_data_w2 <- matrix(0, nrow = ell, ncol= length(file_names))
#for (i in 1:length(file_names)){
#  tmp_ncdf <- nc_open(file_names[i])
#  tmpdata <-  ncvar_get(tmp_ncdf, varid = "theta")
#  extracted_data_w2[,i] <- c(tmpdata)
#  nc_close(tmp_ncdf)
#}

# extracted_data_w2 is (size of output) x (number of wave 2 runs)
# Add the wave 1 NROY runs to the columns (RAW DATA, NOT CENTRED)
SCMsvdW2 <- MakeDataBasis(data = cbind(extracted_data_w2, extracted_data_w1[,Wave1DesignNROYinds]),
                          weightinv = DiscInv, W = Disc, RemoveMean = TRUE)
# If need to scale:
#SCMsvdW2 <- CentreAndBasis(MeanField = cbind(extracted_data_w2, extracted_data_w1[,Wave1DesignNROYinds]),
#                          weightinv = DiscInv, scaling = 1)


# The wave 2 ensemble has a different mean, hence we also need to re-centre the observations so
# that our comparisons are consistent
ObsDatW2 <- dataLES[,2] - SCMsvdW2$EnsembleMean

# Just in case things have gone badly, check the reconstruction error for each. If has increased, may need more runs
RW_W1 <- ReconError(ObsDat, SCMsvd$tBasis, weightinv = DiscInv, scale = FALSE)
RW_W2 <- ReconError(ObsDat, SCMsvd$tBasis, weightinv = DiscInv, scale = FALSE)

# Rotate the wave 2 basis, if necessary
vSVD_W2 <- VarMSEplot(SCMsvdW2, ObsDatW2, weightinv = DiscInv)
rotSCM_W2 <- RotateBasis(SCMsvdW2, ObsDatW2, kmax = 3, weightinv = DiscInv, v = c(0.45,0.15,0.1), vtot = 0.95, MaxTime = 20)
vROT_W2 <- VarMSEplot(rotSCM_W2, ObsDatW2, weightinv = DiscInv)
qROT_W2 <- which(vROT_W2[,2] > 0.95)[1] # or whichever proportion

# Emulate as usual
setwd('CONV7')
load('WAVE2/Wave2.RData')
Wave2Design <- wave_param_US
# Need to add the wave 1 NROY runs to this
Wave2Design <- rbind(Wave2Design, Wave1Design[Wave1DesignNROYinds,])
# Get data for emulators
tDataSCM_W2 <- GetEmulatableDataWeighted(Design = Wave2Design, EnsembleData = rotSCM_W2, HowManyBasisVectors = qROT_W2, weightinv = DiscInv)
# Only use new runs to fit emulators (use rest for validation etc.)
StanEmsSCM_W2 <- InitialBasisEmulators(tDataSCM_W2[1:dim(wave_param_US)[1],], HowManyEmulators=qROT_W2, TryFourier = TRUE)

# Now history match using the wave 2 emulators, centred obs, and using the results from wave 1 history matching
FieldHM_W2 <- PredictAndHM(rotSCM_W2, ObsDatW2, StanEmsSCM_W2, tDataSCM_W2, Error = 0*Disc, Disc = Disc, weightinv = DiscInv,
                         PreviousWave = FieldHM)


# If we had multiple metrics at the previous wave, now it's easier to run at the same design
# For metric 1, do as usual
#FieldHM
# For metric k = 2,...m, do:
#FieldHMk <- PredictAndHM(, ..., Design = FieldHM$Design) # evaluates at same design
# Then if we want to rule out runs based on some multi metric, for example here let's assume that we want all of them to be below threshold
# As we've evaluated all at same design, can find NROY runs as follows (for 3 metrics here):
#inNROY_mm <- which(FieldHM$inNROY == TRUE & FieldHM2$inNROY == TRUE & FieldHM3$inNROY == TRUE)
# For optimal use at wave 2, edit the $inNROY vector in each to indicate whether in the multi-metric NROY space
#FieldHM$inNROY[inNROY_mm] <- FieldHM2$inNROY[inNROY_mm] <- FieldHM3$inNROY[inNROY_mm] <- TRUE
#FieldHM$inNROY[-inNROY_mm] <- FieldHM2$inNROY[-inNROY_mm] <- FieldHM3$inNROY[-inNROY_mm] <- FALSE
# so that when do wave 2, only evaluate at the inNROY_mm design runs



