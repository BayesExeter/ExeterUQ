# ExeterUQ

**ExeterUQ** is an open repository developed by researchers at the University of Exeter for Uncertainty Quantification (UQ) primarily in R and RStan (the Stan interface for R). **ExeterUQ** possesses a number of key and unique features such as 

1. Bayesian GP emulation for a single model output or multivariate outputs (time-series, spatial fields)
2. Calibration via history matching
3. Multi-wave history matching (refocussing)

**Directories**

`BuildEmulator` cointains all the necessary R and Stan files to construct a GP emulator, perform diagnostics (Leave One Out)
and predictions.

`HistoryMatching` contains all the necassary R files to perform and visualize history matching and iterative refocussing.

`Demonstartions` contains all the R files and RData files to perform a simple demonstrations of the functionality of **ExeterUQ** tools.

## Introduction to Gaussian Process Emulators
