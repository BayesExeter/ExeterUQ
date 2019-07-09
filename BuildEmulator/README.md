# README

Please make sure that the following packages are installed on RStudio:

* GenSA
* far
* fields
* lhs
* maps
* mco
* mvtnorm
* ncdf4
* parallel
* rstan
* shape
* tensor

In order to install and load 'lhs' package in RStudio console type:
```
install.packages("lhs")
library("lhs")
```
In order to use the functions from `BayesCalibration`, `Diagnostics`, `SensitivityAnalysis`, `UncertaintyAnalysis` folders, please make sure that your emulator, e.g. `EMULATE.gpstan(...)`, returns a list containing the following components:

* `Names`: a named vector of coefficients corresponding to the linear terms in the model.
* `linModel`: an object of class "formula": a symbolic description of the model to be fitted, together with the coefficients' values.
* `mainEffects`:
* `Interactions`: any interaction terms found.
* `Factors`: any factor terms specified.
* `FactorInteractions`: any factor interactions specified.
* `ThreeWayInters`: any three way interactions specified.
* `Fouriers`:
* `pre.Lists`:
* `DataString`: the naming of the data frame.
* `ResponseString`: the naming of the response variable.
* `Design`: the design matrix.
* `tF`: the vector of response.
* `H`: the regression matrix.
* `ParameterSamples`: the list of posterior samples for GP Emulator parameters `sigma_sq`, `delta_par` and `beta`.
* `StanPredict`:
* `FastParts`:
* `nugget`: the pre-defined nugget term value.




