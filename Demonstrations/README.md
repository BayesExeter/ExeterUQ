**Demonstrations** directory contains RData and R files for testing **ExeterUQ** tools. New users do not need to be experts in
R and Stan to use **ExeterUQ** tools. In order to construct an emulator and perform 2 Waves of History Matching, they are only
required to operate with `UserExperimentDefinition.R` and specify the following

`experiment_name`: a string that corresponds to the first part of your RData file that contains
the computer model runs, in particular a data frame `tData`.

`observation_name`: string that corresponds to the first part of your RData file that contains
the obserbations `tObs` and variances of observation errors `tObsErr`.

`cutoff`: corresponds to the threshold value for implausibility measure.

`tau`: corresponds to the parameter that controls the implausibility measure.

`sample_size`: integer value that determines the size of the input space at which we are planning to 
compute implausibility function.

`Disc`: is a vector of variances of discrepancy terms.

`metric`: corresponds to a vector of names of metrics of interest.
