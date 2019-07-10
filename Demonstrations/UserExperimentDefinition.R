# Input file to construct a GP emulator and perform HM --------------------

experiment_name = "lmdzAWave"
observation_name = "lmdzLES"
# Define the Wave number of iterative refocussing
WAVEN = 1
# Define the value of a threshold
cutoff = 3
# Define a parameter that determines your implausibility measure
tau = 0
# Define the number of points in input space at which you are planning to 
# compute implausibility
sample_size = 10000
# Define your variance of model error
Disc = c(0)
# Define a vector of metrics
metric = c("SCMdata")
