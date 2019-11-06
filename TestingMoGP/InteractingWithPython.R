library(reticulate)
#Sys.setenv(RETICULATE_PYTHON = "/usr/local/bin/python3")
#use_python("/usr/local/bin/python3")
#reticulate::repl_python()
#  from mogp_emulator.tests import benchmark_branin
#exit
#FROM THE MOGP_EMULATOR DIRECTORY
#DW Run: setwd("/Users/danielwilliamson/Dropbox/BayesExeter/mogp_emulator")
source_python("../ExeterUQ/TestingMoGP/GetPythonFunctions.py") #this worked before making use of branin$ obsolete. Find out why it doesnt work now
branin <- import("mogp_emulator.tests.benchmark_branin")
n_emulators = 3
n_simulations = 15
n_testing = 10
training <- branin$generate_training_data(n_emulators,n_simulations)
names(training) <- c("inputs", "outputs", "em_params")
#Note training$inputs is nxd (15x2) and training$outputs is mxn (3x15)

#Fit the GP, first initialising then estimating parameters
a_gp <- branin$MultiOutputGP(training$inputs, training$outputs)
GPlist <- a_gp$learn_hyperparameters()
#Need to find out how these parameters are output and stored. 
#Seems unclear as no method to return hyperparameters without refitting

Validation <- branin$generate_test_data(n_testing, training$em_params)
names(Validation) <- c("Valid_inputs","Valid_outputs")

#Make predictions
predictions <- a_gp$predict(Validation$Valid_inputs, do_deriv = FALSE, do_unc = TRUE)
names(predictions) <- c("Expectation","Variance","Derivatives")
(Validation$Valid_outputs - predictions$Expectation)/sqrt(predictions$Variance)
(Validation$Valid_outputs - mean(training$outputs))/sqrt(predictions$Variance)
#Note doing much better than a prior, so this must be doing something. 

#Try Save and reload
a_gp$save_emulators(filename = "../ExeterUQ/TestingMoGP/FirstPythonEmulators")
rm(a_gp)
newGP <- branin$MultiOutputGP("../ExeterUQ/TestingMoGP/FirstPythonEmulators.npz")
predictions2 <- newGP$predict(Validation$Valid_inputs, do_deriv = FALSE, do_unc = TRUE)
names(predictions2) <- c("Expectation","Variance","Derivatives")
(Validation$Valid_outputs - predictions2$Expectation)/sqrt(predictions2$Variance)
(Validation$Valid_outputs - mean(training$outputs))/sqrt(predictions2$Variance)

