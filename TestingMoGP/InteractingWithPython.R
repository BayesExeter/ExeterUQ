library(reticulate)
use_python("/usr/local/bin/python3")
Sys.setenv(RETICULATE_PYTHON = "/usr/local/bin/python3")
reticulate::repl_python()
  from mogp_emulator.tests import benchmark_branin
exit

