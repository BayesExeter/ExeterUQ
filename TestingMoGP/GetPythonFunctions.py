
import numpy as np
from mogp_emulator import MultiOutputGP
from mogp_emulator import MonteCarloDesign, LatinHypercubeDesign
from scipy.stats import uniform
try:
    import matplotlib.pyplot as plt
    makeplots = True
except ImportError:
    makeplots = False
    
from mogp_emulator.tests import benchmark_branin
def f():
  return benchmark_branin.__file__
  



