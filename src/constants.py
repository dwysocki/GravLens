"""
Physical constant definitions.

Unit convention is:

  length  : mega parsec
  mass    : solar mass
  time    : second
  charge  : coulomb
  current : ampere
"""
import scipy.constants as cnst
from scipy.constants import kilo, mega

_Mpc = mega*cnst.parsec
_Msun = 1.98855e30

e = cnst.e
c = cnst.c / _Mpc
G = cnst.G / (_Mpc**3 / _Msun)
mₑ = cnst.electron_mass / _Msun
mₚ = cnst.proton_mass / _Msun
ε_0 = cnst.epsilon_0 / (_Msun**-1 * _Mpc**-3)
μ_0 = cnst.mu_0 / (_Msun * _Mpc)
