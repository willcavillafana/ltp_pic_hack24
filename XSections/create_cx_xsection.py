import math
import numpy as np
import sys


a = 0.0543
Eref = 1.0
Sref = 55.3e-20


E = np.logspace(-3,4,141)
SIGMA = Sref*(1.0 + a*np.log(np.divide(Eref,E)))**2

OUT = np.zeros((len(E),2))

OUT[:,0] = E
OUT[:,1] = SIGMA

np.savetxt('ar+_ar_cxchange.dat',OUT)
