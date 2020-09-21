import pyDOE
import numpy as np

n_scenarios = 1000

diam_range = [2e-8,1e-6]
sigma_range = [1.4,2.4]
conc_range = [1e5,1e11]
n_modes = 3
lhs_min = np.zeros((3*3))
lhs_max = np.zeros((3*3))
for ii in range(3):
   lhs_min[ii] = np.log10(diam_range[0])
   lhs_max[ii] = np.log10(diam_range[1])
   lhs_min[ii + n_modes] = sigma_range[0]
   lhs_max[ii + n_modes] = sigma_range[1]
   lhs_min[ii + 2*n_modes] = np.log10(conc_range[0])
   lhs_max[ii + 2*n_modes] = np.log10(conc_range[1])

print(lhs_min)
print(lhs_max) 
lhs_prob = pyDOE.lhs(len(lhs_min), n_scenarios)
lhs = lhs_min + (lhs_max-lhs_min) * lhs_prob

lhs[:,0:3] = 10**lhs[:,0:3]
lhs[:,6:] = 10**lhs[:,6:]

np.savetxt('lhs.txt', lhs)
