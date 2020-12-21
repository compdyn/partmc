import pyDOE
import numpy as np

n_scenarios = 1000 # number of scenarios
n_modes = 1 # number of modes

diam_range = [2e-9,50e-6]
sigma_range = [1.4,2.4]
conc_range = [1e5,1e11]

lhs_min = np.zeros((3*n_modes))
lhs_max = np.zeros((3*n_modes))
for ii in range(n_modes):
   lhs_min[ii] = np.log10(diam_range[0])
   lhs_max[ii] = np.log10(diam_range[1])
   lhs_min[ii + n_modes] = sigma_range[0]
   lhs_max[ii + n_modes] = sigma_range[1]
   lhs_min[ii + 2*n_modes] = np.log10(conc_range[0])
   lhs_max[ii + 2*n_modes] = np.log10(conc_range[1])

lhs_prob = pyDOE.lhs(len(lhs_min), n_scenarios)
lhs = lhs_min + (lhs_max-lhs_min) * lhs_prob

# Diameter and number concentrations are in log space
lhs[:,0:n_modes] = 10**lhs[:,0:n_modes]
lhs[:,2*n_modes:] = 10**lhs[:,2*n_modes:]

np.savetxt('lhs.txt', lhs)
