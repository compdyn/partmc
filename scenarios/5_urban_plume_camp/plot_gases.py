from netCDF4 import Dataset
import numpy
import matplotlib
import matplotlib.colors as colors
import matplotlib.patches as pts
import math
import sys
import matplotlib.offsetbox as offsetbox
import cartopy
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

fig, axes = plt.subplots(1,1,
   figsize=(3.5, 3.5))
nt = 60
A = numpy.zeros(nt)
B = numpy.zeros(nt)
C = numpy.zeros(nt)
for ii in range(1,nt+1):
   filename = 'out/urban_plume_0001_%08i.nc' % ii
   f = Dataset(filename, "r", format="NETCDF4")
   gas_mix_rat = f.variables['gas_mixing_ratio'][:]
   A[ii-1] = gas_mix_rat[0]
   B[ii-1] = gas_mix_rat[1]
   C[ii-1] = gas_mix_rat[2]
   f.close()
axes.plot(A, label='A')
axes.plot(B, label='B')
axes.plot(C, label='C')
fig.legend()
fig.savefig("gases.pdf")
