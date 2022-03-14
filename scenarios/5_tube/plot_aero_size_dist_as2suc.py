#!/usr/bin/env python

import sys, os
sys.path.append("../../tool")
import mpl_helper
import matplotlib
import scipy.io
import numpy as np
import pandas as pd

#Mixing experiment
data  = pd.read_csv('exp_data/exp_as2suc_3.txt',sep='\s+|,', skiprows=19, nrows=102, header=None, engine = 'python')
new_data = data.T
new_data.columns = new_data.iloc[0]

diam = np.zeros(100)
dist = np.zeros((new_data.shape[0]-2,100))

#Read the diameter
for i in range(100):
    j = i+2
    diam[i] = new_data.iloc[0,j]
#Read the distribution
for i in range(new_data.shape[0]-2):
    j = i +2
    dist[i,:] = new_data.iloc[j,2:]

print(dist.shape)    

#Premixing AS
data_as  = pd.read_csv('exp_data/exp_as2suc_aspre.txt',sep='\s+|,', skiprows=19, nrows=102, header=None, engine = 'python')
new_data_as = data_as.T
new_data_as.columns = new_data_as.iloc[0]

diam_as = np.zeros(100)
dist_as = np.zeros((new_data.shape[0]-2,100))

#Read the diameter
for i in range(100):
    j = i+2
    diam_as[i] = new_data_as.iloc[0,j]
    
#Read the distribution
for i in range(new_data_as.shape[0]-2):
    j = i + 2
    dist_as[i,:] = new_data_as.iloc[j,2:]


#Premixing Suc
data_suc  = pd.read_csv('exp_data/exp_as2suc_sucpre.txt',sep='\s+|,', skiprows=19, nrows=102, header=None, engine = 'python')
new_data_suc = data_suc.T
new_data_suc.columns = new_data_suc.iloc[0]

diam_suc = np.zeros(100)
dist_suc = np.zeros((new_data_suc.shape[0]-2,100))

#Read the diameter
for i in range(100):
    j = i+2
    diam_suc[i] = new_data_suc.iloc[0,j]
    
#Read the distribution
for i in range(new_data_suc.shape[0]-2):
    j = i + 2
    dist_suc[i,:] = new_data_suc.iloc[j,2:]    

# Model output
num_dist_plot = np.zeros([180,9])
i = 0
for (filename, index) in mpl_helper.get_filename_list('out_pfr_as2suc_3/', r'urban_plume_([0-9]+)_process\.nc'):
    ncf = scipy.io.netcdf_file(filename)
    diam_pmc = ncf.variables["diam"].data.copy() * 1e9
    num_dist = ncf.variables["num_dist"].data.copy() * np.log(10) / 1e6
    print(num_dist.shape, diam.shape)
    num_dist_plot[:,i] = num_dist
    i = i + 1


(figure, axes) = mpl_helper.make_fig(left_margin=1, right_margin=0.5)

#axes.plot(diam_as, dist_as[1,:], label="AS premix")
#axes.plot(diam_as, dist_as[3,:], label="AS premix")

axes.plot(diam_suc, dist_suc[1,:], label="Sucrose premix")
axes.plot(diam_suc, dist_suc[3,:], label="Sucrose premix")

### 1 --- 34
### 2 --- 40
### 3 --- 38
axes.plot(diam, dist[35,:], label="Exp Time:" + new_data.iloc[35+2,0])
axes.plot(diam, dist[36,:], label="Exp Time:" + new_data.iloc[36+2,0])
axes.plot(diam, dist[37,:], label="Exp Time:" + new_data.iloc[37+2,0])
#axes.plot(diam, dist[38,:], label="Exp Time:" + new_data.iloc[38+2,0])
#axes.plot(diam, dist[39,:], label="Exp Time:" + new_data.iloc[39+2,0])

axes.plot(diam_pmc, num_dist_plot[:,0], label="pmc 0 min")
#axes.plot(diam_pmc, num_dist_plot[:,1], label="10 min")
#axes.plot(diam_pmc, num_dist_plot[:,2], label="20 min")
#axes.plot(diam_pmc, num_dist_plot[:,3], label="30 min")
#axes.plot(diam_pmc, num_dist_plot[:,4], label="40 min")
axes.plot(diam_pmc, num_dist_plot[:,5], label="pmc 50 min")
axes.plot(diam_pmc, num_dist_plot[:,6], label="pmc 60 min")
#axes.plot(diam_pmc, num_dist_plot[:,7], label="70 min")
#axes.plot(diam_pmc, num_dist_plot[:,8], label="80 min")

    
axes.set_xscale("log")
axes.set_xlabel(r"dry diameter $D_{\rm dry}$ / nm")
axes.set_xlim(10, 1e3)

axes.set_yscale("linear")
axes.set_ylabel(r"dN/dlogDp / $\rm cm^{-3}$")

axes.grid(True)
axes.legend()
out_filename = "out_pfr_as2suc_3/urban_plume_num_dist.pdf"
figure.savefig(out_filename)
matplotlib.pyplot.close(figure)
