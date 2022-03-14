#Read data
#Test
import numpy as np
from matplotlib import ticker, cm
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import seaborn as sns
mpl.rcParams['font.size'] = 12

#Mixing experiment
data  = pd.read_csv('out_pfr_suc2as_1/exp_suc2as_1.txt',sep='\s+|,', skiprows=19, nrows=102, header=None, engine = 'python')
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


#Premixing AS
data_as  = pd.read_csv('out_pfr_suc2as_1/exp_suc2as_aspre.txt',sep='\s+|,', skiprows=19, nrows=102, header=None, engine = 'python')
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
data_suc  = pd.read_csv('out_pfr_suc2as_1/exp_suc2as_sucpre.txt',sep='\s+|,', skiprows=19, nrows=102, header=None, engine = 'python')
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

    
    
print(dist.shape)    
plt.figure(figsize=(10,7))
plt.xscale('log')
plt.plot(diam_as, dist_as[1,:], label="Time:" + new_data_as.iloc[1+2,0])
plt.plot(diam_as, dist_as[3,:], label="Time:" + new_data_as.iloc[3+2,0])

plt.plot(diam_suc, dist_suc[1,:], label="Time:" + new_data_suc.iloc[1+2,0])
plt.plot(diam_suc, dist_suc[3,:], label="Time:" + new_data_suc.iloc[3+2,0])


plt.plot(diam, dist[10,:], label="Time:" + new_data.iloc[10+2,0])
plt.plot(diam, dist[34,:], label="Time:" + new_data.iloc[34+2,0])
plt.plot(diam, dist[35,:], label="Time:" + new_data.iloc[35+2,0])
plt.plot(diam, dist[36,:], label="Time:" + new_data.iloc[36+2,0])
plt.plot(diam, dist[37,:], label="Time:" + new_data.iloc[37+2,0])
plt.plot(diam, dist[38,:], label="Time:" + new_data.iloc[38+2,0])
plt.plot(diam, dist[39,:], label="Time:" + new_data.iloc[39+2,0])
plt.plot(diam, dist[40,:], label="Time:" + new_data.iloc[40+2,0])
plt.plot(diam, dist[41,:], label="Time:" + new_data.iloc[41+2,0])
plt.plot(diam, dist[42,:], label="Time:" + new_data.iloc[42+2,0])

plt.xlabel('Diameter, nm')
plt.ylabel('Number concentration')
plt.legend(loc=0, ncol=1)
#for i in range(20):
#    plt.xscale('log')
#    plt.plot(diam, dist[i,:], label="Time:" + new_data.iloc[i+2,0])
#    plt.xlabel('Diameter, nm')
#    plt.ylabel('Number concentration')
#    plt.legend(loc=1, ncol=2)
plt.savefig("out_pfr_suc2as_1/size_dist_exp.pdf")    
