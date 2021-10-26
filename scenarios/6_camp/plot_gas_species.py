from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import matplotlib.offsetbox as offsetbox
import matplotlib.font_manager
from matplotlib import ticker


colors = ['#66c2a5','#fc8d62','#8da0cb']

def get_gas_mixing_ratio(fileprefix, n_times, species, n_repeats):
    gas_mix_rat = np.zeros(n_times)
    for i_repeat in range(n_repeats):
        for i_time in range(n_times):
            file = '%s_%04i_%08i.nc' % (fileprefix,i_repeat+1,i_time+1)
            f = Dataset(file, "r", format="NETCDF4")
            gas_names = f.variables['gas_species'].names.split(",")
            gas_values = f.variables['gas_mixing_ratio'][:]
            gas_mix_rat[i_time] += gas_values[gas_names.index(species)]      
    return gas_mix_rat / n_repeats

def get_times(fileprefix, n_times):
    times = np.zeros(n_times)
    for i_time in range(n_times):
        file = '%s_0001_%08i.nc' % (fileprefix,i_time+1)
        f = Dataset(file, "r", format="NETCDF4")
        times[i_time] = (f.variables['time'][:]) / 60
        f.close()
    return times

def get_gas_species(fileprefix):
    i_time = 1
    file = '%s%08i.nc' % (fileprefix,i_time)
    f = Dataset(file, "r", format="NETCDF4")
    
    return(f.variables['gas_species'].names.split(","))


part_label = 'CAMP-part'
part_symbol = 'o'

n_times = 145
n_runs = 1

page_width = 16/2.54
page_height = 2*(page_width)/ 1.618
fig = plt.figure(figsize=(page_width,page_height)) #,edgecolor='black',linewidth=10)
spec2 = matplotlib.gridspec.GridSpec(ncols=1, nrows=2, figure=fig,wspace=.35,
                                 hspace=.1,bottom=.05,top=.975,left=.1,right=.975)
axes = []
f0 = fig.add_subplot(spec2[0, 0])
f1 = fig.add_subplot(spec2[1, 0])
axes.append([f0])
axes.append([f1])

def major_formatter(x, pos):
    return f'{x:.2f}'

def time_formatter(x, pos):
    return f'{x/60: .0f}'

freq = 2
partmc_ind_to_plot = np.arange(0,145,freq)
file = './out/camp' 
times = get_times(file, n_times)
var = "O3"
mix_rat = get_gas_mixing_ratio(file, n_times, var, n_runs)
pos = [0,0]
axes[pos[0]][pos[1]].plot(times[partmc_ind_to_plot], mix_rat[partmc_ind_to_plot] / 1000,
                     color=colors[0], ls='', marker='%s'%(part_symbol),label='%s' %(var))
axes[pos[0]][pos[1]].legend(handletextpad=0,labelspacing=.1,framealpha=1, edgecolor='k')

axes[pos[0]][pos[1]].grid(True,linestyle='--')
axes[pos[0]][pos[1]].set_xlim([0,24*60])
axes[pos[0]][pos[1]].xaxis.set_major_locator(ticker.FixedLocator([0,360,720,1080,1440]))
axes[pos[0]][pos[1]].xaxis.set_major_formatter(time_formatter)
axes[pos[0]][pos[1]].set_ylabel("Mixing ratio (ppm)")
axes[pos[0]][pos[1]].set_ylim([0,.14])
positions = [0,.02,.04,.06,.08,.10,.12,.14]
axes[pos[0]][pos[1]].yaxis.set_major_locator(ticker.FixedLocator(positions))
axes[pos[0]][pos[1]].yaxis.set_major_formatter(major_formatter)

pos = [1,0]
var = 'ISOP'
mix_rat = get_gas_mixing_ratio(file, n_times, var, n_runs)
axes[pos[0]][pos[1]].plot(times[partmc_ind_to_plot], mix_rat[partmc_ind_to_plot] / 1000,
                     color=colors[0], ls='', marker='%s' %(part_symbol), label='ISOP')
var = 'ISOP-P1'
mix_rat = get_gas_mixing_ratio(file, n_times, var, n_runs)
axes[pos[0]][pos[1]].plot(times[partmc_ind_to_plot], mix_rat[partmc_ind_to_plot] / 1000,
                     color=colors[1], ls='', marker='%s'%(part_symbol), label='ISOP-P1')
var = 'ISOP-P2'
mix_rat = get_gas_mixing_ratio(file, n_times, var, n_runs)
axes[pos[0]][pos[1]].plot(times[partmc_ind_to_plot], mix_rat[partmc_ind_to_plot] / 1000,
                     color=colors[2], ls='', marker='%s'%(part_symbol), label='ISOP-P2')

axes[pos[0]][pos[1]].legend(ncol=2, handletextpad=0, labelspacing=.1,
           framealpha=1, edgecolor='k')

ymax = [.005]
axes[pos[0]][pos[1]].grid(True,linestyle='--')
axes[pos[0]][pos[1]].set_xlim([0,24*60])
axes[pos[0]][pos[1]].xaxis.set_major_locator(ticker.FixedLocator([0,360,720,1080,1440]))
axes[pos[0]][pos[1]].xaxis.set_major_formatter(time_formatter)
axes[pos[0]][pos[1]].set_ylabel("Mixing ratio (ppm)")
axes[pos[0]][pos[1]].set_ylim([0,ymax[0]])
axes[pos[0]][pos[1]].yaxis.set_major_formatter(major_formatter)
out_filename =  'out/camp_gases.pdf'
print("Writing %s" % out_filename)
fig.savefig(out_filename)
