from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import matplotlib.offsetbox as offsetbox
import matplotlib.font_manager
from matplotlib import ticker

def get_mass_conc(fileprefix, n_times, species, n_repeats=1):
    mass_conc = np.zeros((n_times,n_repeats))
    for i_repeat in range(n_repeats):
        for i_time in range(n_times):
            file = '%s_%04i_%08i.nc' % (fileprefix,i_repeat+1,i_time+1)
            f = Dataset(file, "r", format="NETCDF4")
            aero_names = f.variables['aero_species'].names.split(",")
            spec_index = aero_names.index(species)
            aero_num_conc = f.variables['aero_num_conc'][:]
            aero_part_mass = f.variables['aero_particle_mass'][:]
            mass_conc[i_time,i_repeat] = np.sum(aero_part_mass[spec_index,:] * aero_num_conc)
            f.close()
    return (np.mean(mass_conc,axis=1),np.std(mass_conc,axis=1))

def get_times(fileprefix, n_times):
    times = np.zeros(n_times)
    for i_time in range(n_times):
        file = '%s_0001_%08i.nc' % (fileprefix,i_time+1)
        f = Dataset(file, "r", format="NETCDF4")
        times[i_time] = (f.variables['time'][:]) / 60
        f.close()
    return times

def get_aerosol_species(fileprefix):
    i_time = 1
    file = '%s%08i.nc' % (fileprefix,i_time)
    f = Dataset(file, "r", format="NETCDF4")
    
    return(f.variables['aero_species'].names.split(","))

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

colors = ['#66c2a5','#fc8d62','#8da0cb']

def major_formatter(x, pos):
    return f'{1e9*x:.1f}'

pos = [0,0]
mass_conc, mass_std = get_mass_conc(file, n_times,'organic_matter.ISOP-P1_aero',
     n_repeats=n_runs)
axes[pos[0]][pos[1]].plot(times[partmc_ind_to_plot], mass_conc[partmc_ind_to_plot],
             color=colors[0], ls='', marker='%s'%(part_symbol), label="ISOP_P1_aero",
             markerfacecolor='none')

axes[pos[0]][pos[1]].legend(framealpha=1, edgecolor='k') 
axes[pos[0]][pos[1]].grid(True,linestyle='--')
axes[pos[0]][pos[1]].set_xlim([0,24*60])
axes[pos[0]][pos[1]].xaxis.set_major_locator(ticker.FixedLocator([0,360,720,1080,1440]))
axes[pos[0]][pos[1]].xaxis.set_major_formatter(time_formatter)
axes[pos[0]][pos[1]].set_ylim([0,3e-9])
axes[pos[0]][pos[1]].yaxis.set_major_formatter(major_formatter)
axes[pos[0]][pos[1]].set_ylabel(r'Mass concentration ($\mu \rm g \, m^{-3}$)')
pos = [1,0]
mass_conc, mass_std = get_mass_conc(file, n_times,'organic_matter.ISOP-P2_aero',
     n_repeats=n_runs)
axes[pos[0]][pos[1]].plot(times[partmc_ind_to_plot], mass_conc[partmc_ind_to_plot],
             color=colors[0], ls='', marker='%s'%(part_symbol), label="ISOP-P2_aero",
             markerfacecolor='none')
axes[pos[0]][pos[1]].legend(framealpha=1, edgecolor='k')
axes[pos[0]][pos[1]].grid(True,linestyle='--')
axes[pos[0]][pos[1]].set_xlabel('Simulation time (h)')
axes[pos[0]][pos[1]].set_xlim([0,24*60])
axes[pos[0]][pos[1]].xaxis.set_major_locator(ticker.FixedLocator([0,360,720,1080,1440]))
axes[pos[0]][pos[1]].xaxis.set_major_formatter(time_formatter)
axes[pos[0]][pos[1]].set_ylabel(r'Mass concentration ($\mu \rm g \, m^{-3}$)')
axes[pos[0]][pos[1]].set_ylim([0,3e-10])
axes[pos[0]][pos[1]].yaxis.set_major_formatter(major_formatter)
out_filename = 'out/timeseries_aero.pdf'
print("Writing %s" % out_filename)
fig.savefig(out_filename)
