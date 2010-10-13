#!/usr/bin/env python

import os, sys
import scipy.io
import numpy as np

sys.path.append("/Users/nriemer/subversion/partmc/trunk/tool")
import partmc
import mpl_helper
import matplotlib
import pickle

class Struct(object): 
	def __init__(self): 
		pass

input = open('particle_set_wc_no_nh3_03.pkl', 'rb')
particle_set = pickle.load(input)
input.close()

emit_diam = np.zeros([len(particle_set)])
time_for_aging = np.zeros([len(particle_set)])
bc_frac_emit = np.zeros([len(particle_set)])
emit_comp_vols = np.zeros([len(particle_set)])
emit_time = np.zeros([len(particle_set)])
emit_s_crit = np.zeros([len(particle_set)])
emit_kappa = np.zeros([len(particle_set)])
aging_kappa = np.zeros([len(particle_set)])
aging_diameter = np.zeros([len(particle_set)])
aging_h2o = np.zeros([len(particle_set)])

i_counter = 0
for id in particle_set.keys():
	emit_diam[i_counter] = particle_set[id].emit_diam
	bc_frac_emit[i_counter] = particle_set[id].emit_bc_fraction
	emit_comp_vols[i_counter] = particle_set[id].emit_comp_vols
	emit_time[i_counter] = particle_set[id].emit_time / 3600.
	emit_s_crit[i_counter] = particle_set[id].emit_s_crit
	emit_kappa[i_counter] = particle_set[id].emit_kappa

#	print "emit_time ", i_counter, emit_time[i_counter]

	if particle_set[id].aging_time != -1:
		time_for_aging[i_counter] = (particle_set[id].aging_time - particle_set[id].emit_time) / 3600.	
	else:
		time_for_aging[i_counter] = -1
	i_counter = i_counter + 1

emit_morning = (emit_time < 6.)
emit_afternoon = ((emit_time > 6.) & (emit_time < 12.))
emit_night = (emit_time > 12)

# Scatter plot, colored by time of emission
(figure, axes_array, cbar_axes_array) = mpl_helper.make_fig_array(3,1,
								  top_margin=1, bottom_margin=0.45,
								  left_margin=1.07, right_margin=0.65,
								  vert_sep=0.3, horiz_sep=0.3,
								  colorbar="shared", colorbar_location="top")
axes = axes_array[2][0]
cbar_axes = cbar_axes_array[0]
p = axes.scatter(emit_diam[emit_morning], time_for_aging[emit_morning], c=emit_time[emit_morning], 
		 norm = matplotlib.colors.Normalize(vmin=0, vmax=12), s=2, linewidths=0)
axes.set_xscale("log")
axes.set_yscale("linear")
axes.set_xlim(1e-9, 1e-6)
axes.set_ylim(0, 15)
axes.set_ylabel(r"per-particle aging time / h")
axes.grid(True)
cbar = figure.colorbar(p, cax=cbar_axes, 
		       orientation='horizontal')
cbar_axes.xaxis.set_label_position('top')
cbar.set_label(r"time of emission / h")

axes = axes_array[1][0]
p = axes.scatter(emit_diam[emit_afternoon], time_for_aging[emit_afternoon], c=emit_time[emit_afternoon], 
		 norm = matplotlib.colors.Normalize(vmin=0, vmax=12), s=2, linewidths=0)
axes.set_xscale("log")
axes.set_yscale("linear")
axes.set_ylim(0, 15)
axes.set_ylabel(r"per-particle aging time / h")
axes.grid(True)

axes = axes_array[0][0]
p = axes.scatter(emit_diam[emit_night], time_for_aging[emit_night], c=emit_time[emit_night],
		 norm = matplotlib.colors.Normalize(vmin=0, vmax=12), s=2, linewidths=0)
axes.set_xscale("log")
axes.set_yscale("linear")
axes.set_ylim(0, 15)
axes.set_ylabel(r"per-particle aging time / h")
axes.set_xlabel(r"dry diameter at emission $D_{\rm dry}$ / $\rm \mu m$")
axes.grid(True)

mpl_helper.remove_fig_array_axes(axes_array)
figure.savefig("aging_wc_no_nh3_03.pdf")

# Scatter plot, colored by s_crit at time of emission
(figure, axes_array, cbar_axes_array) = mpl_helper.make_fig_array(3,1,
								  top_margin=1, bottom_margin=0.45,
								  left_margin=1.07, right_margin=0.65,
								  vert_sep=0.3, horiz_sep=0.3,
								  colorbar="shared", colorbar_location="top")
axes = axes_array[2][0]
cbar_axes = cbar_axes_array[0]
p = axes.scatter(emit_diam[emit_morning], time_for_aging[emit_morning], c=emit_kappa[emit_morning], 
		 norm = matplotlib.colors.Normalize(vmin=0, vmax=1), s=2, linewidths=0)
axes.set_xscale("log")
axes.set_yscale("linear")
axes.set_xlim(1e-9, 1e-6)
axes.set_ylim(0, 15)
axes.set_ylabel(r"per-particle aging time / h")
axes.grid(True)
cbar = figure.colorbar(p, cax=cbar_axes, 
		       orientation='horizontal')
cbar_axes.xaxis.set_label_position('top')
cbar.set_label(r"$\kappa$ at time of emission / $\%$")

axes = axes_array[1][0]
p = axes.scatter(emit_diam[emit_afternoon], time_for_aging[emit_afternoon], c=emit_kappa[emit_afternoon], 
		 norm = matplotlib.colors.Normalize(vmin=0, vmax=1), s=2, linewidths=0)
axes.set_xscale("log")
axes.set_yscale("linear")
axes.set_ylim(0, 15)
axes.set_ylabel(r"per-particle aging time / h")
axes.grid(True)

axes = axes_array[0][0]
p = axes.scatter(emit_diam[emit_night], time_for_aging[emit_night], c=emit_kappa[emit_night],
		 norm = matplotlib.colors.Normalize(vmin=0, vmax=1), s=2, linewidths=0)
axes.set_xscale("log")
axes.set_yscale("linear")
axes.set_ylim(0, 15)
axes.set_ylabel(r"per-particle aging time / h")
axes.set_xlabel(r"dry diameter at emission $D_{\rm dry}$ / $\rm \mu m$")
axes.grid(True)

mpl_helper.remove_fig_array_axes(axes_array)
figure.savefig("aging_kappa_wc_no_nh3_03.pdf")


# 2D Histogram plot
x_axis = partmc.log_grid(min=1e-9,max=1e-5,n_bin=70)
y_axis = partmc.linear_grid(min=0,max=15,n_bin=30)

hist2d = partmc.histogram_2d(emit_diam, time_for_aging, x_axis, y_axis, weights = 1 / emit_comp_vols)

(figure, axes_array, cbar_axes_array) = mpl_helper.make_fig_array(1,1, figure_width=10,
								  top_margin=1, bottom_margin=0.45,
								  left_margin=1.07, right_margin=0.65,
								  vert_sep=0.3, horiz_sep=0.3,
								  colorbar="shared", colorbar_location="top")

axes = axes_array[0][0]
cbar_axes = cbar_axes_array[0]
p = axes.pcolor(x_axis.edges(), y_axis.edges(), hist2d.transpose(), 
		norm = matplotlib.colors.LogNorm(vmin=1e8, vmax=1e11), linewidths = 0.1)
axes.set_xscale("log")
axes.set_yscale("linear")
axes.set_ylabel(r"aging time / h")
axes.set_xlabel(r"dry diameter at emission / $\rm \mu m$")
#axes.set_ylim(1e-3,1e2)
#axes.set_xlim(5e-3, 5)
axes.grid(True)
cbar = figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext(),
		       orientation='horizontal')
cbar_axes.xaxis.set_label_position('top')
cbar.set_label(r"number conc. $n(D,w)$ / $\rm cm^{-3}$")

mpl_helper.remove_fig_array_axes(axes_array)
figure.savefig("aging_wc_hist_no_nh3_03.pdf")

#print "keys ", particle_set.keys()
#for id in particle_set.keys():
#	print id, particle_set[id].emit_time, particle_set[id].aging_time
			
