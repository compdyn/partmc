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

input = open('particle_set_wc_06.pkl', 'rb')
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
	aging_kappa[i_counter] = particle_set[id].aging_kappa
	aging_diameter[i_counter] = particle_set[id].aging_diameter
	aging_h2o[i_counter] = particle_set[id].aging_h2o
	print "aging_h2o ", aging_h2o[i_counter]
	if particle_set[id].aging_time != -1:
		time_for_aging[i_counter] = (particle_set[id].aging_time - particle_set[id].emit_time) / 3600.	
	else:
		time_for_aging[i_counter] = -1
	i_counter = i_counter + 1

emit_morning = (emit_time < 6.)
emit_afternoon = ((emit_time > 6.) & (emit_time < 12.))
emit_night = (emit_time > 12)

# Scatter plot, colored by h2o content at time of emission
(figure, axes_array, cbar_axes_array) = mpl_helper.make_fig_array(3,1,
								  top_margin=1, bottom_margin=0.45,
								  left_margin=1.07, right_margin=0.65,
								  vert_sep=0.3, horiz_sep=0.3,
								  colorbar="shared", colorbar_location="top")
axes = axes_array[2][0]
cbar_axes = cbar_axes_array[0]
p = axes.scatter(emit_diam[emit_morning], time_for_aging[emit_morning], c=aging_h2o[emit_morning], 
		 norm = matplotlib.colors.LogNorm(vmin=1e-21, vmax=1e-18), s=2, linewidths=0)
axes.set_xscale("log")
axes.set_yscale("linear")
axes.set_xlim(1e-9, 1e-6)
axes.set_ylim(0, 15)
axes.set_ylabel(r"per-particle aging time / h")
axes.grid(True)
cbar = figure.colorbar(p, cax=cbar_axes, 
		       orientation='horizontal')
cbar_axes.xaxis.set_label_position('top')
cbar.set_label(r"$\rm H_2O$ at time of aging / $kg$")

axes = axes_array[1][0]
p = axes.scatter(emit_diam[emit_afternoon], time_for_aging[emit_afternoon], c=aging_h2o[emit_afternoon], 
		 norm = matplotlib.colors.LogNorm(vmin=1e-21, vmax=1e-18), s=2, linewidths=0)
axes.set_xscale("log")
axes.set_yscale("linear")
axes.set_ylim(0, 15)
axes.set_ylabel(r"per-particle aging time / h")
axes.grid(True)

axes = axes_array[0][0]
p = axes.scatter(emit_diam[emit_night], time_for_aging[emit_night], c=aging_h2o[emit_night],
		 norm = matplotlib.colors.LogNorm(vmin=1e-21, vmax=1e-18), s=2, linewidths=0)
axes.set_xscale("log")
axes.set_yscale("linear")
axes.set_ylim(0, 15)
axes.set_ylabel(r"per-particle aging time / h")
axes.set_xlabel(r"dry diameter at emission $D_{\rm dry}$ / $\rm \mu m$")
axes.grid(True)

mpl_helper.remove_fig_array_axes(axes_array)
figure.savefig("aging_h2o_at_aging_wc_06.pdf")
			
