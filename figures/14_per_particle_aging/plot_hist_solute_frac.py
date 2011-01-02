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
aging_no3_fraction = np.zeros([len(particle_set)])
aging_so4_fraction = np.zeros([len(particle_set)])
aging_nh4_fraction = np.zeros([len(particle_set)])
aging_solute_fraction = np.zeros([len(particle_set)])

i_counter = 0
for id in particle_set.keys():
	emit_diam[i_counter] = particle_set[id].emit_diam * 1e6
	bc_frac_emit[i_counter] = particle_set[id].emit_bc_fraction
	emit_comp_vols[i_counter] = particle_set[id].emit_comp_vols
	emit_time[i_counter] = particle_set[id].emit_time / 3600.
	emit_s_crit[i_counter] = particle_set[id].emit_s_crit
	emit_kappa[i_counter] = particle_set[id].emit_kappa
	aging_diameter[i_counter] = particle_set[id].aging_diameter
	aging_h2o[i_counter] = particle_set[id].aging_h2o
	aging_no3_fraction[i_counter] = particle_set[id].aging_no3_fraction
	aging_so4_fraction[i_counter] = particle_set[id].aging_so4_fraction
	aging_nh4_fraction[i_counter] = particle_set[id].aging_nh4_fraction
	aging_solute_fraction[i_counter] = particle_set[id].aging_solute_fraction
#	print "emit_time ", i_counter, emit_time[i_counter]

	if particle_set[id].aging_time != -1:
		time_for_aging[i_counter] = (particle_set[id].aging_time - particle_set[id].emit_time) / 3600.	
	else:
		time_for_aging[i_counter] = -1
	i_counter = i_counter + 1

emit_morning = ((emit_time < 6.) & (bc_frac_emit > 0))
emit_afternoon = (((emit_time > 6.) & (emit_time < 12.)) & (bc_frac_emit > 0))
emit_night = ((emit_time > 12) & (bc_frac_emit > 0 ))

bc_containing = (bc_frac_emit > 0)

# 2D Histogram plot
x_axis = partmc.log_grid(min=1e-3,max=1e1,n_bin=70)
y_axis = partmc.linear_grid(min=0,max=1,n_bin=50)

hist2d = partmc.histogram_2d(emit_diam[bc_containing], aging_solute_fraction[bc_containing], 
			     x_axis, y_axis, weights = 1 / emit_comp_vols[bc_containing])

hist2d_morning = partmc.histogram_2d(emit_diam[emit_morning], aging_solute_fraction[emit_morning], 
			     x_axis, y_axis, weights = 1 / emit_comp_vols[emit_morning])

hist2d_afternoon = partmc.histogram_2d(emit_diam[emit_afternoon], aging_solute_fraction[emit_afternoon], 
			     x_axis, y_axis, weights = 1 / emit_comp_vols[emit_afternoon])

hist2d = hist2d * 1e-6
hist2d_morning = hist2d_morning * 1e-6
hist2d_afternoon = hist2d_afternoon * 1e-6

(figure, axes_array, cbar_axes_array) = mpl_helper.make_fig_array(1,1, figure_width=5,
								  top_margin=1, bottom_margin=0.45,
								  left_margin=1.07, right_margin=0.65,
								  vert_sep=0.3, horiz_sep=0.3,
								  colorbar="shared", colorbar_location="top")

axes = axes_array[0][0]
cbar_axes = cbar_axes_array[0]
p = axes.pcolor(x_axis.edges(), y_axis.edges(), hist2d.transpose(), 
		norm = matplotlib.colors.LogNorm(vmin=1e3, vmax=1e6), linewidths = 0.1)
axes.set_xscale("log")
axes.set_yscale("linear")
axes.set_ylabel(r"solute fraction")
axes.set_xlabel(r"dry diameter at emission $D$/ $\rm \mu m$")
axes.set_ylim(0,1)
axes.set_xlim(1e-3, 1e0)
axes.grid(True)
cbar = figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext(),
		       orientation='horizontal')
cbar_axes.xaxis.set_label_position('top')
cbar.set_label(r"number conc. $n(D, w_{\rm sol})$ / $\rm cm^{-3}$")

mpl_helper.remove_fig_array_axes(axes_array)
figure.savefig("aging_wc_hist_solute_fract_06_allday.pdf")

(figure, axes_array, cbar_axes_array) = mpl_helper.make_fig_array(1,1, figure_width=5,
								  top_margin=1, bottom_margin=0.45,
								  left_margin=1.07, right_margin=0.65,
								  vert_sep=0.3, horiz_sep=0.3,
								  colorbar="shared", colorbar_location="top")

axes = axes_array[0][0]
cbar_axes = cbar_axes_array[0]
p = axes.pcolor(x_axis.edges(), y_axis.edges(), hist2d_morning.transpose(), 
		norm = matplotlib.colors.LogNorm(vmin=1e3, vmax=1e6), linewidths = 0.1)
axes.set_xscale("log")
axes.set_yscale("linear")
axes.set_ylabel(r"solute fraction")
axes.set_xlabel(r"dry diameter at emission $D$/ $\rm \mu m$")
axes.set_ylim(0,1)
axes.set_xlim(1e-3, 1e0)
axes.grid(True)
cbar = figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext(),
		       orientation='horizontal')
cbar_axes.xaxis.set_label_position('top')
cbar.set_label(r"number conc. $n(D,w_{\rm sol})$ / $\rm cm^{-3}$")

mpl_helper.remove_fig_array_axes(axes_array)
figure.savefig("aging_wc_hist_solute_frac_06_morning.pdf")

(figure, axes_array, cbar_axes_array) = mpl_helper.make_fig_array(1,1, figure_width=5,
								  top_margin=1, bottom_margin=0.45,
								  left_margin=1.07, right_margin=0.65,
								  vert_sep=0.3, horiz_sep=0.3,
								  colorbar="shared", colorbar_location="top")

axes = axes_array[0][0]
cbar_axes = cbar_axes_array[0]
p = axes.pcolor(x_axis.edges(), y_axis.edges(), hist2d_afternoon.transpose(), 
		norm = matplotlib.colors.LogNorm(vmin=1e3, vmax=1e6), linewidths = 0.1)
axes.set_xscale("log")
axes.set_yscale("linear")
axes.set_ylabel(r"solute fraction")
axes.set_xlabel(r"dry diameter at emission $D$/ $\rm \mu m$")
axes.set_ylim(0,1)
axes.set_xlim(1e-3, 1e0)
axes.grid(True)
cbar = figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext(),
		       orientation='horizontal')
cbar_axes.xaxis.set_label_position('top')
cbar.set_label(r"number conc. $n(D,w_{\rm sol})$ / $\rm cm^{-3}$")

mpl_helper.remove_fig_array_axes(axes_array)
figure.savefig("aging_wc_hist_solute_frac_06_afternoon.pdf")			
