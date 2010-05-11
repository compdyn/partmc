#!/usr/bin/env python

import scipy.io
import sys, math
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import partmc

matplotlib.rc('text', usetex = True)
matplotlib.rc('xtick.major', pad = 8)
matplotlib.rc('ytick.major', pad = 8)
matplotlib.rc('xtick', labelsize = 10)
matplotlib.rc('legend', fontsize = 10, borderpad = 0.7, borderaxespad = 1)
matplotlib.rc('font', size = 10, family = "serif",
              serif = ["Computer Modern Roman"])
matplotlib.rc('lines', linewidth = 0.5)
matplotlib.rc('patch', linewidth = 0.5)
matplotlib.rc('axes', linewidth = 0.5)

aerosol_species_tex = {
    "SO4": "SO$_4$",
    "NO3": "NO$_3$",
    "Cl": "Cl",
    "NH4": "NH$_4$",
    "MSA": "MSA",
    "ARO1": "ARO1",
    "ARO2": "ARO2",
    "ALK1": "ALK1",
    "OLE1": "OLE1",
    "API1": "API1",
    "API2": "API2",
    "LIM1": "LIM1",
    "LIM2": "LIM2",
    "CO3": "CO$_3$",
    "Na": "Na",
    "Ca": "Ca",
    "OIN": "OIN",
    "OC": "OC",
    "BC": "BC",
    "H2O": "H$_2$O",
    }

def make_fig(figure_width = 4,
             figure_height = None,
             axis_ratio = (1 + math.sqrt(5)) / 2, # golden ratio
             left_margin = 0.65,
             right_margin = 0.2,
             bottom_margin = 0.5,
             top_margin = 0.3,
             colorbar = False,
             colorbar_width = 0.15,
             colorbar_height_fraction = 0.8,
             colorbar_offset = 0.2):
    axis_width = figure_width - left_margin - right_margin
    axis_height = axis_width / axis_ratio
    figure_height = bottom_margin + axis_height + top_margin
    left_margin_fraction = left_margin / figure_width
    bottom_margin_fraction = bottom_margin / figure_height
    axis_width_fraction = axis_width / figure_width
    axis_height_fraction = axis_height / figure_height
    figure = plt.figure()
    figure.set_figwidth(figure_width)
    figure.set_figheight(figure_height)
    axes = figure.add_axes([left_margin_fraction,
                            bottom_margin_fraction,
                            axis_width_fraction,
                            axis_height_fraction])
    axes.grid(True)
    axes.grid(True, which = 'minor')
    axes.minorticks_on()
    if colorbar:
        cb_left_fraction = (left_margin + axis_width + colorbar_offset) / figure_width
        cb_bottom_fraction = (bottom_margin + axis_height * (1.0 - colorbar_height_fraction) / 2.0) / figure_height
        cb_width_fraction = colorbar_width / figure_width
        cb_height_fraction = axis_height * colorbar_height_fraction / figure_height
        colorbar_axes = figure.add_axes([cb_left_fraction,
                                         cb_bottom_fraction,
                                         cb_width_fraction,
                                         cb_height_fraction])
    else:
        colorbar_axes = None
    return (figure, axes, colorbar_axes)

#ncf = scipy.io.netcdf.netcdf_file("~/t/urban_plume_wc_100K_flat_0001_00000010.nc", 'r')
#ncf = scipy.io.netcdf.netcdf_file("~/t/urban_plume_wc_100K_flat_0001_00000010.nc", 'r')
ncf = scipy.io.netcdf.netcdf_file("../../scenarios/4_nucleate/out/urban_plume_wc_0001_00000100.nc", 'r')
particles = partmc.aero_particle_array_t(ncf)
ncf.close()

(figure, axes, colorbar_axes) = make_fig(figure_width=12, colorbar=True, right_margin=1)

#i_species = [0, 1, 5, 6, 7, 8, 9, 17, 18]
#i_species = [0, 1, 17, 18]
i_species = [0, 1, 3, 5, 6, 7, 8, 9, 17, 18, 19, 2]
#i_species = range(len(particles.aero_data.names))
#print [(i, particles.aero_data.names[i]) for i in i_species]

#n_particles = 1000
n_particles = len(particles.ids)
i_particles = np.random.randint(0, len(particles.ids), n_particles)
n_species = len(i_species)
#for i in range(n_particles):
#    plt.semilogy(range(n_species), particles.raw_masses[:,i], 'k-')
x = np.tile(np.atleast_2d(np.arange(n_species)).transpose(), n_particles)
y = particles.raw_masses[:, i_particles][i_species, :]
spec_masses = particles.raw_masses[i_species, :]
m = spec_masses / sum(spec_masses, 0)
y = spec_masses[:, i_particles]
y = m[:, i_particles]

m_mean = m.mean(1)
m_mean_tile = np.tile(np.atleast_2d(m_mean).transpose(), np.size(m, 1))
m_nice = np.where(m > 0.0, m, m_mean_tile)

d = np.log(m_nice).transpose()
#d = np.nan_to_num(d)
d = d - d.mean(0)
#(U,L,Vh) = np.linalg.svd(np.dot(d.transpose(), d))
#mt = np.dot(U.transpose(), m)
#mt = mt.transpose()
#mt = mt - mt.min(0)
#mt = mt / mt.max(0)
#mt = mt.transpose()

dry_diameters = particles.dry_diameters()
dry_diameters = dry_diameters[i_particles]
dry_diameters *= 1e6
ddry_min = min(dry_diameters)
ddry_max = max(dry_diameters)

#ddry_min = 10**math.floor(math.log10(ddry_min))
#ddry_max = 10**math.ceil(math.log10(ddry_max))

print ddry_min, ddry_max

cmap = matplotlib.cm.get_cmap("spectral")
norm = matplotlib.colors.LogNorm(vmin=ddry_min, vmax=ddry_max)

cols = np.log(dry_diameters)
cols = cols - min(cols)
cols = cols / max(cols)

for i in range(np.size(x,1)):
    col = cmap(norm(dry_diameters[i]))
    axes.semilogy(x[:,i], y[:,i], color=col, linestyle='-', marker='o', linewidth=1, alpha=0.9)

axes.set_xticks(range(n_species))
axes.set_xticklabels([aerosol_species_tex[particles.aero_data.names[i]] for i in i_species])

axes.set_xlim(0, n_species - 1)
axes.set_ylim(1e-5, 1)

axes.set_xlabel('aerosol species')
axes.set_ylabel('mass fraction')

#figure.colorbar(p, cax = colorbar_axes, format = matplotlib.ticker.LogFormatterMathtext())
#colorbar_axes.set_ylabel(r"number conc. $(\rm cm^{-3})$")
#colorbar_yaxis = colorbar_axes.get_yaxis()

cbar = matplotlib.colorbar.ColorbarBase(colorbar_axes, cmap=cmap, norm=norm, format=matplotlib.ticker.ScalarFormatter())
cbar.set_label(r"dry diameter $(\mu\rm m)$")

#cbar_yaxis = colorbar_axes.get_xaxis()
#cbar_yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
#cbar_yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%f"))

#bc = particles.masses(include = ["BC"])
#dry_mass = particles.masses(exclude = ["H2O"])
#bc_frac = bc / dry_mass

#dry_diameters = particles.dry_diameters()

#    plt.xlabel(r'particle size $i$')
#    plt.ylabel(r'number concentration $E[n^V_i(2)]$')
plt.savefig('figs/parallel_10.pdf')
