#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

from pdb import set_trace

caseName = "out_comp1"
initial_time = 0
final_time = 3600

OutDir = "./" + caseName

def vol2rad(vol):
    return (3 * vol / 4 / np.pi) ** (1/3.0)

def aero_particle_volumn(aero_particle_mass, aero_density, aero_species_name2ind):
    vol = 0
    for spec_name in aero_species_name2ind:
        vol += aero_particle_mass[aero_species_name2ind[spec_name], :] / aero_density[aero_species_name2ind[spec_name]]
    return vol

def aero_particle_dry_volumn(aero_particle_mass, aero_density, aero_species_name2ind):
    vol = 0
    for spec_name in aero_species_name2ind:
        if spec_name == "H2O":
            continue
        vol += aero_particle_mass[aero_species_name2ind[spec_name], :] / aero_density[aero_species_name2ind[spec_name]]
    return vol

def aero_particle_diameter(aero_particle_mass, aero_density, aero_species_name2ind):
    vol = aero_particle_volumn(aero_particle_mass, aero_density, aero_species_name2ind)
    radius = vol2rad(vol)
    diameter = 2 * radius
    return diameter

def aero_particle_dry_diameter(aero_particle_mass, aero_density, aero_species_name2ind):
    vol = aero_particle_dry_volumn(aero_particle_mass, aero_density, aero_species_name2ind)
    radius = vol2rad(vol)
    diameter = 2 * radius
    return diameter

def plot_bar(hist, edges, ax, **kwargs):
    #mid = (edges[:-1] + edges[1:]) / 2
    """
    if log_show:
        min_interval = np.min(np.log10(edges[1:]) - np.log10(edges[:-1]))
        ax.bar(np.log10(mid), hist, width = min_interval * 3)
        ax.set_xlabel(xlabel + " (10^* " + xunits + ")")
    else:
        min_interval = np.min(edges[1:] - edges[:-1])
        ax.bar(mid, hist, width = min_interval * 1)
        ax.set_xlabel(xlabel + " (" + xunits + ")")
    """
    min_interval =  np.min(edges[1:] - edges[:-1])
    ax.bar(edges[:-1], hist / (edges[1:] - edges[:-1]), align = "edge", width = min_interval * 0.9, **kwargs)

 
def plot_bar_bck(hist, edges, ax, xlabel, xunits, log_show = False):
    mid = (edges[:-1] + edges[1:]) / 2
    if log_show:
        min_interval = np.min(np.log10(edges[1:]) - np.log10(edges[:-1]))
        ax.bar(np.log10(mid), hist, width = min_interval * 3)
        ax.set_xlabel(xlabel + " (10^* " + xunits + ")")
    else:
        min_interval = np.min(edges[1:] - edges[:-1])
        ax.bar(mid, hist, width = min_interval * 1)
        ax.set_xlabel(xlabel + " (" + xunits + ")")
    
    #ax.stairs(hist, np.log10(edges))

def draw_hist(time, ax1, ax2, ax3, bins = 100, log_show = False):
    time_in_fileName = int(time / 100) + 1
    fileName = OutDir + "/freezing_part_0001_" + str(time_in_fileName).zfill(8) + ".nc"
    ncf = nc.Dataset(fileName)
    time_ = float(ncf.variables["time"][:].filled(np.nan))
    assert np.abs(time - time_) < 10e-6
    aero_particle_mass = ncf.variables["aero_particle_mass"][:].filled(np.nan) # (aero_species, aero_particle)
    aero_density = ncf.variables["aero_density"][:].filled(np.nan) #(aero_species) 
    aero_num_conc = ncf.variables["aero_num_conc"][:].filled(np.nan) # (aero_particle)
    p_frozen = ncf.variables["aero_frozen_probability"][:].filled(np.nan)
    aero_particle_vol = (aero_particle_mass.T / aero_density).T
    aero_species_ind = ncf.variables["aero_species"][:].filled(np.nan) - 1
    aero_species_names = ncf.variables["aero_species"].names.split(",")
    aero_species_ind2name = {ind:name for ind, name in zip(aero_species_ind, aero_species_names)}
    aero_species_name2ind = {name:ind for ind, name in zip(aero_species_ind, aero_species_names)}

    diameter = aero_particle_diameter(aero_particle_mass, aero_density, aero_species_name2ind)
    dry_diameter = aero_particle_dry_diameter(aero_particle_mass, aero_density, aero_species_name2ind)

    diameter *= 10**6
    dry_diameter *= 10**6
    hists_dict = {}
    if log_show:
        hist, bin_edges = np.histogram(np.log10(diameter), weights=aero_num_conc, bins = bins)
        hists_dict["diameter"] = (hist, bin_edges)
        hist, bin_edges = np.histogram(np.log10(dry_diameter), weights=aero_num_conc, bins = bins)
        hists_dict["dry_diameter"] = (hist, bin_edges)
        hist, bin_edges = np.histogram(np.log10(dry_diameter), weights = aero_num_conc * p_frozen, bins = bins)
        hists_dict["ice_dry_diameter"] = (hist, bin_edges)

    else:
        hist, bin_edges = np.histogram(diameter, weights=aero_num_conc, bins = bins)
        hists_dict["diameter"] = (hist, bin_edges)
        hist, bin_edges = np.histogram(dry_diameter, weights=aero_num_conc, bins = bins)
        hists_dict["dry_diameter"] = (hist, bin_edges)
        hist, bin_edges = np.histogram(dry_diameter, weights=aero_num_conc * p_frozen, bins = bins)
        hists_dict["ice_dry_diameter"] = (hist, bin_edges)



    #ax1.bar((hists_dict["diameter"][1][:-1] + hists_dict["diameter"][1][1:]) / 2, hists_dict["diameter"][0])
    #ax2.bar((hists_dict["dry_diameter"][1][:-1] + hists_dict["dry_diameter"][1][1:]) / 2, hists_dict["dry_diameter"][0])
    plot_bar(hists_dict["diameter"][0], hists_dict["diameter"][1], ax1, color = "green")
    plot_bar(hists_dict["dry_diameter"][0], hists_dict["dry_diameter"][1], ax2, color = "green")
    plot_bar(hists_dict["ice_dry_diameter"][0], hists_dict["ice_dry_diameter"][1], ax3, color = "blue")
    ax1.grid()
    ax2.grid()
    ax3.grid()

    if log_show:
        ax1.set_xlabel("diameter (10^* um)")
        ax2.set_xlabel("dry diameter (10^* um)")
        ax3.set_xlabel("dry diameter (10^* um)")
    else:
        ax1.set_xlabel("diameter (um)")
        ax2.set_xlabel("dry diameter (um)")
        ax3.set_xlabel("dry diameter (um)")

    ax1.set_ylabel("Number concentration density (m^-3)")
    ax2.set_ylabel("Number concentration density (m^-3)")
    ax3.set_ylabel("Number concentration density of ice (m^-3)")



    #set_trace()

if __name__ == "__main__":

    fig = plt.figure(figsize = (15, 16))
    ax1 = fig.add_subplot(2, 3, 1)
    ax2 = fig.add_subplot(2, 3, 2)
    ax3 = fig.add_subplot(2, 3, 3)

    ax4 = fig.add_subplot(2, 3, 4)
    ax5 = fig.add_subplot(2, 3, 5)
    ax6 = fig.add_subplot(2, 3, 6)

    draw_hist(initial_time, ax1, ax2, ax3, log_show = True)
    ax1.set_title("Initial distribution")
    ax2.set_title("Initial distribution")
    ax3.set_title("Ice initial distribution")


    draw_hist(final_time, ax4, ax5, ax6, log_show = True)
    ax4.set_title("Final distribution")
    ax5.set_title("Final distribution")
    ax6.set_title("Ice final distribution")

    plt.tight_layout(pad=5, w_pad=5, h_pad=10)
    #plt.tight_layout()

    plt.show()




