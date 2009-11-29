#!/usr/bin/env python
# Copyright (C) 2007-2009 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys, re, textwrap
import copy as module_copy
import numpy, math
import random as py_random
from numpy import *
import scipy.optimize
from Scientific.IO.NetCDF import *

class aero_data_t:

    def __init__(self, ncf):
        if "aero_species" not in ncf.variables.keys():
            raise Exception("aero_species variable not found in NetCDF file")
        if "names" not in dir(ncf.variables["aero_species"]):
            raise Exception("aero_species variable does not have names attribute")
        self.name = ncf.variables["aero_species"].names.split(",")
        if "aero_mosaic_index" not in ncf.variables.keys():
            raise Exception("aero_mosaic_index variable not found in NetCDF file")
        self.mosaic_index = ncf.variables["aero_mosaic_index"][:]
        if "aero_density" not in ncf.variables.keys():
            raise Exception("aero_density variable not found in NetCDF file")
        self.density = ncf.variables["aero_density"][:]
        if "aero_num_ions" not in ncf.variables.keys():
            raise Exception("aero_num_ions variable not found in NetCDF file")
        self.num_ions = ncf.variables["aero_num_ions"][:]
        if "aero_solubility" not in ncf.variables.keys():
            raise Exception("aero_solubility variable not found in NetCDF file")
        self.solubility = ncf.variables["aero_solubility"][:]
        if "aero_molec_weight" not in ncf.variables.keys():
            raise Exception("aero_molec_weight variable not found in NetCDF file")
        self.molec_weight = ncf.variables["aero_molec_weight"][:]
        if "aero_kappa" not in ncf.variables.keys():
            raise Exception("aero_kappa variable not found in NetCDF file")
        self.kappa = ncf.variables["aero_kappa"][:]

class env_state_t:

    def __init__(self, ncf):
        if "temperature" not in ncf.variables.keys():
            raise Exception("temperature variable not found in NetCDF file")
        self.temperature = float(ncf.variables["temperature"].getValue())
        if "relative_humidity" not in ncf.variables.keys():
            raise Exception("relative_humidity variable not found in NetCDF file")
        self.relative_humidity = float(ncf.variables["relative_humidity"].getValue())
        if "pressure" not in ncf.variables.keys():
            raise Exception("pressure variable not found in NetCDF file")
        self.pressure = float(ncf.variables["pressure"].getValue())
        if "longitude" not in ncf.variables.keys():
            raise Exception("longitude variable not found in NetCDF file")
        self.longitude = float(ncf.variables["longitude"].getValue())
        if "latitude" not in ncf.variables.keys():
            raise Exception("latitude variable not found in NetCDF file")
        self.latitude = float(ncf.variables["latitude"].getValue())
        if "altitude" not in ncf.variables.keys():
            raise Exception("altitude variable not found in NetCDF file")
        self.altitude = float(ncf.variables["altitude"].getValue())
        if "start_time_of_day" not in ncf.variables.keys():
            raise Exception("start_time_of_day variable not found in NetCDF file")
        self.start_time_of_day = float(ncf.variables["start_time_of_day"].getValue())
        if "start_day_of_year" not in ncf.variables.keys():
            raise Exception("start_day_of_year variable not found in NetCDF file")
        self.start_day_of_year = int(ncf.variables["start_day_of_year"].getValue())
        if "elapsed_time" not in ncf.variables.keys():
            raise Exception("elapsed_time variable not found in NetCDF file")
        self.elapsed_time = float(ncf.variables["elapsed_time"].getValue())
        if "height" not in ncf.variables.keys():
            raise Exception("height variable not found in NetCDF file")
        self.height = float(ncf.variables["height"].getValue())

    def A(self, const):
        return 4.0 * const["water_surf_eng"] * const["water_molec_weight"] \
               / (const["univ_gas_const"] * self.temperature \
                  * const["water_density"])

class gas_data_t:

    def __init__(self, ncf):
        if "gas_species" not in ncf.variables.keys():
            raise Exception("gas_species variable not found in NetCDF file")
        if "names" not in dir(ncf.variables["gas_species"]):
            raise Exception("gas_species variable does not have names attribute")
        self.name = ncf.variables["gas_species"].names.split(",")
        if "gas_mosaic_index" not in ncf.variables.keys():
            raise Exception("gas_mosaic_index variable not found in NetCDF file")
        self.mosaic_index = ncf.variables["gas_mosaic_index"][:]
        if "gas_molec_weight" not in ncf.variables.keys():
            raise Exception("gas_molec_weight variable not found in NetCDF file")
        self.molec_weight = ncf.variables["gas_molec_weight"][:]

class gas_state_t:

    def __init__(self, ncf):
        self.gas_data = gas_data_t(ncf)
        if "gas_mixing_ratio" not in ncf.variables.keys():
            raise Exception("gas_mixing_ratio variable not found in NetCDF file")
        self.mixing_ratio = ncf.variables["gas_mixing_ratio"][:]

    def mixing_ratio_by_species(self, species):
        if species not in self.gas_data.name:
            raise Exception("unknown species: %s" % species)
        index = self.gas_data.name.index(species)
        return self.mixing_ratio[index]

class aero_particle_array_t:

    def __init__(self, ncf = None, n_particles = None, aero_data = None,
                 include_ids = None, exclude_ids = None):
        if ncf == None:
            self.aero_data = aero_data
            self.n_particles = n_particles
            self.masses = zeros([len(aero_data.name), n_particles])
            self.n_orig_part = zeros(n_particles, 'int32')
            self.absorb_cross_sect = zeros(n_particles)
            self.scatter_cross_sect = zeros(n_particles)
            self.asymmetry = zeros(n_particles)
            self.refract_shell_real = zeros(n_particles)
            self.refract_shell_imag = zeros(n_particles)
            self.refract_core_real = zeros(n_particles)
            self.refract_core_imag = zeros(n_particles)
            self.core_vol = zeros(n_particles)
            self.water_hyst_leg = zeros(n_particles, 'int32')
            self.comp_vol = zeros(n_particles)
            self.id = zeros(n_particles, 'int32')
            self.least_create_time = zeros(n_particles)
            self.greatest_create_time = zeros(n_particles)
            return

        self.aero_data = aero_data_t(ncf)

        for (ncf_var, self_var) in [
            ("aero_particle_mass", "masses"),
            ("aero_n_orig_part", "n_orig_part"),
            ("aero_absorb_cross_sect", "absorb_cross_sect"),
            ("aero_scatter_cross_sect", "scatter_cross_sect"),
            ("aero_asymmetry", "asymmetry"),
            ("aero_refract_shell_real", "refract_shell_real"),
            ("aero_refract_shell_imag", "refract_shell_imag"),
            ("aero_refract_core_real", "refract_core_real"),
            ("aero_refract_core_imag", "refract_core_imag"),
            ("aero_core_vol", "core_vol"),
            ("aero_water_hyst_leg", "water_hyst_leg"),
            ("aero_comp_vol", "comp_vol"),
            ("aero_id", "id"),
            ("aero_least_create_time", "least_create_time"),
            ("aero_greatest_create_time", "greatest_create_time"),
            ]:
            if ncf_var not in ncf.variables.keys():
                raise Exception("%s variable not found in NetCDF file" % ncf_var)
            self.__dict__[self_var] = asarray(ncf.variables[ncf_var].getValue())

        if include_ids != None or exclude_ids != None:
            keep_indexes = [i for i in range(size(self.id)) \
                            if (include_ids != None and self.id[i] in include_ids) \
                            or (exclude_ids != None and self.id[i] not in exclude_ids)]
            self.masses = self.masses[:, keep_indexes]
            self.n_orig_part = self.n_orig_part[keep_indexes]
            self.absorb_cross_sect = self.absorb_cross_sect[keep_indexes]
            self.scatter_cross_sect = self.scatter_cross_sect[keep_indexes]
            self.asymmetry = self.asymmetry[keep_indexes]
            self.refract_shell_real = self.refract_shell_real[keep_indexes]
            self.refract_shell_imag = self.refract_shell_imag[keep_indexes]
            self.refract_core_real = self.refract_core_real[keep_indexes]
            self.refract_core_imag = self.refract_core_imag[keep_indexes]
            self.core_vol = self.core_vol[keep_indexes]
            self.water_hyst_leg = self.water_hyst_leg[keep_indexes]
            self.comp_vol = self.comp_vol[keep_indexes]
            self.id = self.id[keep_indexes]
            self.least_create_time = self.least_create_time[keep_indexes]
            self.greatest_create_time = self.greatest_create_time[keep_indexes]

        self.n_particles = size(self.masses, 1)

    def sum_mass_by_species(self, include = None, exclude = None,
                            species_weights = None):
        if include != None:
            for species in include:
                if species not in self.aero_data.name:
                    raise Exception("unknown species: %s" % species)
            species_list = set(include)
        else:
            species_list = set(self.aero_data.name)
        if exclude != None:
            for species in exclude:
                if species not in self.aero_data.name:
                    raise Exception("unknown species: %s" % species)
            species_list -= set(exclude)
        species_list = list(species_list)
        if len(species_list) == 0:
            raise Exception("no species left to sum over")
        index = self.aero_data.name.index(species_list[0])
        if species_weights != None:
            val = self.masses[index,:].copy() * species_weights[index]
        else:
            val = array(self.masses[index,:].copy())
        for i in range(len(species_list) - 1):
            index = self.aero_data.name.index(species_list[i + 1])
            if species_weights != None:
                val += self.masses[index,:] * species_weights[index]
            else:
                val += array(self.masses[index,:])
        return val
    
    def mass(self, include = None, exclude = None):
        return self.sum_mass_by_species(include = include, exclude = exclude)

    def volume(self, include = None, exclude = None):
        species_weights = 1.0 / array(self.aero_data.density)
        return self.sum_mass_by_species(include = include, exclude = exclude,
                                        species_weights = species_weights)

    def moles(self, include = None, exclude = None):
        species_weights = self.aero_data.molec_weight \
                          / self.aero_data.density
        return self.sum_mass_by_species(include = include, exclude = exclude,
                                        species_weights = species_weights)

    def radius(self):
        return (self.volume() * 3.0/4.0 / math.pi)**(1.0/3.0)

    def dry_radius(self):
        return (self.volume(exclude = ["H2O"]) * 3.0/4.0 / math.pi)**(1.0/3.0)

    def diameter(self):
        return 2.0 * self.radius()

    def dry_diameter(self):
        return 2.0 * self.dry_radius()

    def surface_area(self):
        return 4.0 * math.pi * self.radius()**2

    def dry_surface_area(self):
        return 4.0 * math.pi * self.dry_radius()**2

    def solute_kappa_simple(self):
        species_weights = array(self.aero_data.kappa) \
            / array(self.aero_data.density)
        solute_volume_kappa = self.sum_mass_by_species(exclude = ["H2O"],
                                                       species_weights = species_weights)
        solute_volume = self.volume(exclude = ["H2O"])
        return solute_volume_kappa / solute_volume

    def solute_kappa(self):
        i_water = self.aero_data.name.index("H2O")
        M_w = self.aero_data.molec_weight[i_water]
        rho_w = self.aero_data.density[i_water]
        species_weights = zeros([len(self.aero_data.name)])
        for i_spec in range(size(species_weights)):
            if i_spec == self.aero_data.name == "H2O":
                continue
            if self.aero_data.num_ions[i_spec] > 0:
                if self.aero_data.kappa[i_spec] != 0:
                    raise Exception("species has nonzero num_ions and kappa: %s" % self.name[i_spec])
                M_a = self.aero_data.molec_weight[i_spec]
                rho_a = self.aero_data.density[i_spec]
                species_weights[i_spec] = M_w * rho_a / (M_a * rho_w) \
                                          * self.aero_data.num_ions[i_spec]
            else:
                species_weights[i_spec] = self.aero_data.kappa[i_spec]
        species_weights /= array(self.aero_data.density)
        solute_volume_kappa = self.sum_mass_by_species(exclude = ["H2O"],
                                                       species_weights = species_weights)
        solute_volume = self.volume(exclude = ["H2O"])
        return solute_volume_kappa / solute_volume

    def critical_rh_approx(self, env_state, const):
        A = env_state.A(const)
        C = sqrt(4.0 * A**3 / 27.0)
        diam = self.dry_diameter()
        kappa = self.solute_kappa()
        S = C / sqrt(kappa * diam**3) + 1.0
        return S

    def critical_rh(self, env_state, const):
        kappa = self.solute_kappa()
        dry_diam = self.dry_diameter()
        return critical_rh(env_state, const, kappa, dry_diam)

    def critical_diameter(self, env_state, const):
        kappa = self.solute_kappa()
        dry_diam = self.dry_diameter()
        return critical_diameter(env_state, const, kappa, dry_diam)

    def bin_average(self, diameter_axis, dry_diameter = True):
        averaged_particles = aero_particle_array_t(n_particles = diameter_axis.n_bin,
                                                   aero_data = self.aero_data)
        if dry_diameter:
            diameter = self.dry_diameter()
        else:
            diameter = self.diameter()
        diameter_bin = diameter_axis.find(diameter)
        num_conc = zeros(diameter_axis.n_bin)
        masses_conc = zeros([self.masses.shape[0], diameter_axis.n_bin])
        for i in range(self.n_particles):
            b = diameter_bin[i]
            if diameter_axis.valid_bin(b):
                num_conc[b] += 1.0 / self.comp_vol[i]
                masses_conc[:,b] += self.masses[:,i] / self.comp_vol[i]
        for b in range(averaged_particles.n_particles):
            averaged_particles.comp_vol[b] = 1.0 / num_conc[b]
            averaged_particles.masses[:,b] = masses_conc[:,b] * averaged_particles.comp_vol[b]
        return averaged_particles

def critical_rh(env_state, const, kappa, dry_diam):
    A = env_state.A(const)
    crit_diam = critical_diameter(env_state, const, kappa, dry_diam)
    return (crit_diam**3 - dry_diam**3) \
        / (crit_diam**3 - dry_diam**3 * (1 - kappa)) \
        * exp(A / crit_diam)

def critical_diameter(env_state, const, kappa, dry_diam):
    A = env_state.A(const)
    c6 = ones_like(dry_diam)
    c5 = zeros_like(dry_diam)
    c4 = - 3.0 * dry_diam**3 * kappa / A
    c3 = - dry_diam**3 * (2.0 - kappa)
    c2 = zeros_like(dry_diam)
    c1 = zeros_like(dry_diam)
    c0 = dry_diam**6 * (1.0 - kappa)
    def f(d):
        return c6 * d**6 + c5 * d**5 + c4 * d**4 \
            + c3 * d**3 + c2 * d**2 + c1 * d + c0

    d1 = dry_diam
    f1 = f(d1)
    if any(f1 >= 0.0):
        raise Exception("initialization failure for d1")
    d2 = 2.0 * d1
    f2 = f(d2)
    for iteration in range(100):
        if all(f2 >= 0.0):
            break
        d2 = where(f2 <= 0, 2.0 * d2, d2)
        f2 = f(d2)
    else:
        raise Exception("initialization failure for d2")
    dc = zeros_like(d1)
    for i in range(len(kappa)):
        def fi(d):
            return c6[i] * d**6 + c5[i] * d**5 + c4[i] * d**4 \
                + c3[i] * d**3 + c2[i] * d**2 + c1[i] * d + c0[i]
        dc[i] = scipy.optimize.brentq(fi, d1[i], d2[i])
    return dc

class aero_removed_info_t:

    AERO_INFO_NONE = 0
    AERO_INFO_DILUTION = 1
    AERO_INFO_COAG = 2
    AERO_INFO_HALVED = 3

    def __init__(self, ncf):
        for (ncf_var, self_var) in [
            ("aero_removed_id", "id"),
            ("aero_removed_action", "action"),
            ("aero_removed_other_id", "other_id"),
            ]:
            if ncf_var not in ncf.variables.keys():
                raise Exception("%s variable not found in NetCDF file" % ncf_var)
            self.__dict__[self_var] = asarray(ncf.variables[ncf_var].getValue())

        #self.aero_removed_id = [int(i) for i in self.aero_removed_id]
        #self.aero_removed_action = [int(i) for i in self.aero_removed_action]
        #self.aero_removed_other_id = [int(i) for i in self.aero_removed_other_id]

        if (len(self.aero_removed_id) == 1) and (self.aero_removed_id[0] == 0):
            self.id = array([],'int32')
            self.action = array([],'int32')
            self.other_id = array([],'int32')

class pmc_axis:

    def __init__(self, min, max, n_bin):
        self.min = float(min)
        self.max = float(max)
        self.n_bin = n_bin

    def find_clipped(self, values):
        indices = self.find(values)
        indices = indices.clip(0, self.n_bin - 1)
        return indices

    def find_clipped_outer(self, values):
        indices = self.find(values)
        indices += 1
        indices = indices.clip(0, self.n_bin + 1)
        return indices

    def closest_edge(self, value):
        i = self.find_clipped(value)
        lower_edge = self.edge(i)
        upper_edge = self.edge(i + 1)
        if abs(value - lower_edge) < abs(value - upper_edge):
            return i
        else:
            return i + 1

    def edges(self):
        return array([self.edge(i) for i in range(self.n_bin + 1)])

    def centers(self):
        return array([self.center(i) for i in range(self.n_bin)])

class pmc_linear_axis(pmc_axis):

    def __init__(self, min, max, n_bin):
        self.min = float(min)
        self.max = float(max)
        self.n_bin = n_bin

    def scale(self, factor):
        self.min = self.min * factor
        self.max = self.max * factor

    def grid_size(self, index):
        return (self.max - self.min) / float(self.n_bin)

    def valid_bin(self, bin):
        if (bin >= 0) and (bin < self.n_bin):
            return True
        return False

    def find(self, values):
        indices = (floor((values - self.min) * self.n_bin
                         / (self.max - self.min))).astype(int)
        #indices = indices.clip(0, self.n_bin - 1)
        return indices

    def edge(self, index):
        if (index < 0) or (index > self.n_bin):
            raise Exception("index out of range: %d" % index)
        if index == self.n_bin:
            return self.max
        elif index == 0:
            return self.min
        else:
            return float(index) / float(self.n_bin) * (self.max - self.min) \
                   + self.min

    def center(self, index):
        if (index < 0) or (index >= self.n_bin):
            raise Exception("index out of range: %d" % index)
        return (float(index) + 0.5) / float(self.n_bin) \
               * (self.max - self.min) + self.min

    def half_sample(self):
        if self.n_bin % 2 != 0:
            raise Exception("n_bin must be an even number")
        return pmc_linear_axis(min = self.min, max = self.max,
                               n_bin = self.n_bin / 2)
        
class pmc_log_axis(pmc_axis):

    def __init__(self, min, max, n_bin):
        self.min = float(min)
        self.max = float(max)
        self.n_bin = n_bin

    def scale(self, factor):
        self.min = self.min * factor
        self.max = self.max * factor

    def grid_size(self, index):
        return (log10(self.max) - log10(self.min)) / float(self.n_bin)

    def valid_bin(self, bin):
        if (bin >= 0) and (bin < self.n_bin):
            return True
        return False

    def find(self, value):
        indices = (floor((log(value) - log(self.min)) * self.n_bin
                   / (log(self.max) - log(self.min)))).astype(int)
        #indices = indices.clip(0, self.n_bin - 1)
        return indices

    def edge(self, index):
        if (index < 0) or (index > self.n_bin):
            raise Exception("index out of range: %d" % index)
        if index == self.n_bin:
            return self.max
        elif index == 0:
            return self.min
        else:
            return math.exp(float(index) / float(self.n_bin)
                            * (log(self.max) - log(self.min))
                            + log(self.min))
        
    def center(self, index):
        if (index < 0) or (index >= self.n_bin):
            raise Exception("index out of range: %d" % index)
        return math.exp((float(index) + 0.5) / float(self.n_bin)
                        * (log(self.max) - log(self.min))
                        + log(self.min))

    def half_sample(self):
        if self.n_bin % 2 != 0:
            raise Exception("n_bin must be an even number")
        return pmc_log_axis(min = self.min, max = self.max,
                            n_bin = self.n_bin / 2)

def histogram_1d(x_values, x_axis, weights = None):
    """Make a 1D histogram.

    The histogram is of points at positions x_values[i] for each i.
    Example:
    >>> x_axis = pmc_data_nc.pmc_log_axis(min = 1e-8, max = 1e-5, n_bin = 70)
    >>> hist = histogram_1d(diam, x_axis, weights = 1 / comp_vol)
    >>> semilogx(x_axis.centers(), hist)
    """
    if weights is not None:
        if len(x_values) != len(weights):
            raise Exception("x_values and weights have different lengths")
    x_bins = x_axis.find(x_values)
    hist = numpy.zeros([x_axis.n_bin])
    for i in range(len(x_values)):
        if x_axis.valid_bin(x_bins[i]):
            value = 1.0 / x_axis.grid_size(x_bins[i])
            if weights is not None:
                value *= weights[i]
            hist[x_bins[i]] += value
    return hist

def histogram_2d(x_values, y_values, x_axis, y_axis, weights = None, only_positive = True):
    """Make a 2D histogram.

    The histogram is of points at positions (x_values[i], y_values[i])
    for each i. Example:
    >>> x_axis = pmc_data_nc.pmc_log_axis(min = 1e-8, max = 1e-5, n_bin = 70)
    >>> y_axis = pmc_data_nc.pmc_linear_axis(min = 0, max = 1, n_bin = 50)
    >>> hist = histogram_2d(diam, bc_frac, x_axis, y_axis, weights = 1 / comp_vol)
    >>> pcolor(x_axis.edges(), y_axis.edges(), hist.transpose(),
               norm = matplotlib.colors.LogNorm(), linewidths = 0.1)
    """
    if len(x_values) != len(y_values):
        raise Exception("x_values and y_values have different lengths")
    if weights is not None:
        if len(x_values) != len(weights):
            raise Exception("x_values and weights have different lengths")
    x_bins = x_axis.find(x_values)
    y_bins = y_axis.find(y_values)
    hist = numpy.zeros([x_axis.n_bin, y_axis.n_bin])
    for i in range(len(x_values)):
        if x_axis.valid_bin(x_bins[i]) and y_axis.valid_bin(y_bins[i]):
            value = 1.0 / (x_axis.grid_size(x_bins[i]) * y_axis.grid_size(y_bins[i]))
            if weights is not None:
                value *= weights[i]
            hist[x_bins[i], y_bins[i]] += value
    if only_positive:
        mask = numpy.ma.make_mask(hist <= 0.0)
        hist = numpy.ma.array(hist, mask = mask)
    return hist

def multival_2d(x_values, y_values, z_values, x_axis, y_axis, rand_arrange = True):
    """Make a 2D matrix with 0%/33%/66%/100% percentile values.

    The returned matrix represents z_values[i] at position
    (x_values[i], y_values[i]) for each i. Example:
    >>> x_axis = pmc_data_nc.pmc_log_axis(min = 1e-8, max = 1e-5, n_bin = 140)
    >>> y_axis = pmc_data_nc.pmc_linear_axis(min = 0, max = 1, n_bin = 100)
    >>> vals = multival_2d(diam, bc_frac, h2o, x_axis, y_axis)
    >>> pcolor(x_axis.edges(), y_axis.edges(), vals.transpose(),
               norm = matplotlib.colors.LogNorm(), linewidths = 0.1)
    """
    if len(x_values) != len(y_values):
        raise Exception("x_values and y_values have different lengths")
    if len(x_values) != len(z_values):
        raise Exception("x_values and z_values have different lengths")

    low_x_axis = x_axis.half_sample()
    low_y_axis = y_axis.half_sample()
    x_bins = low_x_axis.find(x_values)
    y_bins = low_y_axis.find(y_values)
    z = [[[] for j in range(low_y_axis.n_bin)]
         for i in range(low_x_axis.n_bin)]
    for i in range(len(x_values)):
        if low_x_axis.valid_bin(x_bins[i]) and low_y_axis.valid_bin(y_bins[i]):
            z[x_bins[i]][y_bins[i]].append(z_values[i])
    for x_bin in range(low_x_axis.n_bin):
        for y_bin in range(low_y_axis.n_bin):
            z[x_bin][y_bin].sort()
    grid = numpy.zeros([x_axis.n_bin, y_axis.n_bin])
    mask = numpy.zeros([x_axis.n_bin, y_axis.n_bin], bool)
    for x_bin in range(low_x_axis.n_bin):
        for y_bin in range(low_y_axis.n_bin):
            if len(z[x_bin][y_bin]) > 0:
                subs = [(0,0),(0,1),(1,0),(1,1)]
                if rand_arrange:
                    py_random.shuffle(subs)
                (sub_min, sub_max, sub_low, sub_high) = subs
                val_min = min(z[x_bin][y_bin])
                val_max = max(z[x_bin][y_bin])
                val_low = percentile(z[x_bin][y_bin], 0.3333)
                val_high = percentile(z[x_bin][y_bin], 0.6666)
                for (sub, val) in [(sub_min, val_min),
                                   (sub_max, val_max),
                                   (sub_low, val_low),
                                   (sub_high, val_high)]:
                    sub_i, sub_j = sub
                    i = x_bin * 2 + sub_i
                    j = y_bin * 2 + sub_j
                    grid[i,j] = val
                    mask[i,j] = True
    mask = logical_not(mask)
    vals = numpy.ma.array(grid, mask = mask)
    return vals

def pmc_histogram_2d(array, x_axis, y_axis, mask = None, inv_mask = None):
    data = []
    for i in range(x_axis.n_bin):
        for j in range(y_axis.n_bin):
            if ((mask == None) and (inv_mask == None) and (array[i,j] > 0.0)) \
                   or ((mask != None) and (mask[i,j] != 0)) \
                   or ((inv_mask != None) and (inv_mask[i,j] == 0)):
                data.append([x_axis.edge(i), x_axis.edge(i + 1),
                             y_axis.edge(j), y_axis.edge(j + 1),
                             array[i,j]])
    return data

def pmc_histogram_2d_multi(array_list, x_axis, y_axis, mask = None,
                           inv_mask = None):
    data = []
    for i in range(x_axis.n_bin):
        for j in range(y_axis.n_bin):
            if ((mask == None) and (inv_mask == None)
                and (array_list[0][i,j] > 0.0)) \
                or ((mask != None) and (mask[i,j] != 0)) \
                or ((inv_mask != None) and (inv_mask[i,j] == 0)):
                data_item = [x_axis.edge(i), x_axis.edge(i + 1),
                             y_axis.edge(j), y_axis.edge(j + 1)]
                for array in array_list:
                    data_item.append(array[i,j])
                data.append(data_item)
    return data

def time_of_day_string(time_seconds, separator = ":"):
    time_of_day = time_seconds % (24 * 3600.0)
    hours = int(time_of_day / 3600.0)
    minutes = int(time_of_day / 60.0) % 60
    seconds = int(time_of_day) % 60
    return "%02d%s%02d" % (hours, separator, minutes)

def read_history(constructor, directory, filename_pattern, print_progress = False):
    filenames = os.listdir(directory)
    filenames.sort()
    data = []
    filename_re = re.compile(filename_pattern)
    for filename in filenames:
        if filename_re.search(filename):
            if print_progress:
                print filename
            netcdf_filename = os.path.join(directory, filename)
            ncf = NetCDFFile(netcdf_filename)
            env_state = env_state_t(ncf)
            data.append([env_state.elapsed_time, constructor(ncf)])
            ncf.close()
    data.sort()
    return data

def read_any(constructor, directory, filename_pattern):
    filenames = os.listdir(directory)
    filename_re = re.compile(filename_pattern)
    for filename in filenames:
        if filename_re.search(filename):
            netcdf_filename = os.path.join(directory, filename)
            ncf = NetCDFFile(netcdf_filename)
            data = constructor(ncf)
            ncf.close()
            return data
    raise Exception("no NetCDF file found in %s matching %s"
                    % (directory, filename_pattern))

def get_filename_list(dir, file_pattern):
    filename_list = []
    filenames = os.listdir(dir)
    if len(filenames)  == 0:
        raise Exception("No files in %s match %s" % (dir, file_pattern))
    file_re = re.compile(file_pattern)
    for filename in filenames:
        match = file_re.search(filename)
        if match:
            output_key = match.group(1)
            full_filename = os.path.join(dir, filename)
            filename_list.append(full_filename)
    filename_list.sort()
    if len(filename_list) == 0:
        raise Exception("No files found in %s matching %s"
                        % (dir, file_pattern))
    return filename_list

def get_time_filename_list(dir, file_pattern):
    time_filename_list = []
    filenames = os.listdir(dir)
    if len(filenames) == 0:
        raise Exception("No files in %s match %s" % (dir, file_pattern))
    file_re = re.compile(file_pattern)
    for filename in filenames:
        match = file_re.search(filename)
        if match:
            output_key = match.group(1)
            netcdf_filename = os.path.join(dir, filename)
            ncf = NetCDFFile(netcdf_filename)
            env_state = env_state_t(ncf)
            time_filename_list.append([env_state.elapsed_time,
                                       netcdf_filename,
                                       output_key])
            ncf.close()
    time_filename_list.sort()
    if len(time_filename_list) == 0:
        raise Exception("No files found in %s matching %s"
                        % (dir, file_pattern))
    return time_filename_list

def find_nearest_time(time_indexed_data, search_time):
    min_diff = abs(search_time - time_indexed_data[0][0])
    min_i = 0
    for i in range(1,len(time_indexed_data)):
        diff = abs(search_time - time_indexed_data[i][0])
        if diff < min_diff:
            min_diff = diff
            min_i = i
    return min_i

def file_filename_at_time(time_filename_list, search_time):
    i = find_nearest_time(time_filename_list, search_time)
    return time_filename_list[i][1]

def load_constants(filename):
    const = {}
    consts_file = open(filename) # filename of "constants.f90"
    in_const_t = False
    found_const_t = False
    start_re = re.compile("^ *type const_t *$")
    end_re = re.compile("^ *end type const_t *$")
    const_re = re.compile("^ *real[(]kind=dp[)] :: ([^ ]+) = ([-0-9.]+)d([-0-9]+) *$")
    for line in consts_file:
        if in_const_t:
            match = const_re.search(line)
            if match:
                name = match.group(1)
                mantissa = float(match.group(2))
                exponent = float(match.group(3))
                const[name] = mantissa * 10.0**exponent
            if end_re.search(line):
                in_const_t = False
        else:
            if start_re.search(line):
                in_const_t = True
                found_const_t = True
    if not found_const_t:
        raise Exception("constants.f90 ended without finding const_t")
    if in_const_t:
        raise Exception("constants.f90 ended without finding end of const_t")
    return const

def nan_to_value(x, value):
    values = x.copy()
    values.fill(value)
    return where(isnan(x), values, x)

def inf_to_value(x, value):
    values = x.copy()
    values.fill(value)
    return where(isinf(x), values, x)

def delta(arr):
    return (arr[1:] - arr[:-1])

def mid(arr):
    return (arr[1:] + arr[:-1]) / 2.0

def filter_inf(plot_data):
    return [[x, y] for [x, y] in plot_data if isfinite(y)]

#def sign(x):
#    if x > 0:
#        return 1
#    elif x < 0:
#        return -1
#    return 0

def mean(x):
    return float(sum(x)) / len(x)

def chop_sign_data_helper(plot_data, chopped_data):
    if len(plot_data) == 0:
        return chopped_data
    if sign(plot_data[0][1]) == 0:
        return chop_sign_data_helper(plot_data[1:], chopped_data)
    i = 0
    while (i < len(plot_data)) \
              and (sign(plot_data[i][1]) == sign(plot_data[0][1])):
        i += 1
    chopped_data.append(plot_data[0:i])
    return chop_sign_data_helper(plot_data[i:], chopped_data)

def chop_sign_data(plot_data):
    return chop_sign_data_helper(plot_data, [])

def chop_sign_data_iterative(plot_data):
    chopped_data = []
    start = 0
    while start < len(plot_data):
        while (start < len(plot_data)) and (sign(plot_data[start][1]) == 0):
            start += 1
        i = 0
        while (start + i < len(plot_data)) \
                and (sign(plot_data[start + i][1])
                     == sign(plot_data[start][1])):
            i += 1
        if i > 0:
            chopped_data.append(plot_data[start:(start + i)])
        start += i
    return chopped_data

def smooth(x,window_len=10,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string   
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."

    if window_len<3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is not one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"

    s=numpy.r_[2*x[0]-x[window_len:1:-1],x,2*x[-1]-x[-1:-window_len:-1]]
    if window == 'flat':
        w=ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')

    y=numpy.convolve(w/w.sum(),s,mode='same')
    return y[window_len-1:-window_len+1]

def cumulative_plot_data(x, y_inc, start = 0.0, final = None):
    plot_data = []
    y = start
    for i in range(x.size):
        plot_data.append([x[i], y])
        y += y_inc[i]
        if (i == x.size - 1) and (final != None):
            y = final
        plot_data.append([x[i], y])
    return plot_data

def cumulative_hi_res(x, y_inc, start = 0.0, final = None,
                      min_x_step = None, min_y_step = None,
                      min_x_factor = None, min_y_factor = None):
    plot_data = []
    i = 0
    y = start
    plot_data.append([x[i], y])
    last_x, last_y = x[i], y
    for i in range(1,x.size):
        y += y_inc[i]
        if (i == x.size - 1) and (final != None):
            y = final
        if ((min_x_step != None) and (x[i] - last_x > min_x_step)) \
                or ((min_y_step != None) and (y - last_y > min_y_step)) \
                or ((min_x_factor != None) and (x[i]/last_x > min_x_factor)) \
                or ((min_y_factor != None) and (y/last_y > min_y_factor)) \
                or (i == x.size - 1):
            plot_data.append([x[i], y])
            last_x = x[i]
            last_y = y
    return plot_data

def percentile(data, p):
    # p in [0,1]
    # data must be sorted
    i = int(floor((len(data) - 1) * p + 0.5))
    return data[i]

def grid_plot_data(x, y, z, x_axis, y_axis, hi_x_axis, hi_y_axis):
    x_bin = x_axis.find(x)
    y_bin = y_axis.find(y)
    values = [[[] for j in range(y_axis.n_bin)]
              for i in range(x_axis.n_bin)]
    for i in range(len(x)):
        if x_axis.valid_bin(x_bin[i]) and y_axis.valid_bin(y_bin[i]):
            values[x_bin[i]][y_bin[i]].append(z[i])
    for x_bin in range(x_axis.n_bin):
        for y_bin in range(y_axis.n_bin):
            values[x_bin][y_bin].sort()
    grid = numpy.zeros([hi_x_axis.n_bin, hi_y_axis.n_bin])
    mask = numpy.zeros([hi_x_axis.n_bin, hi_y_axis.n_bin], int)
    for x_bin in range(x_axis.n_bin):
        for y_bin in range(y_axis.n_bin):
            if len(values[x_bin][y_bin]) > 0:
                subs = [(0,0),(0,1),(1,0),(1,1)]
                py_random.shuffle(subs)
                (sub_min, sub_max, sub_low, sub_high) = subs
                val_min = min(values[x_bin][y_bin])
                val_max = max(values[x_bin][y_bin])
                val_low = percentile(values[x_bin][y_bin], 0.3333)
                val_high = percentile(values[x_bin][y_bin], 0.6666)
                for (sub, val) in [(sub_min, val_min),
                                   (sub_max, val_max),
                                   (sub_low, val_low),
                                   (sub_high, val_high)]:
                    sub_i, sub_j = sub
                    i = x_bin * 2 + sub_i
                    j = y_bin * 2 + sub_j
                    grid[i,j] = val
                    mask[i,j] = 1
    plot_data = pmc_histogram_2d(grid, hi_x_axis, hi_y_axis, mask = mask)
    return plot_data
