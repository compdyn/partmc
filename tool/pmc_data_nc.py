#!/usr/bin/env python
# Copyright (C) 2007-2008 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys, re, textwrap
import copy as module_copy
import numpy, math
from numpy import *

class reducer:
    def __init__(self):
	pass

    def reduce(self, d):
	raise Exception("not implemented")

class select(reducer):
    def __init__(self, dim_name, value):
	self.dim_name = dim_name # dimension name to reduce
	self.value = value       # value to select

    def reduce(self, d):
	i_dim = d.find_dim_by_name(self.dim_name)
	i_val = d.dims[i_dim].find_grid_by_value(self.value)
	reduce_slice = [slice(None, None, None) for i in range(len(d.dims))]
	reduce_slice[i_dim] = i_val
	reduce_slice = tuple(reduce_slice)
	del d.dims[i_dim]
	d.data = d.data[reduce_slice]

def sum_range(data, widths, i_dim, i_range):
    new_data_shape = list(shape(data))
    del new_data_shape[i_dim]
    new_data = zeros(new_data_shape, dtype = float64)
    for i_val in i_range:
	reduce_slice = [s_[:] for i in range(len(shape(data)))]
	reduce_slice[i_dim] = s_[i_val]
	reduce_slice = tuple(reduce_slice)
	new_data += data[reduce_slice] * widths[i_val]
    return new_data

class sum(reducer):
    def __init__(self, dim_name, above = None, below = None, only = None,
                 without = None):
	self.dim_name = dim_name      # dimension name to reduce
	self.above = above            # value to sum above
	self.below = below            # value to sum below
	self.only = only              # list of values to sum
        self.without = without        # list of values to leave out of sum
	if only != None:
	    if (above != None) or (below != None) or (without != None):
		raise Exception("cannot provide above, below or without if"
				+ " only is is given: %s" % dim_name)

    def reduce(self, d):
	i_dim = d.find_dim_by_name(self.dim_name)
	if d.dims[i_dim].grid_widths == None:
	    raise Exception("cannot sum dimension without widths: %s"
			    % self.dim_name)
	if self.only != None:
	    i_range = []
	    for val in self.only:
		i_range.append(d.dims[i_dim].find_grid_by_value(val))
	else:
	    if self.above != None:
		i_val_low = d.dims[i_dim].find_grid_by_value(self.above)
	    else:
		i_val_low = 0
	    if self.below != None:
		i_val_high = d.dims[i_dim].find_grid_by_value(self.below)
	    else:
		i_val_high = size(d.data, i_dim) - 1
	    i_range = range(i_val_low, i_val_high + 1)
            if self.without != None:
                for val in self.without:
                    del i_range[d.dims[i_dim].find_grid_by_value(val)]
	d.data = sum_range(d.data, d.dims[i_dim].grid_widths, i_dim, i_range)
	del d.dims[i_dim]

class pmc_dim:
    def __init__(self):
	self.name = None         # axis name
	self.data_type = None    # "string" or "number"
	self.unit = None         # unit of data in grid cells
	self.grid_centers = None # grid cell centers on axis
	self.grid_edges = None   # grid cell edges on axis
	self.grid_widths = None  # widths of grid cells
	self.grid_units = None   # units of data values or None

    def __init__(self, ncf, name):
	if name not in ncf.dimensions.keys():
	    raise Exception("dimension name not found: %s" % name)
	self.name = name
	length = ncf.dimensions[name]
	if name not in ncf.variables.keys():
	    raise Exception("dimension variable not found: %s" % name)
	if "unit" not in dir(ncf.variables[name]):
	    raise Exception("unit attribute not found for dimension: %s" % name)
	self.unit = ncf.variables[name].unit
	edges_name = name + "_edges"
	widths_name = name + "_widths"
	self.grid_centers = array(ncf.variables[name][:])
	if length:
	    if len(self.grid_centers) != length:
		raise Exception("incorrect length of variable: %s" % name)
	else:
	    length = len(self.grid_centers)
	self.data_type = 'number'
	if "names" in dir(ncf.variables[name]):
	    self.grid_centers = ncf.variables[name].names.split(',')
	    self.data_type = 'string'
	    if edges_name in ncf.variables.keys():
		raise Exception("edges variable cannot be present for "
				" string dimension: %s" % name)
	self.grid_edges = None
	if edges_name in ncf.variables.keys():
	    self.grid_edges = array(ncf.variables[edges_name][:])
	    if len(self.grid_edges) != length + 1:
		raise Exception("incorrect length of edges variable: %s" % name)
	self.grid_widths = None
	if widths_name in ncf.variables.keys():
	    self.grid_widths = array(ncf.variables[widths_name][:])
	    if len(self.grid_widths) != length:
		raise Exception("incorrect length of widths variable: %s"
				% name)
	self.grid_units = None
	if "data_units" in dir(ncf.variables[name]):
	    self.grid_units = ncf.variables[name].data_units.split(',')
	    if len(self.grid_units) != length:
		raise Exception("incorrect length of units variable: %s" % name)

    def __cmp__(self, other):
	val = cmp(self.name, other.name)
	if val != 0:
	    return val
	val = cmp(self.data_type, other.data_type)
	if val != 0:
	    return val
	val = cmp(self.unit, other.unit)
	if val != 0:
	    return val
	val = cmp(self.grid_centers, other.grid_centers)
	if val != 0:
	    return val
	val = cmp(self.grid_edges, other.grid_edges)
	if val != 0:
	    return val
	val = cmp(self.grid_widths, other.grid_widths)
	if val != 0:
	    return val
	val = cmp(self.grid_units, other.grid_units)
	if val != 0:
	    return val
	return 0

    def write_summary(self, f):
	output = []
	output.append("dim: %s, type %s"
		% (self.name, self.data_type))
	if self.unit:
	    output.append(" with unit %s" % self.unit)
	all_grid_widths_same = True
	if self.grid_widths != None:
	    for i in range(1,len(self.grid_widths)):
		if self.grid_widths[i] != self.grid_widths[0]:
		    all_grid_widths_same = False
	output.append(", ")
	if self.data_type == "number":
	    if self.grid_units != None:
		raise Exception("grid_units not allowed for data_type number")
	    if self.grid_edges != None:
		min_grid = self.grid_edges[0]
		max_grid = self.grid_edges[-1]
	    else:
		min_grid = self.grid_centers[0]
		max_grid = self.grid_centers[-1]
	    if self.grid_widths != None:
		if all_grid_widths_same:
		    output.append("%d grid cells from %g to %g "
				  "(all with width %g)\n"
				  % (len(self.grid_widths), min_grid,
				     max_grid, self.grid_widths[-1]))
		else:
		    output.append("%d grid cells from %g (width %g) "
				  "to %g (width %g)\n"
				  % (len(self.grid_widths),
				     min_grid, self.grid_widths[0],
				     max_grid, self.grid_widths[-1]))
	    else:
		output.append("%d grid cells from %g to %g\n"
			      % (len(self.grid_centers), min_grid, max_grid))
	elif self.data_type == "string":
	    output.append("%d grid values: " % len(self.grid_centers))
	    grid_strs = []
	    for i in range(len(self.grid_centers)):
		suffixes = []
		if self.grid_units != None:
		    suffixes.append("unit %s" % self.grid_units[i])
		if self.grid_widths != None and not all_grid_widths_same:
		    suffixes.append("width %g" % self.grid_widths[i])
		str = self.grid_centers[i]
		if suffixes:
		    str += " (" + ", ".join(suffixes) + ")"
		grid_strs.append(str)
	    output.append(", ".join(grid_strs))
	    if self.grid_widths != None and all_grid_widths_same:
		output.append(" (all with width %g)" % self.grid_widths[0])
	else:
	    raise Exception("unknown data_type: %s" % self.data_type)
	wrap = textwrap.TextWrapper(subsequent_indent = "     ")
	f.write(wrap.fill("".join(output)))
	f.write("\n")

    def find_grid_by_value(self, value):
	if self.data_type == "number":
	    if self.grid_edges != None:
		if value < self.grid_edges[0] or value > self.grid_edges[-1]:
		    raise Exception("value %g out of range [%g,%g]"
				    % value, self.grid_edges[0],
				    self.grid_edges[-1])
		i_val = 0
		while i_val < len(self.grid_edges) - 2:
		    if value < self.grid_edges[i_val + 1]:
			break
		    i_val += 1
	    else:
		i_val = 0
		while i_val < len(self.grid_centers) - 1:
		    if abs(self.grid_centers[i_val + 1] - value) \
			    > abs(self.grid_centers[i_val] - value):
			break
		    i_val += 1
	elif self.data_type == "string":
	    i_val = 0
	    for val in self.grid_centers:
		if value == val:
		    break
		i_val += 1
	    else:
		raise Exception("value not found: %s" % str(value))
	else:
	    raise Exception("unknown data_type: %s" % self.data_type)
	return i_val
    
    def scale(self, factor):
	if self.grid_centers != None:
	    self.grid_centers = [x * factor for x in self.grid_centers]
	if self.grid_edges != None:
	    self.grid_edges = [x * factor for x in self.grid_edges]

class pmc_var:
    def __init__(self):
	self.name = None         # name of the data-set
	self.dims = []           # list of data_dim() objects
	self.unit = None         # unit of data values
	self.data = None         # n-D array of data

    def __init__(self, ncf, name, reducers):
	if name not in ncf.variables.keys():
	    raise Exception("variable not found: %s" % name)
	self.name = name
	self.unit = None
	if "unit" in dir(ncf.variables[name]):
	    self.unit = ncf.variables[name].unit
	self.dims = []
	for d in ncf.variables[name].dimensions:
	    self.dims.append(pmc_dim(ncf, d))
	self.data = array(ncf.variables[name][:])
	self.reduce(reducers)

    def write_summary(self, f):
	f.write("name: %s (%d dims)\n" % (self.name, len(self.dims)))
	for dim in self.dims:
	    dim.write_summary(f)
	f.write("data: %d elements in %s array\n"
		% (size(self.data),
		   " x ".join([str(i) for i in shape(self.data)])))

    def reduce(self, reducers):
	for r in reducers:
	    r.reduce(self)

    def find_dim_by_name(self, dim_name):
	i_dim = 0
	for dim in self.dims:
	    if dim_name == dim.name:
		break
	    i_dim += 1
	else:
	    raise Exception("dimension not found: %s" % dim_name)
	return i_dim

    def dim_by_name(self, dim_name):
        i_dim = self.find_dim_by_name(dim_name)
        return self.dims[i_dim]

    def data_center_list(self, strip_zero = False):
	if len(self.dims) != 1:
	    raise Exception("can only generate list with exactly one dim")
	if self.dims[0].grid_centers == None:
	    raise Exception("cannot generate center_list without grid_centers")
	data_list = []
	for i in range(size(self.data)):
	    d = self.data[i]
	    if strip_zero and d <= 0.0:
		d = None
	    data_list.append([self.dims[0].grid_centers[i], d])
	return data_list

    def data_edge_list(self, strip_zero = False):
	if len(self.dims) != 1:
	    raise Exception("can only generate list with exactly one dim")
	if self.dims[0].grid_edges == None:
	    raise Exception("cannot generate edge_list without grid_edges")
	data_list = []
	for i in range(size(self.data)):
	    d = self.data[i]
	    if strip_zero and d <= 0.0:
		d = None
	    data_list.append([self.dims[0].grid_edges[i],
			      self.dims[0].grid_edges[i+1], d])
	return data_list

    def data_2d_list(self, strip_zero = False, flip_axes = False,
		     min = None, max = None):
	if len(self.dims) != 2:
	    raise Exception("can only generate 2d_list with exactly two dims")
	if self.dims[0].grid_edges == None:
	    raise Exception("2d_list error: no grid_edges for dim 0")
	if self.dims[1].grid_edges == None:
	    raise Exception("2d_list error: no grid_edges for dim 1")
	data_list = []
	if max:
	    max_val = max
	else:
	    max_val = self.data.max()
	if min:
	    min_val = min
	else:
	    min_val = 0.0
	for i in range(size(self.data, 0)):
	    for j in range(size(self.data, 1)):
		if (not strip_zero) or (self.data[i,j] > 0):
		    scaled_data = (self.data[i,j] - min_val) \
			/ (max_val - min_val)
		    if scaled_data < 0.0:
			scaled_data = 0.0
		    if scaled_data > 1.0:
			scaled_data = 1.0
		    if flip_axes:
			data_list.append([self.dims[1].grid_edges[j],
					  self.dims[1].grid_edges[j+1],
					  self.dims[0].grid_edges[i],
					  self.dims[0].grid_edges[i+1],
					  scaled_data])
		    else:
			data_list.append([self.dims[0].grid_edges[i],
					  self.dims[0].grid_edges[i+1],
					  self.dims[1].grid_edges[j],
					  self.dims[1].grid_edges[j+1],
					  scaled_data])
	return data_list

    def scale_dim(self, dim_name, factor):
	i_dim = self.find_dim_by_name(dim_name)
	self.dims[i_dim].scale(factor)

    def scale(self, factor):
	self.data = self.data * factor

class aero_data_t:

    def __init__(self, ncf):
        if "aero_species" not in ncf.variables.keys():
            raise Exception("aero_species variable not found in NetCDF file")
        if "names" not in dir(ncf.variables["aero_species"]):
            raise Exception("aero_species variable does not have names attribute")
        self.name = ncf.variables["aero_species"].names.split(",")
        if "mosaic_index" not in ncf.variables.keys():
            raise Exception("mosaic_index variable not found in NetCDF file")
        self.mosaic_index = ncf.variables["mosaic_index"][:]
        if "density" not in ncf.variables.keys():
            raise Exception("density variable not found in NetCDF file")
        self.density = ncf.variables["density"][:]
        if "num_ions" not in ncf.variables.keys():
            raise Exception("num_ions variable not found in NetCDF file")
        self.num_ions = ncf.variables["num_ions"][:]
        if "solubility" not in ncf.variables.keys():
            raise Exception("solubility variable not found in NetCDF file")
        self.solubility = ncf.variables["solubility"][:]
        if "molec_weight" not in ncf.variables.keys():
            raise Exception("molec_weight variable not found in NetCDF file")
        self.molec_weight = ncf.variables["molec_weight"][:]
        if "kappa" not in ncf.variables.keys():
            raise Exception("kappa variable not found in NetCDF file")
        self.kappa = ncf.variables["kappa"][:]

class aero_particle_t:

    def __init__(self, ncf, index, aero_data):
        self.aero_data = aero_data
        if "aero_comp_mass" not in ncf.variables.keys():
            raise Exception("aero_comp_mass variable not found in NetCDF file")
        self.masses = ncf.variables["aero_comp_mass"][:,index]
        if "n_orig_part" not in ncf.variables.keys():
            raise Exception("n_orig_part variable not found in NetCDF file")
	self.absorb_cross_sect = ncf.variables["n_orig_part"][index]
        if "absorb_cross_sect" not in ncf.variables.keys():
            raise Exception("absorb_cross_sect variable not found in NetCDF file")
	self.absorb_cross_sect = ncf.variables["absorb_cross_sect"][index]
        if "scatter_cross_sect" not in ncf.variables.keys():
            raise Exception("scatter_cross_sect variable not found in NetCDF file")
	self.scatter_cross_sect = ncf.variables["scatter_cross_sect"][index]
        if "asymmetry" not in ncf.variables.keys():
            raise Exception("asymmetry variable not found in NetCDF file")
	self.asymmetry = ncf.variables["asymmetry"][index]
        if "refract_shell_real" not in ncf.variables.keys():
            raise Exception("refract_shell_real variable not found in NetCDF file")
	self.refract_shell_real = ncf.variables["refract_shell_real"][index]
        if "refract_shell_imag" not in ncf.variables.keys():
            raise Exception("refract_shell_imag variable not found in NetCDF file")
	self.refract_shell_imag = ncf.variables["refract_shell_imag"][index]
        if "refract_core_real" not in ncf.variables.keys():
            raise Exception("refract_core_real variable not found in NetCDF file")
	self.refract_core_real = ncf.variables["refract_core_real"][index]
        if "refract_core_imag" not in ncf.variables.keys():
            raise Exception("refract_core_imag variable not found in NetCDF file")
	self.refract_core_imag = ncf.variables["refract_core_imag"][index]
        if "core_vol" not in ncf.variables.keys():
            raise Exception("core_vol variable not found in NetCDF file")
	self.core_vol = ncf.variables["core_vol"][index]
        if "water_hyst_leg" not in ncf.variables.keys():
            raise Exception("water_hyst_leg variable not found in NetCDF file")
	self.water_hyst_leg = ncf.variables["water_hyst_leg"][index]
        if "comp_vol" not in ncf.variables.keys():
            raise Exception("comp_vol variable not found in NetCDF file")
	self.comp_vol = ncf.variables["comp_vol"][index]
        if "aero_id" not in ncf.variables.keys():
            raise Exception("aero_id variable not found in NetCDF file")
	self.id = ncf.variables["aero_id"][index]
        if "least_create_time" not in ncf.variables.keys():
            raise Exception("least_create_time variable not found in NetCDF file")
	self.least_create_time = ncf.variables["least_create_time"][index]
        if "greatest_create_time" not in ncf.variables.keys():
            raise Exception("greatest_create_time variable not found in NetCDF file")
	self.greatest_create_time = ncf.variables["greatest_create_time"][index]

    def species_mass(self, species):
        return self.masses[self.aero_data.name.index(species)]
        
    def mass(self):
        return numpy.sum(self.masses)

    def volume(self):
        return numpy.sum(self.masses / self.aero_data.density)

    def radius(self):
        return (self.volume() * 3.0/4.0 / math.pi)**(1.0/3.0)

    def comp_frac(self, a_species, b_species, frac_type):
        for s in a_species:
            if s not in self.aero_data.name:
                raise Exception("unknown species name: %s" % s)
        for s in b_species:
            if s not in self.aero_data.name:
                raise Exception("unknown species name: %s" % s)
        a_val = 0.0
        b_val = 0.0
        for (i, s) in enumerate(self.aero_data.name):
            if s in a_species:
                if frac_type == "mass":
                    a_val = a_val + self.masses[i]
                elif frac_type == "volume":
                    a_val = a_val + self.masses[i] / self.aero_density[i]
                elif frac_type == "mole":
                    a_val = a_val + self.masses[i] / self.molec_weight[i]
                else:
                    raise Exception("unknown frac_type: %s" % frac_type)
            if s in b_species:
                if frac_type == "mass":
                    b_val = b_val + self.masses[i]
                elif frac_type == "volume":
                    b_val = b_val + self.masses[i] / self.aero_density[i]
                elif frac_type == "mole":
                    b_val = b_val + self.masses[i] / self.molec_weight[i]
                else:
                    raise Exception("unknown frac_type: %s" % frac_type)
        if (a_val == 0.0) and (b_val == 0.0):
            return 0.0
        else:
            return b_val / (a_val + b_val)

def read_particles(ncf):
    if "aero_particle" not in ncf.dimensions.keys():
        raise Exception("aero_particle dimension not found in NetCDF file")
    n_part = ncf.dimensions["aero_particle"]
    aero_data = aero_data_t(ncf)
    particles = []
    for i in range(n_part):
        particles.append(aero_particle_t(ncf, i, aero_data))
    return particles

class pmc_axis:

    def __init__(self, min, max, n_bin):
        self.min = float(min)
        self.max = float(max)
        self.n_bin = n_bin

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

    def find(self, value):
        return int((value - self.min) * self.n_bin / (self.max - self.min))

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
        
class pmc_log_axis:

    def __init__(self, min, max, n_bin):
        self.min = float(min)
        self.max = float(max)
        self.n_bin = n_bin

    def scale(self, factor):
        self.min = self.min * factor
        self.max = self.max * factor

    def grid_size(self, index):
        return (math.log(self.max) - math.log(self.min)) / float(self.n_bin)

    def find(self, value):
        return int((math.log(value) - math.log(self.min)) * self.n_bin
                     / (math.log(self.max) - math.log(self.min)))

    def edge(self, index):
        if (index < 0) or (index > self.n_bin):
            raise Exception("index out of range: %d" % index)
        if index == self.n_bin:
            return self.max
        elif index == 0:
            return self.min
        else:
            return math.exp(float(index) / float(self.n_bin)
                            * (math.log(self.max) - math.log(self.min))
                            + math.log(self.min))

def pmc_histogram_2d(array, x_axis, y_axis):
    data = []
    for i in range(x_axis.n_bin):
        for j in range(y_axis.n_bin):
            if array[i,j] > 0.0:
                data.append([x_axis.edge(i), x_axis.edge(i + 1),
                             y_axis.edge(j), y_axis.edge(j + 1),
                             array[i,j]])
    return data
