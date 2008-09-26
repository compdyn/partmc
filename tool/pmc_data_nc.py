#!/usr/bin/env python
# Copyright (C) 2007-2008 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys, re, textwrap
import copy as module_copy
import numpy, math
from numpy import *
from Scientific.IO.NetCDF import *

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

class aero_particle_t:

    def __init__(self, ncf, index, aero_data):
        self.aero_data = aero_data
        if "aero_comp_mass" not in ncf.variables.keys():
            raise Exception("aero_comp_mass variable not found in NetCDF file")
        self.masses = ncf.variables["aero_comp_mass"][:,index]
        if "aero_n_orig_part" not in ncf.variables.keys():
            raise Exception("aero_n_orig_part variable not found in NetCDF file")
	self.n_orig_part = ncf.variables["aero_n_orig_part"][index]
        if "aero_absorb_cross_sect" not in ncf.variables.keys():
            raise Exception("aero_absorb_cross_sect variable not found in NetCDF file")
	self.absorb_cross_sect = ncf.variables["aero_absorb_cross_sect"][index]
        if "aero_scatter_cross_sect" not in ncf.variables.keys():
            raise Exception("aero_scatter_cross_sect variable not found in NetCDF file")
	self.scatter_cross_sect = ncf.variables["aero_scatter_cross_sect"][index]
        if "aero_asymmetry" not in ncf.variables.keys():
            raise Exception("aero_asymmetry variable not found in NetCDF file")
	self.asymmetry = ncf.variables["aero_asymmetry"][index]
        if "aero_refract_shell_real" not in ncf.variables.keys():
            raise Exception("aero_refract_shell_real variable not found in NetCDF file")
	self.refract_shell_real = ncf.variables["aero_refract_shell_real"][index]
        if "aero_refract_shell_imag" not in ncf.variables.keys():
            raise Exception("aero_refract_shell_imag variable not found in NetCDF file")
	self.refract_shell_imag = ncf.variables["aero_refract_shell_imag"][index]
        if "aero_refract_core_real" not in ncf.variables.keys():
            raise Exception("aero_refract_core_real variable not found in NetCDF file")
	self.refract_core_real = ncf.variables["aero_refract_core_real"][index]
        if "aero_refract_core_imag" not in ncf.variables.keys():
            raise Exception("aero_refract_core_imag variable not found in NetCDF file")
	self.refract_core_imag = ncf.variables["aero_refract_core_imag"][index]
        if "aero_core_vol" not in ncf.variables.keys():
            raise Exception("aero_core_vol variable not found in NetCDF file")
	self.core_vol = ncf.variables["aero_core_vol"][index]
        if "aero_water_hyst_leg" not in ncf.variables.keys():
            raise Exception("aero_water_hyst_leg variable not found in NetCDF file")
	self.water_hyst_leg = ncf.variables["aero_water_hyst_leg"][index]
        if "aero_comp_vol" not in ncf.variables.keys():
            raise Exception("aero_comp_vol variable not found in NetCDF file")
	self.comp_vol = ncf.variables["aero_comp_vol"][index]
        if "aero_id" not in ncf.variables.keys():
            raise Exception("aero_id variable not found in NetCDF file")
	self.id = ncf.variables["aero_id"][index]
        if "aero_least_create_time" not in ncf.variables.keys():
            raise Exception("aero_least_create_time variable not found in NetCDF file")
	self.least_create_time = ncf.variables["aero_least_create_time"][index]
        if "aero_greatest_create_time" not in ncf.variables.keys():
            raise Exception("aero_greatest_create_time variable not found in NetCDF file")
	self.greatest_create_time = ncf.variables["aero_greatest_create_time"][index]

    def sum_by_species(self, array, include = None, exclude = None):
        if include:
            for species in include:
                if species not in self.aero_data.name:
                    raise Exception("unknown species: %s" % species)
            species_list = set(include)
        else:
            species_list = set(self.aero_data.name)
        if exclude:
            for species in exclude:
                if species not in self.aero_data.name:
                    raise Exception("unknown species: %s" % species)
            species_list -= set(exclude)
        if len(species_list) == len(self.aero_data.name):
            return numpy.sum(array)
        else:
            val = 0.0
            for species in species_list:
                val += array[self.aero_data.name.index(species)]
        return val
    
    def mass(self, include = None, exclude = None):
        return self.sum_by_species(self.masses, include, exclude)

    def volume(self, include = None, exclude = None):
        return self.sum_by_species(self.masses / self.aero_data.density,
                                   include, exclude)

    def moles(self, include = None, exclude = None):
        return self.sum_by_species(self.masses / self.aero_data.density
                                   * self.aero_density.molec_weight,
                                   include, exclude)

    def radius(self):
        return (self.volume() * 3.0/4.0 / math.pi)**(1.0/3.0)

    def dry_radius(self):
        return (self.volume(exclude = ["H2O"]) * 3.0/4.0 / math.pi)**(1.0/3.0)

    def diameter(self):
        return 2.0 * self.radius()

    def dry_diameter(self):
        return 2.0 * self.dry_radius()

def read_particles(ncf, ids = None):
    if "aero_particle" not in ncf.dimensions.keys():
        raise Exception("aero_particle dimension not found in NetCDF file")
    n_part = ncf.dimensions["aero_particle"]

    if ids == None:
        indices = range(n_part)
    else:
        particle_ids = ncf.variables["aero_id"].getValue()
        indices = []
        for i in ids:
            ind = (particle_ids == i).argmax()
            if (ind > 0) or (particle_ids[0] == i):
                indices.append(ind)

    aero_data = aero_data_t(ncf)
    particles = []
    for i in indices:
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

    def find(self, values):
        indices = ((values - self.min) * self.n_bin
                   / (self.max - self.min)).astype(int)
        indices = indices.clip(0, self.n_bin - 1)
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
        
class pmc_log_axis:

    def __init__(self, min, max, n_bin):
        self.min = float(min)
        self.max = float(max)
        self.n_bin = n_bin

    def scale(self, factor):
        self.min = self.min * factor
        self.max = self.max * factor

    def grid_size(self, index):
        return (log10(self.max) - log10(self.min)) / float(self.n_bin)

    def find(self, value):
        indices = ((log(value) - log(self.min)) * self.n_bin
                   / (log(self.max) - log(self.min))).astype(int)
        indices = indices.clip(0, self.n_bin - 1)
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

class aero_particle_array_t:

    def __init__(self, ncf):
        self.aero_data = aero_data_t(ncf)
        if "aero_comp_mass" not in ncf.variables.keys():
            raise Exception("aero_comp_mass variable not found in NetCDF file")
        self.masses = ncf.variables["aero_comp_mass"].getValue()
        self.n_particles = size(self.masses, 1)
        if "aero_n_orig_part" not in ncf.variables.keys():
            raise Exception("aero_n_orig_part variable not found in NetCDF file")
	self.n_orig_part = ncf.variables["aero_n_orig_part"].getValue()
        if "aero_absorb_cross_sect" not in ncf.variables.keys():
            raise Exception("aero_absorb_cross_sect variable not found in NetCDF file")
	self.absorb_cross_sect = ncf.variables["aero_absorb_cross_sect"].getValue()
        if "aero_scatter_cross_sect" not in ncf.variables.keys():
            raise Exception("aero_scatter_cross_sect variable not found in NetCDF file")
	self.scatter_cross_sect = ncf.variables["aero_scatter_cross_sect"].getValue()
        if "aero_asymmetry" not in ncf.variables.keys():
            raise Exception("aero_asymmetry variable not found in NetCDF file")
	self.asymmetry = ncf.variables["aero_asymmetry"].getValue()
        if "aero_refract_shell_real" not in ncf.variables.keys():
            raise Exception("aero_refract_shell_real variable not found in NetCDF file")
	self.refract_shell_real = ncf.variables["aero_refract_shell_real"].getValue()
        if "aero_refract_shell_imag" not in ncf.variables.keys():
            raise Exception("aero_refract_shell_imag variable not found in NetCDF file")
	self.refract_shell_imag = ncf.variables["aero_refract_shell_imag"].getValue()
        if "aero_refract_core_real" not in ncf.variables.keys():
            raise Exception("aero_refract_core_real variable not found in NetCDF file")
	self.refract_core_real = ncf.variables["aero_refract_core_real"].getValue()
        if "aero_refract_core_imag" not in ncf.variables.keys():
            raise Exception("aero_refract_core_imag variable not found in NetCDF file")
	self.refract_core_imag = ncf.variables["aero_refract_core_imag"].getValue()
        if "aero_core_vol" not in ncf.variables.keys():
            raise Exception("aero_core_vol variable not found in NetCDF file")
	self.core_vol = ncf.variables["aero_core_vol"].getValue()
        if "aero_water_hyst_leg" not in ncf.variables.keys():
            raise Exception("aero_water_hyst_leg variable not found in NetCDF file")
	self.water_hyst_leg = ncf.variables["aero_water_hyst_leg"].getValue()
        if "aero_comp_vol" not in ncf.variables.keys():
            raise Exception("aero_comp_vol variable not found in NetCDF file")
	self.comp_vol = ncf.variables["aero_comp_vol"].getValue()
        if "aero_id" not in ncf.variables.keys():
            raise Exception("aero_id variable not found in NetCDF file")
	self.id = ncf.variables["aero_id"].getValue()
        if "aero_least_create_time" not in ncf.variables.keys():
            raise Exception("aero_least_create_time variable not found in NetCDF file")
	self.least_create_time = ncf.variables["aero_least_create_time"].getValue()
        if "aero_greatest_create_time" not in ncf.variables.keys():
            raise Exception("aero_greatest_create_time variable not found in NetCDF file")
	self.greatest_create_time = ncf.variables["aero_greatest_create_time"].getValue()

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
            val = self.masses[index,:].copy()
        for i in range(len(species_list) - 1):
            index = self.aero_data.name.index(species_list[i + 1])
            if species_weights != None:
                val += self.masses[index,:] * species_weights[index]
            else:
                val += self.masses[index,:]
        return val
    
    def mass(self, include = None, exclude = None):
        return self.sum_mass_by_species(include = include, exclude = exclude)

    def volume(self, include = None, exclude = None):
        species_weights = 1.0 / self.aero_data.density
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

    def solute_kappa(self):
        species_weights = self.aero_data.kappa \
            / self.aero_data.density
        solute_volume_kappa = self.sum_mass_by_species(exclude = ["H2O"],
                                                       species_weights = species_weights)
        solute_volume = self.volume(exclude = ["H2O"])
        return solute_volume_kappa / solute_volume

    def kappa_rh(self, env_state):
        A = 4.0 * const["water_surf_eng"] * const["water_molec_weight"] \
            / (const["univ_gas_const"] * env_state.temperature \
               * const["water_density"])
        C = sqrt(4.0 * A**3 / 27.0)
        diam = self.diameter()
        kappa = self.solute_kappa()
        return C / sqrt(kappa * diam**3) + 1.0

def time_of_day_string(time_seconds):
    time_of_day = time_seconds % (24 * 3600.0)
    hours = int(time_of_day / 3600.0)
    minutes = int(time_of_day / 60.0) % 60
    seconds = int(time_of_day) % 60
    return "%02d:%02d" % (hours, minutes)

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
        if "gas_concentration" not in ncf.variables.keys():
            raise Exception("gas_concentration variable not found in NetCDF file")
        self.concentration = ncf.variables["gas_concentration"][:]

    def concentration_by_species(self, species):
        if species not in self.gas_data.name:
            raise Exception("unknown species: %s" % species)
        index = self.gas_data.name.index(species)
        return self.concentration[index]

def read_history(constructor, directory, filename_pattern):
    filenames = os.listdir(directory)
    data = []
    filename_re = re.compile(filename_pattern)
    for filename in filenames:
        if filename_re.search(filename):
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

def get_time_filename_list(dir, file_pattern):
    time_filename_list = []
    filenames = os.listdir(dir)
    if filenames == []:
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

const = {}
consts_file = open("../src/constants.f90")
in_const_t = False
found_const_t = False
start_re = re.compile("^ *type const_t *$")
end_re = re.compile("^ *end type const_t *$")
const_re = re.compile("^ *real[*]8 :: ([^ ]+) = ([-0-9.]+)d([-0-9]+) *$")
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
