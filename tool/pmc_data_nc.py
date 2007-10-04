#!/usr/bin/env python
# Copyright (C) 2007 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys, re, textwrap
import copy as module_copy
import numpy
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

def sum_range(data, widths, i_dim, i_start, i_end):
    new_data_shape = list(shape(data))
    del new_data_shape[i_dim]
    new_data = zeros(new_data_shape, dtype = float64)
    for i_val in range(i_start, i_end):
	reduce_slice = [s_[:] for i in range(len(shape(data)))]
	reduce_slice[i_dim] = s_[i_val]
	reduce_slice = tuple(reduce_slice)
	new_data += array(data[reduce_slice]) * widths[i_val]
    return new_data

class sum(reducer):
    def __init__(self, dim_name):
	self.dim_name = dim_name # dimension name to reduce

    def reduce(self, d):
	i_dim = d.find_dim_by_name(self.dim_name)
	if d.dims[i_dim].grid_widths == None:
	    raise Exception("cannot sum dimension without widths: %s"
			    % self.dim_name)
	d.data = sum_range(d.data, d.dims[i_dim].grid_widths, i_dim,
			   0, size(d.data, i_dim))
	del d.dims[i_dim]

class sum_above(reducer):
    def __init__(self, dim_name, value):
	self.dim_name = dim_name # dimension name to reduce
	self.value = value       # value to sum above

    def reduce(self, d):
	i_dim = d.find_dim_by_name(self.dim_name)
	if d.dims[i_dim].grid_widths == None:
	    raise Exception("cannot sum dimension without widths: %s"
			    % self.dim_name)
	i_val = d.dims[i_dim].find_grid_by_value(self.value)
	d.data = sum_range(d.data, d.dims[i_dim].grid_widths, i_dim,
			   i_val, size(d.data, i_dim))
	del d.dims[i_dim]

class sum_below(reducer):
    def __init__(self, dim_name, value):
	self.dim_name = dim_name # dimension name to reduce
	self.value = value       # value to sum below

    def reduce(self, d):
	i_dim = d.find_dim_by_name(self.dim_name)
	if d.dims[i_dim].grid_widths == None:
	    raise Exception("cannot sum dimension without widths: %s"
			    % self.dim_name)
	i_val = d.dims[i_dim].find_grid_by_value(self.value)
	d.data = sum_range(d.data, d.dims[i_dim].grid_widths, i_dim,
			   0, i_val + 1)
	del d.dims[i_dim]

class sum_between(reducer):
    def __init__(self, dim_name, value_low, value_high):
	self.dim_name = dim_name      # dimension name to reduce
	self.value_low = value_low    # value to sum above
	self.value_high = value_high  # value to sum below

    def reduce(self, d):
	i_dim = d.find_dim_by_name(self.dim_name)
	if d.dims[i_dim].grid_widths == None:
	    raise Exception("cannot sum dimension without widths: %s"
			    % self.dim_name)
	i_val_low = d.dims[i_dim].find_grid_by_value(self.value_low)
	i_val_high = d.dims[i_dim].find_grid_by_value(self.value_high)
	d.data = sum_range(d.data, d.dims[i_dim].grid_widths, i_dim,
			   i_val_low, i_val_high + 1)
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
	self.grid_centers = ncf.variables[name][:]
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
	    self.grid_edges = ncf.variables[edges_name][:]
	    if len(self.grid_edges) != length + 1:
		raise Exception("incorrect length of edges variable: %s" % name)
	self.grid_widths = None
	if widths_name in ncf.variables.keys():
	    self.grid_widths = ncf.variables[widths_name][:]
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
		raise Exception("value not found: %s" % str(self.value))
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
	self.data = ncf.variables[name][:]
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

    def data_2d_list(self, strip_zero = False, flip_axes = False):
	if len(self.dims) != 2:
	    raise Exception("can only generate 2d_list with exactly two dims")
	if self.dims[0].grid_edges == None:
	    raise Exception("cannot generate 2d_list without grid_edges")
	if self.dims[1].grid_edges == None:
	    raise Exception("cannot generate 2d_list without grid_edges")
	data_list = []
	max_val = self.data.max()
	for i in range(size(self.data, 0)):
	    for j in range(size(self.data, 1)):
		if (not strip_zero) or (self.data[i,j] > 0):
		    if flip_axes:
			data_list.append([self.dims[1].grid_edges[j],
					  self.dims[1].grid_edges[j+1],
					  self.dims[0].grid_edges[i],
					  self.dims[0].grid_edges[i+1],
					  self.data[i,j] / max_val])
		    else:
			data_list.append([self.dims[0].grid_edges[i],
					  self.dims[0].grid_edges[i+1],
					  self.dims[1].grid_edges[j],
					  self.dims[1].grid_edges[j+1],
					  self.data[i,j] / max_val])
	return data_list

    def scale_dim(self, dim_name, factor):
	i_dim = self.find_dim_by_name(dim_name)
	self.dims[i_dim].scale(factor)

    def scale(self, factor):
	self.data = self.data * factor
