#!/usr/bin/env python
# Copyright (C) 2007 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys, re, textwrap
import copy as module_copy
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

class sum(reducer):
    def __init__(self, dim_name):
	self.dim_name = dim_name # dimension name to reduce

    def reduce(self, d):
	i_dim = d.find_dim_by_name(self.dim_name)
	new_data_shape = list(shape(d.data))
	del new_data_shape[i_dim]
	new_data = zeros(new_data_shape, dtype = float64)
	for i_val in range(len(d.dims[i_dim].grid_widths)):
	    reduce_slice = [slice(None, None, None) for i in range(len(d.dims))]
	    reduce_slice[i_dim] = i_val
	    reduce_slice = tuple(reduce_slice)
	    new_data += d.data[reduce_slice] * d.dims[i_dim].grid_widths[i_val]
	del d.dims[i_dim]
	d.data = new_data

class data_dim:
    def __init__(self):
	self.name = None         # axis name
	self.grid_type = None    # "center", "edge", or "center_edge"
	self.data_type = None    # "string", "integer", or "real"
	self.unit = None         # unit of data in grid cells
	self.grid_centers = None # grid cell centers on axis
	self.grid_edges = None   # grid cell edges on axis
	self.grid_widths = []    # widths of grid cells
	self.grid_units = None   # units of data values or None

    def __cmp__(self, other):
	val = cmp(self.name, other.name)
	if val != 0:
	    return val
	val = cmp(self.grid_type, other.grid_type)
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

    def read(self, f):
	self.name = read_string(f, 'name')
	self.grid_type = read_string(f, 'grid_type')
	if self.grid_type not in ["center", "edge", "center_edge"]:
	    raise Exception("unknown grid_type: %s" % self.grid_type)
	self.data_type = read_string(f, 'data_type')
	self.unit = read_string(f, 'unit')
	have_grid_units = read_string(f, 'have_grid_units')
	if have_grid_units not in ["yes", "no"]:
	    raise Exception("unknown value for have_grid_units: %s"
			    % have_grid_units)
	length = read_integer(f, 'length')
	if self.grid_type in ["center", "center_edge"]:
	    self.grid_centers = []
	    for i_grid in range(length):
		if self.data_type == "real":
		    grid_center = read_indexed_real(f, i_grid + 1,
						    'grid_center')
		elif self.data_type == "string":
		    grid_center = read_indexed_string(f, i_grid + 1,
						      'grid_center')
		elif self.data_type == "integer":
		    grid_center = read_indexed_integer(f, i_grid + 1,
						       'grid_center')
		else:
		    raise Exception("unknown data_type: %s" % self.data_type)
		self.grid_centers.append(grid_center)
	if self.grid_type in ["edge", "center_edge"]:
	    self.grid_edges = []
	    for i_grid in range(length + 1):
		if self.data_type == "real":
		    grid_edge = read_indexed_real(f, i_grid + 1,
						  'grid_edge')
		elif self.data_type == "string":
		    grid_edge = read_indexed_string(f, i_grid + 1,
						    'grid_edge')
		elif self.data_type == "integer":
		    grid_edge = read_indexed_integer(f, i_grid + 1,
						     'grid_edge')
		else:
		    raise Exception("unknown data_type: %s" % self.data_type)
		self.grid_edges.append(grid_edge)
	for i_grid in range(length):
	    self.grid_widths.append(read_indexed_real(f, i_grid + 1,
						     'grid_width'))
	if have_grid_units == "yes":
	    self.grid_units = []
	    for i_grid in range(length):
		self.grid_units.append(read_indexed_string(f, i_grid + 1,
							  'grid_unit'))

    def write_summary(self, f):
	output = []
	output.append("dim: %s, %s grid of type %s"
		% (self.name, self.grid_type, self.data_type))
	if self.unit:
	    output.append(" with unit %s" % self.unit)
	else:
	    output.append(" with no unit")
	if self.grid_type not in ["center", "edge", "center_edge"]:
	    raise Exception("unknown grid_type: %s" % self.grid_type)
	if self.grid_type in ["center", "center_edge"]:
	    if len(self.grid_centers) != len(self.grid_widths):
		raise Exception("grid data length mismatch for center")
	    if len(self.grid_centers) < 1:
		raise Exception("grid_centers empty")
	if self.grid_type in ["edge", "center_edge"]:
	    if len(self.grid_edges) != len(self.grid_widths) + 1:
		raise Exception("grid data length mismatch for edge")
	    if len(self.grid_edges) < 2:
		raise Exception("grid_edges too short")
	if self.grid_units:
	    if len(self.grid_widths) != len(self.grid_units):
		raise Exception("grid unit length mismatch")
	all_grid_widths_same = True
	for i in range(1,len(self.grid_widths)):
	    if self.grid_widths[i] != self.grid_widths[0]:
		all_grid_widths_same = False
	output.append(", ")
	if self.data_type == "integer":
	    if self.grid_units:
		raise Exception("grid_units not allowed for data_type int")
	    if self.grid_type == "center":
		min_grid = self.grid_centers[0]
		max_grid = self.grid_centers[-1]
	    elif self.grid_type in ["edge", "center_edge"]:
		min_grid = self.grid_edges[0]
		max_grid = self.grid_edges[-1]
	    else:
		raise Exception("unknown grid_type: %s" % self.grid_type)
	    if all_grid_widths_same:
		output.append("%d grid cells from %d to %d "
			      "(all with width %g)\n"
			      % (len(self.grid_widths), min_grid,
				 max_grid, self.grid_widths[-1]))
	    else:
		output.append("%d grid cells from %d "
			      "(width %g) to %d (width %g)\n"
			      % (len(self.grid_widths),
				 min_grid, self.grid_widths[0],
				 max_grid, self.grid_widths[-1]))
	elif self.data_type == "real":
	    if self.grid_units:
		raise Exception("grid_units not allowed for data_type float")
	    if self.grid_type == "center":
		min_grid = self.grid_centers[0]
		max_grid = self.grid_centers[-1]
	    elif self.grid_type in ["edge", "center_edge"]:
		min_grid = self.grid_edges[0]
		max_grid = self.grid_edges[-1]
	    else:
		raise Exception("unknown grid_type: %s" % self.grid_type)
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
	elif self.data_type == "string":
	    if self.grid_type != "center":
		raise Exception("only center grid allowed for data_type string")
	    output.append("%d grid values: " % len(self.grid_widths))
	    if self.grid_units and all_grid_widths_same:
		output.append(", ".join(["%s (unit %s)"
				   % (self.grid_centers[i],
				      self.grid_units[i])
				   for i in range(len(self.grid_centers))]))
	    elif self.grid_units and not all_grid_widths_same:
		output.append(", ".join(["%s (width %g, unit %s)"
				   % (self.grid_centers[i],
				      self.grid_widths[i],
				      self.grid_units[i])
				   for i in range(len(self.grid_centers))]))
	    elif not self.grid_units and all_grid_widths_same:
		output.append(", ".join(["%s" % (self.grid_centers[i])
				   for i in range(len(self.grid_centers))]))
	    elif not self.grid_units and not all_grid_widths_same:
		output.append(", ".join(["%s (width %g)"
				   % (self.grid_centers[i],
				      self.grid_widths[i])
				   for i in range(len(self.grid_centers))]))
	    else:
		raise Exception("internal error")
	    if all_grid_widths_same:
		output.append(" (all with width %g)" % self.grid_widths[0])
	else:
	    raise Exception("unknown data_type: %s" % self.data_type)
	wrap = textwrap.TextWrapper(subsequent_indent = "     ")
	f.write(wrap.fill("".join(output)))
	f.write("\n")

    def find_grid_by_value(self, value):
	if self.data_type in ["integer", "real"]:
	    if self.grid_type in ["edge", "center_edge"]:
		if value < self.grid_edges[0] or value > self.grid_edges[-1]:
		    raise Exception("value %g out of range [%g,%g]"
				    % value, self.grid_edges[0],
				    self.grid_edges[-1])
		i_val = 0
		while i_val < len(self.grid_edges) - 2:
		    if value < self.grid_edges[i_val + 1]:
			break
		    i_val += 1
	    elif self.grid_type == "center":
		i_val = 0
		while i_val < len(self.grid_centers) - 1:
		    if abs(self.grid_centers[i_val + 1] - value) \
			    > abs(self.grid_centers[i_val] - value):
			break
		    i_val += 1
	    else:
		raise Exception("unknown grid_type: %s" % self.grid_type)
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
	if self.grid_centers:
	    self.grid_centers = [x * factor for x in self.grid_centers]
	if self.grid_edges:
	    self.grid_edges = [x * factor for x in self.grid_edges]

class data_set:
    def __init__(self):
	self.name = None         # name of the data-set
	self.dims = []           # list of data_dim() objects
	self.data = None         # n-D array of data

    def read(self, f):
	self.name = read_string(f, 'name')
	n_dim = read_integer(f, 'n_dim')
	self.dims = []
	for i_dim in range(n_dim):
	    check_dim = read_integer(f, 'dim')
	    if check_dim != i_dim + 1:
		raise Exception("ERROR: expected dim %d but got dim: %d" \
				% (i_dim + 1, check_dim))
	    dim = data_dim()
	    dim.read(f)
	    self.dims.append(dim)
	read_comment(f, 'data values follow, row major order')
	data_shape = [len(dim.grid_widths) for dim in self.dims]
	self.data = zeros(data_shape, dtype = float64)
	for index in ndindex(*data_shape):
	    self.data[index] = read_unnamed_real(f)

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
	if self.dims[0].grid_type not in ["center", "center_edge"]:
	    raise Exception("can only generate data_center_list for center "
			    "or center_edge data")
	data = []
	for i in range(size(self.data)):
	    d = self.data[i]
	    if strip_zero and d <= 0.0:
		d = None
	    data.append([self.dims[0].grid_centers[i], d])
	return data

    def data_edge_list(self, strip_zero = False):
	if len(self.dims) != 1:
	    raise Exception("can only generate list with exactly one dim")
	if self.dims[0].grid_type not in ["edge", "center_edge"]:
	    raise Exception("can only generate data_center_list for edge "
			    "or center_edge data")
	data = []
	for i in range(size(self.data)):
	    d = self.data[i]
	    if strip_zero and d <= 0.0:
		d = None
	    data.append([self.dims[0].grid_edges[i],
			 self.dims[0].grid_edges[i+1], d])
	return data

    def scale_dim(self, dim_name, factor):
	i_dim = self.find_dim_by_name(dim_name)
	self.dims[i_dim].scale(factor)

class timed_data_set(data_set):
    def __init__(self):
	self.time = None         # time (s) of the data-set
	self.index = None        # index of the data-set
	data_set.__init__(self)

    def read(self, f):
	self.time = read_real(f, 'time')
	self.index = read_integer(f, 'index')
	data_set.read(self, f)

    def write_summary(self, f):
	f.write("time: %g (index %d)\n" % (self.time, self.index))
	data_set.write_summary(self, f)

def read_named(f, name):
    line = f.readline()
    data = line.split()
    if len(data) != 2:
	raise Exception("expected 2 items on line: %s" % line)
    if name != data[0]:
	raise Exception("expected name %s on line: %s" % (name, line))
    return data[1]

def read_indexed_named(f, index, name):
    line = f.readline()
    data = line.split()
    if len(data) != 3:
	raise Exception("expected 3 items on line: %s" % line)
    if name != data[0]:
	raise Exception("expected name %s on line: %s" % (name, line))
    try:
	check_index = int(data[1])
    except:
	raise Exception("ERROR: unable to convert index to integer: %s"
			% data[1])
    if index != check_index:
	raise Exception("ERROR: expected index %d but got: %d" %
			(index, check_index))
    return data[2]

def read_integer(f, name):
    data = read_named(f, name)
    try:
	return int(data)
    except:
	raise Exception("unable to convert to integer: %s" % data)

def read_real(f, name):
    data = read_named(f, name)
    try:
	return float(data)
    except:
	raise Exception("unable to convert to real: %s" % data)

def read_string(f, name):
    return read_named(f, name)

def read_indexed_integer(f, index, name):
    data = read_indexed_named(f, index, name)
    try:
	return int(data)
    except:
	raise Exception("unable to convert to integer: %s" % data)

def read_indexed_real(f, index, name):
    data = read_indexed_named(f, index, name)
    try:
	return float(data)
    except:
	raise Exception("unable to convert to float: %s" % data)

def read_indexed_string(f, index, name):
    return read_indexed_named(f, index, name)

def read_comment(f, comment):
    line = f.readline()
    if line != ("# %s\n" % comment):
	raise Exception("expected comment %s but got %s: " (comment, line))

def read_unnamed_real(f):
    line = f.readline()
    try:
	return float(line)
    except:
# HACK/FIXME: dies of numbers like 0.2341398-101 due to lack of E, so return 0
	return 0
# end HACK/FIXME
	raise Exception("unable to convert to float: %s" % line)

def read_data_set(files, reducers):
    if not files:
	raise Exception("need some files")
    data = []
    for file in files:
	f = open(file)
	d = timed_data_set()
	d.read(f)
	f.close()
	if d.index in [prev_d.index for prev_d in data]:
	    raise Exception("index %d read twice" % d.index)
	d.reduce(reducers)
	data.append(d)
    for i in range(1,len(data)):
	if data[i].name != data[0].name:
	    raise Exception("name mismatch btween indices %d and %d"
			    % (data[i].index, data[0].index))
	if data[i].dims != data[0].dims:
	    raise Exception("dimension mismatch between indices %d and %d"
			    % (data[i].index, data[0].index))
    data.sort(lambda a, b: cmp(a.time, b.time))
    total_data = data_set()
    total_data.name = data[0].name
    total_data.dims = module_copy.deepcopy(data[0].dims)
    time_dim = data_dim()
    time_dim.name = "time"
    time_dim.grid_type = "center"
    time_dim.data_type = "real"
    time_dim.unit = "s"
    time_dim.grid_centers = [d.time for d in data]
    time_dim.grid_widths = [1.0 for d in data]
    time_dim.grid_units = None
    total_data.dims.insert(0, time_dim)
    data_shape = [len(data)] + list(shape(data[0].data))
    total_data.data = zeros(data_shape, dtype = float64)
    for i, d in enumerate(data):
	total_data.data[i,...] = d.data
    return total_data
