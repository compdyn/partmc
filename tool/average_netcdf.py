#!/usr/bin/env python
# Copyright (C) 2007-2008 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys, datetime
import copy as module_copy
from optparse import OptionParser
from Scientific.IO.NetCDF import *
from numpy import *

var_dir_base = ["assignValue", "getValue", "typecode"]

def common_prefix(strings):
    if len(strings) < 1:
	return ""
    min_len = min([len(s) for s in strings])
    i = 0
    while i < min_len:
	char = strings[0][i]
	done = False
	for s in strings[1:]:
	    if s[i] != char:
		done = True
	if done:
	    break
	i += 1
    return strings[0][0:i]

def process_args():
    parser = OptionParser(usage = "usage: %prog [options] <datafiles>")
    parser.add_option("-o", "--output", metavar = "FILE",
		      help = "Filename to write output to. [default: " \
		      + "<datafiles_prefix>avg.nc]")
    (options, args) = parser.parse_args()
    if len(args) < 1:
	parser.print_help()
	print "ERROR: must give at least one datafile argument"
	sys.exit(1)
    data_filenames = args
    return (parser, data_filenames)

def add_netcdf(total, add, var_totals, var_counts):
    # don't handle global attributes, as they are special
    # make sure we have all the dimensions
    for dim_name in add.dimensions.keys():
	print "add dim: %s" % dim_name
	dim_size = add.dimensions[dim_name]
	print "dim size: %s" % dim_size
	if dim_name not in total.dimensions.keys():
	    print "dimension not present in total, adding..."
	    total.createDimension(dim_name, dim_size)
	else:
	    print "dimension present in total, checking size..."
	    total_dim_size = total.dimensions[dim_name]
	    if dim_size != total_dim_size:
		raise Exception("dimension size mismatch: %d versus %d" \
				% (dim_size, total_dim_size))
    # make sure we have all the variables
    for var_name in add.variables.keys():
	print "add var: %s" % var_name
	var = add.variables[var_name]
	var_type = var.typecode()
	print "var_type: %s" % var_type
	var_shape = var.shape
	print "var_shape: %s" % string(var_shape)
	var_dims = var.dimensions
	print "var_dims: %s" % string(var_dims)
	if var_name not in total.variables.keys():
	    # variable not present in total, add it
	    print "variable not present in total, adding..."
	    total_var = total.createVariable(var_name, var_type, var_dims)
	    print "adding attributes to new variable"
	    for att_name in dir(var):
		if att_name not in var_dir_base:
		    print "att_name: %s" % att_name
		    att_value = getattr(var, att_name)
		    print "att_value: %s" % att_value
		    setattr(total_var, att_name, att_value)
	    if var_name in add.dimensions.keys():
		print "variable is also a dimension, setting values..."
		var_value = var.getValue()
		total_var.assignValue(var_value)
	    else:
		print "variable not a dimension, adding to hashes..."
		var_totals[var_name] = array(var.getValue())
		var_counts[var_name] = 1
	else:
	    # variable already present in total, check it matches
	    print "variable already in total, checking..."
	    total_var = total.variables[var_name]
	    total_var_type = total_var.typecode()
	    if var_type != total_var_type:
		raise Exception("variable type mismatch: %s, %s versus %s"
				% (var_name, var_type, total_var_type))
	    total_var_shape = var.shape
	    if var_shape != total_var_shape:
		raise Exception("variable shape mismatch: %s, %s versus %s"
				% (var_name, string(var_shape),
				   string(total_var_shape)))
	    total_var_dims = var.dimensions
	    if var_dims != total_var_dims:
		raise Exception("variable dims mismatch: %s, %s versus %s"
				% (var_name, string(var_dims),
				   string(total_var_dims)))
	    for att_name in dir(var):
		if att_name not in var_dir_base:
		    if att_name not in dir(total_var):
			raise Exception("variable attribute missing: %s, %s"
					% (var_name, att_name))
		    att_value = getattr(var, att_name)
		    total_att_value = getattr(total_var, att_name)
		    if att_value != total_att_value:
			raise Exception("variable attribute mismatch: "
					+ "%s, %s, %s versus %s"
					% (var_name, att_name, att_value,
					   total_att_value))
	    if var_name in add.dimensions.keys():
		print "variable is also a dimension, checking values..."
		var_value = var.getValue()
		total_var_value = total_var.getValue()
		if any(array(var_value) != array(total_var_value)):
		    raise Exception("dimension variable value mismatch: "
				    + "%s, %s versus %s"
				    % (var_name, string(var_value),
				       string(total_var_value)))
	    else:
		print "variable not a dimension, adding to hashes..."
		var_totals[var_name] += array(var.getValue())
		var_counts[var_name] += 1
		print "new var_counts[%s] = %d" % (var_name,
						   var_counts[var_name])

def average_netcdfs(total_filename, sum_filenames):
    print "opening total netcdf file: %s" % total_filename
    total = NetCDFFile(total_filename, "w")
    setattr(total, "title", "PartMC averaged output file")
    setattr(total, "history", "%s created by average_netcdf.py"
	    % datetime.datetime.now().isoformat())
    var_totals = {}
    var_counts = {}
    for add_filename in sum_filenames:
	print "opening add netcdf file: %s" % add_filename
	add = NetCDFFile(add_filename, "r")
	add_netcdf(total, add, var_totals, var_counts)
	add.close()
    for var_name in var_totals.keys():
	var = total.variables[var_name]
	var.assignValue(var_totals[var_name] / var_counts[var_name])
    total.close()

def main():
    (parser, data_filenames) = process_args()
    if parser.values.output:
	total_filename = parser.values.output
    else:
	total_filename = common_prefix(data_filenames) + "avg.nc"
    average_netcdfs(total_filename, data_filenames)

if __name__ == "__main__":
    main()


#mc_data = pmc_var(NetCDFFile("out/golovin_mc_state_0001.nc"),
#		  "aero",
#		  [sum("aero_species")])
#data_slice = module_copy.deepcopy(mc_data)
