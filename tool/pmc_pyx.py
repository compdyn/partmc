#!/usr/bin/env python
# Copyright (C) 2007 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys
sys.path.append(os.path.expanduser("~/.python"))
from pyx import *

color_list = [color.hsb(0/3.0, 1, 1),
	      color.hsb(1/3.0, 1, 1),
	      color.hsb(2/3.0, 1, 1),
	      color.hsb(1/6.0, 1, 1),
	      color.hsb(3/6.0, 1, 1),
	      color.hsb(5/6.0, 1, 1),
	      color.hsb(1/12.0, 1, 1),
	      color.hsb(3/12.0, 1, 1),
	      color.hsb(5/12.0, 1, 1),
	      color.hsb(7/12.0, 1, 1),
	      color.hsb(9/12.0, 1, 1),
	      color.hsb(11/12.0, 1, 1),
	      ]

grid_painter = graph.axis.painter.regular(gridattrs = [style.linestyle.dotted])

aerosol_species_tex = {
    "SO4_a": "SO$_4$",
    "NO3_a": "NO$_3$",
    "Cl_a": "Cl",
    "NH4_a": "NH$_4$",
    "MSA_a": "MSA",
    "ARO1_a": "ARO1",
    "ARO2_a": "ARO2",
    "ALK1_a": "ALK1",
    "OLE1_a": "OLE1",
    "API1_a": "API1",
    "API2_a": "API2",
    "LIM1_a": "LIM1",
    "LIM2_a": "LIM2",
    "CO3_a": "CO$_3$",
    "Na_a": "Na",
    "Ca_a": "Ca",
    "OIN_a": "OIN",
    "OC_a": "Organic carbon",
    "BC_a": "Black carbon",
    "H2O": "H$_2$O",
    }

gas_species_tex = {
    "H2SO4": "H$_2$SO$_4$",
    "HNO3": "HNO$_3$",
    "HCl": "HC$\ell$",
    "NH3": "NH$_3$",
    "NO2": "NO$_2$",
    "NO3": "NO$_3$",
    "N2O5": "N$_2$O$_5$",
    "HNO4": "HNO$_4$",
    "O3": "O$_3$",
    "O1D": "O$_1$D",
    "O3P": "O$_3$P",
    "HO2": "HO$_2$",
    "H2O2": "H$_2$O$_2$",
    "SO2": "SO$_2$",
    "CH4": "CH$_4$",
    "C2H6": "C$_2$H$_6$",
    "CH3O2": "CH$_3$O$_2$",
    "CH3OH": "CH$_3$OH",
    "CH3OOH": "CH$_3$OOH",
    "C2O3": "C$_2$O$_3$",
    "CH3SO2H": "CH$_3$SO$_2$H",
    "CH3SCH2OO": "CH$_3$SCH$_2$OO",
    "CH3SO2": "CH$_3$SO$_2$",
    "CH3SO3": "CH$_3$SO$_3$",
    "CH3SO2OO": "CH$_3$SO$_2$OO",
    "CH3SO2CH2OO": "CH$_3$SO$_2$CH$_2$OO",
    }

def tex_species(species):
    if species in aerosol_species_tex.keys():
	return aerosol_species_tex[species]
    if species in gas_species_tex.keys():
	return gas_species_tex[species]
    return species
