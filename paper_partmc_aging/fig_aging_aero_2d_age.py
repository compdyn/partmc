#!/usr/bin/env python
# Copyright (C) 2007-2009 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys, math, random
import copy as module_copy
from Scientific.IO.NetCDF import *
from pyx import *
sys.path.append(".")
from fig_helper import *
sys.path.append("../tool")
from pmc_data_nc import *
from pmc_pyx import *

out_prefix = "figs_aging/aging_aero_2d_age"

ss_axis_min = 1e-2
ss_axis_max = 1e2
num_ss_bins = 40

bc_max_val = 4.0
ss_max_val = 100.0

const = load_constants("../src/constants.f90")

class colorpoints(graph.style.symbol):
    
    def __init__(self, sizecolumnname="size", colorcolumnname="color",
                 gradient=color.gradient.Rainbow,
                 symbol=graph.style.symbol.circle,
                 symbolattrs=[deco.filled, deco.stroked([color.gray.black])],
                 **kwargs):
        # add some configuration parameters and modify some other
        self.sizecolumnname = sizecolumnname
        self.colorcolumnname = colorcolumnname
        self.gradient = gradient
        graph.style.symbol.__init__(self, symbol=symbol, symbolattrs=symbolattrs, **kwargs)

    def columnnames(self, privatedata, sharedata, agraph, columnnames):
        # register the new column names
        if self.sizecolumnname not in columnnames:
            raise ValueError("column '%s' missing" % self.sizecolumnname)
        if self.colorcolumnname not in columnnames:
            raise ValueError("column '%s' missing" % self.colorcolumnname)
        return ([self.sizecolumnname, self.colorcolumnname] +
                graph.style.symbol.columnnames(self, privatedata,
                                               sharedata, agraph, columnnames))

    def drawpoint(self, privatedata, sharedata, graph, point):
        # replace the original drawpoint method by a slightly revised one
        if sharedata.vposvalid and privatedata.symbolattrs is not None:
            x_pt, y_pt = graph.vpos_pt(*sharedata.vpos)
            color = self.gradient.getcolor(point[self.colorcolumnname])
            privatedata.symbol(privatedata.symbolcanvas, x_pt, y_pt,
                               privatedata.size_pt*point[self.sizecolumnname],
                               privatedata.symbolattrs + [color])

time_filename_list = get_time_filename_list(netcdf_dir_wc, netcdf_pattern_wc)
for color in [True, False]:
    g = graph.graphxy(
        width = grid_graph_width,
        x = graph.axis.log(min = diameter_axis_min,
                           max = diameter_axis_max,
                           title = diameter_axis_label),
        y = graph.axis.log(min = ss_axis_min,
                              max = ss_axis_max,
                              title = r"critical supersaturation $S_{\rm c}$ (\%)"))

    time = 24 * 3600.0
    filename = file_filename_at_time(time_filename_list, time)

    ncf = NetCDFFile(filename)
    particles = aero_particle_array_t(ncf)
    env_state = env_state_t(ncf)
    ncf.close()

    diameter = particles.dry_diameter() * 1e6
    critical_ss = (particles.kappa_rh(env_state, const) - 1.0) * 100.0
    age = (env_state.elapsed_time - particles.least_create_time) / 3600.0
    max_age = max(age)
    scaled_age = age / max_age

    plot_data = zip(diameter, critical_ss, ones_like(diameter), scaled_age)
    #plot_data.sort(key = lambda x: x[3], reverse = True)
    random.shuffle(plot_data)

    if color:
        palette = rainbow_palette
    else:
        palette = gray_palette

    g.plot(graph.data.points(plot_data,
                             x=1, y=2, size=3, color=4),
           [colorpoints(gradient=palette,
                        symbolattrs=[deco.filled],
                        size = 0.05 * unit.v_mm)])

    write_time(g, env_state)

    g.dolayout()
    for axisname in ["x", "y"]:
        for t in g.axes[axisname].data.ticks:
            if t.ticklevel is not None:
                g.stroke(g.axes[axisname].positioner.vgridpath(t.temp_v),
                         [style.linestyle.dotted])
    g.dodata()
    g.doaxes()

    add_canvas_color_bar(
        g,
        min = 0.0,
        max = max_age,
        xpos = g.xpos + g.width + grid_h_space,
        ybottom = g.ypos,
        ytop = g.ypos + g.height,
        title = r"age (hours)",
        palette = palette)

    if color:
        out_filename = "%s_color.pdf" % out_prefix
    else:
        out_filename = "%s_bw.pdf" % out_prefix
    g.writePDFfile(out_filename)
    if not color:
        print "figure height = %.1f cm" % unit.tocm(g.bbox().height())
        print "figure width = %.1f cm" % unit.tocm(g.bbox().width())
