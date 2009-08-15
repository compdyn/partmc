#!/usr/bin/env python
# Copyright (C) 2007-2009 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys
import copy as module_copy
from Scientific.IO.NetCDF import *
from pyx import *
sys.path.append("../tool")
from pmc_data_nc import *
from pmc_pyx import *
sys.path.append(".")
from fig_helper import *
from pyx import canvas, color, attr, text, style, unit, box, path
from pyx import trafo as trafomodule
from pyx.graph.axis import tick

out_prefix = "figs_aging/aging_aero_tau"

plot_info = {
    "num_low": {"label": r"$\tau_{\rm N}$",
                 "label_time": 12, "label_pos": [0, 0],
                 "label_offset": 1 * unit.v_mm,
                 "linewidth": style.linewidth.Thick,
                "grey_level": grey_level_dark,
                 "color": color_list[0], "pattern": line_style_list[0],
                 "graph": "g11"},
    "num_low_cond": {"label": r"$\tau^{\rm cond}_{\rm N}$",
                 "label_time": 14, "label_pos": [1, 1],
                 "label_offset": 1 * unit.v_mm,
                 "linewidth": style.linewidth.Thick,
                "grey_level": grey_level_light,
                 "color": color_list[2], "pattern": line_style_list[1],
                 "graph": "g11"},
    "num_mid": {"label": r"$\tau_{\rm N}$",
                 "label_time": 12, "label_pos": [0, 0],
                 "label_offset": 1.9 * unit.v_mm,
                 "linewidth": style.linewidth.Thick,
                "grey_level": grey_level_dark,
                 "color": color_list[0], "pattern": line_style_list[0],
                 "graph": "g21"},
    "num_mid_cond": {"label": r"$\tau^{\rm cond}_{\rm N}$",
                 "label_time": 18.3, "label_pos": [1, 1],
                 "label_offset": 7 * unit.v_mm,
                 "linewidth": style.linewidth.Thick,
                "grey_level": grey_level_light,
                 "color": color_list[2], "pattern": line_style_list[1],
                 "graph": "g21"},
    "num_high": {"label": r"$\tau_{\rm N}$",
                 "label_time": 12, "label_pos": [0, 0],
                 "label_offset": 2 * unit.v_mm,
                 "linewidth": style.linewidth.Thick,
                "grey_level": grey_level_dark,
                 "color": color_list[0], "pattern": line_style_list[0],
                 "graph": "g31"},
    "num_high_cond": {"label": r"$\tau^{\rm cond}_{\rm N}$",
                 "label_time": 16, "label_pos": [1, 1],
                 "label_offset": 5 * unit.v_mm,
                 "linewidth": style.linewidth.Thick,
                "grey_level": grey_level_light,
                 "color": color_list[2], "pattern": line_style_list[1],
                 "graph": "g31"},
    "mass_low": {"label": r"$\tau_{\rm M}$",
                 "label_time": 11, "label_pos": [0, 0],
                 "label_offset": 2.5 * unit.v_mm,
                 "linewidth": style.linewidth.Thick,
                "grey_level": grey_level_dark,
                 "color": color_list[0], "pattern": line_style_list[0],
                 "graph": "g12"},
    "mass_low_cond": {"label": r"$\tau^{\rm cond}_{\rm M}$",
                 "label_time": 13, "label_pos": [1, 1],
                 "label_offset": 1 * unit.v_mm,
                 "linewidth": style.linewidth.Thick,
                "grey_level": grey_level_light,
                 "color": color_list[2], "pattern": line_style_list[1],
                 "graph": "g12"},
    "mass_mid": {"label": r"$\tau_{\rm M}$",
                 "label_time": 12, "label_pos": [0, 0],
                 "label_offset": 1.5 * unit.v_mm,
                 "linewidth": style.linewidth.Thick,
                "grey_level": grey_level_dark,
                 "color": color_list[0], "pattern": line_style_list[0],
                 "graph": "g22"},
    "mass_mid_cond": {"label": r"$\tau^{\rm cond}_{\rm M}$",
                 "label_time": 14.5, "label_pos": [1, 1],
                 "label_offset": 4 * unit.v_mm,
                 "linewidth": style.linewidth.Thick,
                "grey_level": grey_level_light,
                 "color": color_list[2], "pattern": line_style_list[1],
                 "graph": "g22"},
    "mass_high": {"label": r"$\tau_{\rm M}$",
                 "label_time": 12, "label_pos": [0, 0],
                 "label_offset": 1 * unit.v_mm,
                 "linewidth": style.linewidth.Thick,
                "grey_level": grey_level_dark,
                 "color": color_list[0], "pattern": line_style_list[0],
                 "graph": "g32"},
    "mass_high_cond": {"label": r"$\tau^{\rm cond}_{\rm M}$",
                 "label_time": 13.1, "label_pos": [1, 1],
                 "label_offset": 2 * unit.v_mm,
                 "linewidth": style.linewidth.Thick,
                "grey_level": grey_level_light,
                 "color": color_list[2], "pattern": line_style_list[1],
                 "graph": "g32"},
    }

time = loadtxt("%s/aging_wc_num_time.txt" % aging_data_dir)
comp_vol = loadtxt("%s/aging_wc_num_comp_vol.txt" % aging_data_dir)

num_tau_transfer = loadtxt("%s/aging_wc_num_tau_transfer.txt" % aging_data_dir)
num_tau_transfer_cond = loadtxt("%s/aging_wc_num_tau_transfer_cond.txt" % aging_data_dir)
mass_tau_transfer = loadtxt("%s/aging_wc_mass_tau_transfer.txt" % aging_data_dir)
mass_tau_transfer_cond = loadtxt("%s/aging_wc_mass_tau_transfer_cond.txt" % aging_data_dir)

num_tau_transfer_smooth = loadtxt("%s/aging_wc_num_tau_transfer_smooth.txt" % aging_data_dir)
num_tau_transfer_cond_smooth = loadtxt("%s/aging_wc_num_tau_transfer_cond_smooth.txt" % aging_data_dir)
mass_tau_transfer_smooth = loadtxt("%s/aging_wc_mass_tau_transfer_smooth.txt" % aging_data_dir)
mass_tau_transfer_cond_smooth = loadtxt("%s/aging_wc_mass_tau_transfer_cond_smooth.txt" % aging_data_dir)

print "low level %d = %f%%" % (level_low, level_low_value * 100)
print "mid level %d = %f%%" % (level_mid, level_mid_value * 100)
print "high level %d = %f%%" % (level_high, level_high_value * 100)

print ss_active_axis.edge(level_low-1), ss_active_axis.edge(level_low), ss_active_axis.edge(level_low+1)
print ss_active_axis.edge(level_mid-1), ss_active_axis.edge(level_mid), ss_active_axis.edge(level_mid+1)
print ss_active_axis.edge(level_high-1), ss_active_axis.edge(level_high), ss_active_axis.edge(level_high+1)


print ss_active_axis.edge(level_low-1) / ss_active_axis.edge(level_low), ss_active_axis.edge(level_low) / ss_active_axis.edge(level_low+1)
print ss_active_axis.edge(level_mid-1) / ss_active_axis.edge(level_mid), ss_active_axis.edge(level_mid) / ss_active_axis.edge(level_mid+1)
print ss_active_axis.edge(level_high-1) / ss_active_axis.edge(level_high), ss_active_axis.edge(level_high) / ss_active_axis.edge(level_high+1)

num_low = num_tau_transfer[:,level_low] / 3600 # s to hour
num_mid = num_tau_transfer[:,level_mid] / 3600 # s to hour
num_high = num_tau_transfer[:,level_high] / 3600 # s to hour
num_low_cond = num_tau_transfer_cond[:,level_low] / 3600 # s to hour
num_mid_cond = num_tau_transfer_cond[:,level_mid] / 3600 # s to hour
num_high_cond = num_tau_transfer_cond[:,level_high] / 3600 # s to hour
mass_low = mass_tau_transfer[:,level_low] / 3600 # s to hour
mass_mid = mass_tau_transfer[:,level_mid] / 3600 # s to hour
mass_high = mass_tau_transfer[:,level_high] / 3600 # s to hour
mass_low_cond = mass_tau_transfer_cond[:,level_low] / 3600 # s to hour
mass_mid_cond = mass_tau_transfer_cond[:,level_mid] / 3600 # s to hour
mass_high_cond = mass_tau_transfer_cond[:,level_high] / 3600 # s to hour

num_low_smooth = num_tau_transfer_smooth[:,level_low] / 3600 # s to hour
num_mid_smooth = num_tau_transfer_smooth[:,level_mid] / 3600 # s to hour
num_high_smooth = num_tau_transfer_smooth[:,level_high] / 3600 # s to hour
num_low_cond_smooth = num_tau_transfer_cond_smooth[:,level_low] / 3600 # s to hour
num_mid_cond_smooth = num_tau_transfer_cond_smooth[:,level_mid] / 3600 # s to hour
num_high_cond_smooth = num_tau_transfer_cond_smooth[:,level_high] / 3600 # s to hour
mass_low_smooth = mass_tau_transfer_smooth[:,level_low] / 3600 # s to hour
mass_mid_smooth = mass_tau_transfer_smooth[:,level_mid] / 3600 # s to hour
mass_high_smooth = mass_tau_transfer_smooth[:,level_high] / 3600 # s to hour
mass_low_cond_smooth = mass_tau_transfer_cond_smooth[:,level_low] / 3600 # s to hour
mass_mid_cond_smooth = mass_tau_transfer_cond_smooth[:,level_mid] / 3600 # s to hour
mass_high_cond_smooth = mass_tau_transfer_cond_smooth[:,level_high] / 3600 # s to hour

env_state = read_any(env_state_t, netcdf_dir_wc, netcdf_pattern_wc)
start_time_of_day_min = env_state.start_time_of_day / 60
max_time_min = max(time) / 60

class regular_offset(graph.axis.painter._title):
    """class for painting the ticks and labels of an axis"""

    defaulttickattrs = []
    defaultgridattrs = []
    defaultbasepathattrs = [style.linecap.square]
    defaultlabelattrs = [text.halign.center, text.vshift.mathaxis]

    def __init__(self, innerticklength=graph.axis.painter.ticklength.normal,
                       outerticklength=None,
                       tickattrs=[],
                       gridattrs=None,
                       basepathattrs=[],
                       labeldist=0.3*unit.v_cm,
                       labelattrs=[],
                       labeldirection=None,
                       labelhequalize=0,
                       labelvequalize=1,
                       first_offset = None,
                       last_offset = None,
                       **kwargs):
        self.innerticklength = innerticklength
        self.outerticklength = outerticklength
        self.tickattrs = tickattrs
        self.gridattrs = gridattrs
        self.basepathattrs = basepathattrs
        self.labeldist = labeldist
        self.labelattrs = labelattrs
        self.labeldirection = labeldirection
        self.labelhequalize = labelhequalize
        self.labelvequalize = labelvequalize
        self.first_offset = first_offset
        self.last_offset = last_offset
        graph.axis.painter._title.__init__(self, **kwargs)

    def paint(self, canvas, data, axis, axispos):
        for t in data.ticks:
            t.temp_v = axis.convert(data, t)
            t.temp_x_pt, t.temp_y_pt = axispos.vtickpoint_pt(t.temp_v)
            t.temp_dx, t.temp_dy = axispos.vtickdirection(t.temp_v)
        maxticklevel, maxlabellevel = tick.maxlevels(data.ticks)
        labeldist_pt = unit.topt(self.labeldist)

        # create & align t.temp_labelbox
        for (tick_index,t) in enumerate(data.ticks):
            if t.labellevel is not None:
                labelattrs = attr.selectattrs(self.labelattrs, t.labellevel, maxlabellevel)
                if labelattrs is not None:
                    labelattrs = self.defaultlabelattrs + labelattrs
                    if self.labeldirection is not None:
                        labelattrs.append(self.labeldirection.trafo(t.temp_dx, t.temp_dy))
                    if t.labelattrs is not None:
                        labelattrs.extend(t.labelattrs)
                    if (tick_index == 0) and self.first_offset:
                        t.temp_labelbox = canvas.texrunner.text_pt(t.temp_x_pt + self.first_offset, t.temp_y_pt, t.label, labelattrs)
                    if (tick_index == len(data.ticks) - 1) and self.last_offset:
                        t.temp_labelbox = canvas.texrunner.text_pt(t.temp_x_pt + self.last_offset, t.temp_y_pt, t.label, labelattrs)
                    t.temp_labelbox = canvas.texrunner.text_pt(t.temp_x_pt, t.temp_y_pt, t.label, labelattrs)
        if len(data.ticks) > 1:
            equaldirection = 1
            for t in data.ticks[1:]:
                if t.temp_dx != data.ticks[0].temp_dx or t.temp_dy != data.ticks[0].temp_dy:
                    equaldirection = 0
        else:
            equaldirection = 0
        if equaldirection and ((not data.ticks[0].temp_dx and self.labelvequalize) or
                               (not data.ticks[0].temp_dy and self.labelhequalize)):
            if self.labelattrs is not None:
                box.linealignequal_pt([t.temp_labelbox for t in data.ticks if t.labellevel is not None],
                                      labeldist_pt, -data.ticks[0].temp_dx, -data.ticks[0].temp_dy)
        else:
            for t in data.ticks:
                if t.labellevel is not None and self.labelattrs is not None:
                    t.temp_labelbox.linealign_pt(labeldist_pt, -t.temp_dx, -t.temp_dy)

        for t in data.ticks:
            if t.ticklevel is not None and self.tickattrs is not None:
                tickattrs = attr.selectattrs(self.defaulttickattrs + self.tickattrs, t.ticklevel, maxticklevel)
                if tickattrs is not None:
                    innerticklength = attr.selectattr(self.innerticklength, t.ticklevel, maxticklevel)
                    outerticklength = attr.selectattr(self.outerticklength, t.ticklevel, maxticklevel)
                    if innerticklength is not None or outerticklength is not None:
                        if innerticklength is None:
                            innerticklength = 0
                        if outerticklength is None:
                            outerticklength = 0
                        innerticklength_pt = unit.topt(innerticklength)
                        outerticklength_pt = unit.topt(outerticklength)
                        x1 = t.temp_x_pt + t.temp_dx * innerticklength_pt
                        y1 = t.temp_y_pt + t.temp_dy * innerticklength_pt
                        x2 = t.temp_x_pt - t.temp_dx * outerticklength_pt
                        y2 = t.temp_y_pt - t.temp_dy * outerticklength_pt
                        canvas.stroke(path.line_pt(x1, y1, x2, y2), tickattrs)
                        if outerticklength_pt > canvas.extent_pt:
                            canvas.extent_pt = outerticklength_pt
                        if -innerticklength_pt > canvas.extent_pt:
                            canvas.extent_pt = -innerticklength_pt
            if self.gridattrs is not None:
                gridattrs = attr.selectattrs(self.defaultgridattrs + self.gridattrs, t.ticklevel, maxticklevel)
                if gridattrs is not None:
                    canvas.stroke(axispos.vgridpath(t.temp_v), gridattrs)
            if t.labellevel is not None and self.labelattrs is not None:
                canvas.insert(t.temp_labelbox)
                canvas.labels.append(t.temp_labelbox)
                extent_pt = t.temp_labelbox.extent_pt(t.temp_dx, t.temp_dy) + labeldist_pt
                if extent_pt > canvas.extent_pt:
                    canvas.extent_pt = extent_pt

        if self.labelattrs is None:
            canvas.labels = None

        if self.basepathattrs is not None:
            canvas.stroke(axispos.vbasepath(), self.defaultbasepathattrs + self.basepathattrs)

        # for t in data.ticks:
        #     del t.temp_v    # we've inserted those temporary variables ... and do not care any longer about them
        #     del t.temp_x_pt
        #     del t.temp_y_pt
        #     del t.temp_dx
        #     del t.temp_dy
        #     if t.labellevel is not None and self.labelattrs is not None:
        #         del t.temp_labelbox

        graph.axis.painter._title.paint(self, canvas, data, axis, axispos)

for use_color in [True, False]:
    c = canvas.canvas()
    
    x_axis_left = graph.axis.linear(
        min = 0.,
        max = max_time_min,
        parter = graph.axis.parter.linear(tickdists
                                          = [6 * 60, 3 * 60]),
        texter = time_of_day(base_time
                             = start_time_of_day_min,
                             last_suffix = r"\ \ \ \ \ \ "),
        title = "local standard time (LST) (hours:minutes)")

    x_axis_right = graph.axis.linear(
        min = 0.,
        max = max_time_min,
        parter = graph.axis.parter.linear(tickdists
                                          = [6 * 60, 3 * 60]),
        texter = time_of_day(base_time
                             = start_time_of_day_min,
                             first_prefix = r"\ \ \ \ \ \ "),
        title = "local standard time (LST) (hours:minutes)")

    y_axis = graph.axis.log(
        min = 1e-1,
        max = 1e4,
        title = r"aging timescale $\tau$ (hours)")

    g31 = c.insert(graph.graphxy(
        width = 6.8,
        x = x_axis_left,
        y = y_axis))
    g21 = c.insert(graph.graphxy(
        width = 6.8,
        ypos = g31.ypos + g31.height + grid_v_space,
        x = graph.axis.linkedaxis(g31.axes["x"],
                                  painter = graph.axis.painter.linked()),
        y = y_axis))
    g11 = c.insert(graph.graphxy(
        width = 6.8,
        ypos = g21.ypos + g21.height + grid_v_space,
        x = graph.axis.linkedaxis(g31.axes["x"],
                                  painter = graph.axis.painter.linked()),
        y = y_axis))
    g32 = c.insert(graph.graphxy(
        width = 6.8,
        xpos = g31.xpos + g31.width + grid_h_space,
        x = x_axis_right,
        y = graph.axis.linkedaxis(g31.axes["y"],
                                  painter = graph.axis.painter.linked())))
    g22 = c.insert(graph.graphxy(
        width = 6.8,
        xpos = g32.xpos,
        ypos = g21.ypos,
        x = graph.axis.linkedaxis(g32.axes["x"],
                                  painter = graph.axis.painter.linked()),
        y = graph.axis.linkedaxis(g21.axes["y"],
                                  painter = graph.axis.painter.linked())))
    g12 = c.insert(graph.graphxy(
        width = 6.8,
        xpos = g32.xpos,
        ypos = g11.ypos,
        x = graph.axis.linkedaxis(g32.axes["x"],
                                  painter = graph.axis.painter.linked()),
        y = graph.axis.linkedaxis(g11.axes["y"],
                                  painter = graph.axis.painter.linked())))

    graphs = {"g11": g11, "g21": g21, "g31": g31,
              "g12": g12, "g22": g22, "g32": g32}

    for (key, y_data) \
            in [("num_low_cond", num_low_cond),
                ("num_low", num_low),
                ("num_mid_cond", num_mid_cond),
                ("num_mid", num_mid),
                ("num_high_cond", num_high_cond),
                ("num_high", num_high),
                ("mass_low_cond", mass_low_cond),
                ("mass_low", mass_low),
                ("mass_mid_cond", mass_mid_cond),
                ("mass_mid", mass_mid),
                ("mass_high_cond", mass_high_cond),
                ("mass_high", mass_high)]:
        g = graphs[plot_info[key]["graph"]]
        if use_color:
            grey_color = color.hsb(plot_info[key]["color"].hsb().color["h"], plot_info[key]["grey_level"], 1)
            style_attrs = [plot_info[key]["linewidth"],
                           grey_color]
        else:
            grey_color = color.grey(1 - plot_info[key]["grey_level"])
            style_attrs = [plot_info[key]["linewidth"],
                           grey_color]
        plot_data = zip(time[1:] / 60, y_data)
        g.plot(
            graph.data.points(plot_data, x = 1, y = 2),
            styles = [graph.style.line(lineattrs = style_attrs)])

    for (key, y_data) \
            in [("num_low_cond", num_low_cond_smooth),
                ("num_low", num_low_smooth),
                ("num_mid_cond", num_mid_cond_smooth),
                ("num_mid", num_mid_smooth),
                ("num_high_cond", num_high_cond_smooth),
                ("num_high", num_high_smooth),
                ("mass_low_cond", mass_low_cond_smooth),
                ("mass_low", mass_low_smooth),
                ("mass_mid_cond", mass_mid_cond_smooth),
                ("mass_mid", mass_mid_smooth),
                ("mass_high_cond", mass_high_cond_smooth),
                ("mass_high", mass_high_smooth)]:
        g = graphs[plot_info[key]["graph"]]
        if use_color:
            style_attrs = [plot_info[key]["linewidth"],
                           plot_info[key]["color"]]
        else:
            style_attrs = [plot_info[key]["linewidth"],
                           plot_info[key]["pattern"]]
        plot_data = zip(time[1:] / 60, y_data)
        g.plot(
            graph.data.points(plot_data, x = 1, y = 2),
            styles = [graph.style.line(lineattrs = style_attrs)])

    for (g_name, g) in graphs.iteritems():
        g.dolayout()
        for axisname in ["x"]:
            for t in g.axes[axisname].data.ticks:
                if t.ticklevel is not None:
                    g.stroke(g.axes[axisname].positioner.vgridpath(t.temp_v),
                             [style.linestyle.dotted])
        for axisname in ["y"]:
            for t in g.axes[axisname].data.ticks:
                if t.ticklevel == 0:
                    g.stroke(g.axes[axisname].positioner.vgridpath(t.temp_v),
                             [style.linestyle.dotted])

    for (key, y_data) \
            in [("num_low_cond", num_low_cond_smooth),
                ("num_low", num_low_smooth),
                ("num_mid_cond", num_mid_cond_smooth),
                ("num_mid", num_mid_smooth),
                ("num_high_cond", num_high_cond_smooth),
                ("num_high", num_high_smooth),
                ("mass_low_cond", mass_low_cond_smooth),
                ("mass_low", mass_low_smooth),
                ("mass_mid_cond", mass_mid_cond_smooth),
                ("mass_mid", mass_mid_smooth),
                ("mass_high_cond", mass_high_cond_smooth),
                ("mass_high", mass_high_smooth)]:
        g = graphs[plot_info[key]["graph"]]
        plot_data = zip(time[1:] / 60, y_data)
        label_plot_line_boxed(g, plot_data,
                              plot_info[key]["label_time"] * 60,
                              plot_info[key]["label"],
                              plot_info[key]["label_pos"],
                              draw_text = False,
                              label_offset = plot_info[key]["label_offset"])

    for (g_name, g) in graphs.iteritems():
        g.dodata()
        g.doaxes()

    for (key, y_data) \
            in [("num_low_cond", num_low_cond_smooth),
                ("num_low", num_low_smooth),
                ("num_mid_cond", num_mid_cond_smooth),
                ("num_mid", num_mid_smooth),
                ("num_high_cond", num_high_cond_smooth),
                ("num_high", num_high_smooth),
                ("mass_low_cond", mass_low_cond_smooth),
                ("mass_low", mass_low_smooth),
                ("mass_mid_cond", mass_mid_cond_smooth),
                ("mass_mid", mass_mid_smooth),
                ("mass_high_cond", mass_high_cond_smooth),
                ("mass_high", mass_high_smooth)]:
        g = graphs[plot_info[key]["graph"]]
        plot_data = zip(time[1:] / 60, y_data)
        label_plot_line_boxed(g, plot_data,
                              plot_info[key]["label_time"] * 60,
                              plot_info[key]["label"],
                              plot_info[key]["label_pos"],
                              draw_box = False,
                              label_offset = plot_info[key]["label_offset"])

    write_text_outside(g11, r"critical supersaturation $S_{\rm c} = %.1f\%%$" % (level_low_value * 100))
    write_text_outside(g21, r"critical supersaturation $S_{\rm c} = %.1f\%%$" % (level_mid_value * 100))
    write_text_outside(g31, r"critical supersaturation $S_{\rm c} = %.1f\%%$" % (level_high_value * 100))
    write_text_outside(g12, r"critical supersaturation $S_{\rm c} = %.1f\%%$" % (level_low_value * 100))
    write_text_outside(g22, r"critical supersaturation $S_{\rm c} = %.1f\%%$" % (level_mid_value * 100))
    write_text_outside(g32, r"critical supersaturation $S_{\rm c} = %.1f\%%$" % (level_high_value * 100))

    boxed_text(g11, "number", point = [1, 1], anchor_point_rel = [1, 1])
    boxed_text(g21, "number", point = [1, 1], anchor_point_rel = [1, 1])
    boxed_text(g31, "number", point = [1, 1], anchor_point_rel = [1, 1])
    boxed_text(g12, "mass", point = [1, 1], anchor_point_rel = [1, 1])
    boxed_text(g22, "mass", point = [1, 1], anchor_point_rel = [1, 1])
    boxed_text(g32, "mass", point = [1, 1], anchor_point_rel = [1, 1])

    if use_color:
        out_filename = "%s_color.pdf" % out_prefix
    else:
        out_filename = "%s_bw.pdf" % out_prefix
    c.writePDFfile(out_filename)
    if not use_color:
        print "figure height = %.1f cm" % unit.tocm(c.bbox().height())
        print "figure width = %.1f cm" % unit.tocm(c.bbox().width())
