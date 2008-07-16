#!/usr/bin/env python
# Copyright (C) 2007-2008 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys
sys.path.append(os.path.expanduser("~/.python"))
from pyx import *
import pyx.bbox as bbox
import numpy
import pmc_data_nc

text.set(mode="latex",usefiles=["spam.aux"],texdebug="spam.debug")
#text.set(mode="latex")
#text.set(docopt="10pt")
#text.set(fontmaps="download35.map")
#text.preamble(r"\usepackage{times}")

#text.preamble(r"""\usepackage{times}
#\usepackage{sfmath}
#\renewcommand{\familydefault}{\sfdefault}
#\renewcommand{\normalsize}{\fontsize{8}{8}\selectfont}""")

text.preamble(r"""\renewcommand{\sfdefault}{phv}
\renewcommand{\familydefault}{\sfdefault}
\renewcommand{\normalsize}{\fontsize{9}{9}\selectfont}
\usepackage{sfmath}""")

#text.preamble(r"""\usepackage{amsmath}
#\renewcommand{\normalsize}{\fontfamily{phvv}\selectfont}""")

#text.preamble(r"""\usepackage{amsmath}
#\usepackage{times}
#\usepackage{sfmath}
#\renewcommand{\familydefault}{\sfdefault}
#\renewcommand{\normalsize}{\fontfamily{phvv}\fontsize{3}{3}\selectfont}""")

color_list = [color.hsb(2/3.0, 1, 1),
	      color.hsb(1/3.0, 1, 1),
	      color.hsb(0/3.0, 1, 1),
	      color.hsb(5/6.0, 1, 1),
	      color.hsb(3/6.0, 1, 1),
	      color.hsb(1/6.0, 1, 1),
	      color.hsb(11/12.0, 1, 1),
	      color.hsb(9/12.0, 1, 1),
	      color.hsb(7/12.0, 1, 1),
	      color.hsb(5/12.0, 1, 1),
	      color.hsb(3/12.0, 1, 1),
	      color.hsb(1/12.0, 1, 1),
	      ]

line_style_list = [style.linestyle.solid,
                   style.linestyle(style.linecap.round, # dash
                                   style.dash([2, 2])),
                   style.linestyle(style.linecap.round, # dot
                                   style.dash([0, 2])),
                   style.linestyle(style.linecap.round, # dot-dash
                                   style.dash([0, 2, 2, 2])),
                   style.linestyle(style.linecap.round, # dot-dot-dash
                                   style.dash([0, 2, 0, 2, 2, 2])),
                   style.linestyle(style.linecap.round, # dot-dash-dash
                                   style.dash([0, 2, 2, 2, 2, 2])),
                   style.linestyle(style.linecap.round, # dot-dot-dot-dash
                                   style.dash([0, 2, 0, 2, 0, 2, 2, 2])),
                   style.linestyle(style.linecap.round, # dot-dot-dash-dash
                                   style.dash([0, 2, 0, 2, 2, 2, 2, 2])),
                   style.linestyle(style.linecap.round, # dot-dash-dash-dash
                                   style.dash([0, 2, 2, 2, 2, 2, 2, 2])),
                   ]

class listpalette(color.palette):

    def __init__(self, colorlist):
        self.colorclass = colorlist[0][1].__class__
        self.colorlist = colorlist

    def getcolor(self, param):
	for i in range(len(self.colorlist)):
	    if self.colorlist[i][0] >= param:
		break
	else:
	    raise ValueError
	if i == 0:
	    i = 1
	# list[i-1] < param < list[i]
	alpha = (param - self.colorlist[i-1][0]) \
	    / (self.colorlist[i][0] - self.colorlist[i-1][0])
        colordict = {}
        for key in self.colorlist[0][1].color.keys():
            colordict[key] = alpha * self.colorlist[i][1].color[key] \
		+ (1 - alpha) * self.colorlist[i-1][1].color[key]
        return self.colorclass(**colordict)

rainbow_palette = listpalette([[0.0, color.rgb(0, 0, 1)],  # blue
			       [0.2, color.rgb(0, 1, 1)],  # cyan
			       [0.4, color.rgb(0, 1, 0)],  # green
			       [0.6, color.rgb(1, 1, 0)],  # yellow
			       [0.8, color.rgb(1, 0, 0)],  # red
                               [1.0, color.rgb(1, 0, 1)]]) # magenta

gray_palette = listpalette([[0, color.gray(0.8)],
                            [1, color.gray(0)]])

grid_painter = graph.axis.painter.regular(gridattrs = [style.linestyle.dotted])
major_grid_painter = graph.axis.painter.regular(gridattrs = [attr.changelist([style.linestyle.dotted, None])])

aerosol_species_tex = {
    "SO4": "SO$_4$",
    "NO3": "NO$_3$",
    "Cl": "Cl",
    "NH4": "NH$_4$",
    "MSA": "MSA",
    "ARO1": "ARO1",
    "ARO2": "ARO2",
    "ALK1": "ALK1",
    "OLE1": "OLE1",
    "API1": "API1",
    "API2": "API2",
    "LIM1": "LIM1",
    "LIM2": "LIM2",
    "CO3": "CO$_3$",
    "Na": "Na",
    "Ca": "Ca",
    "OIN": "OIN",
    "OC": "POM",
    "BC": "soot",
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

def path_rounded_rect(x1, y1, x2, y2, r):
    return path.path(path.arc(x1 + r, y1 + r, r, 180, 270),
                     path.arc(x2 - r, y1 + r, r, 270, 360),
                     path.arc(x2 - r, y2 - r, r, 0, 90),
                     path.arc(x1 + r, y2 - r, r, 90, 180),
                     path.closepath())

def label_point(g, x, y, label_x, label_y, label, radius = 1 * unit.t_mm):
    (x_g, y_g) = g.pos(x, y)
    (label_x_g, label_y_g) = g.pos(label_x, label_y)
    c = canvas.canvas()
    c.text(0, 0, label, 
           [text.halign.boxcenter, text.halign.flushcenter,
                               text.valign.middle])
    width = c.bbox().width()
    height = c.bbox().height()
    x_off = width / 2.0 + radius
    y_off = height / 2.0 + radius
    g.stroke(path.line(label_x_g, label_y_g, x_g, y_g),
             [color.gray.white, style.linewidth.THIck])
    g.stroke(path.line(label_x_g, label_y_g, x_g, y_g))
    g.draw(path_rounded_rect(label_x_g - x_off, label_y_g - y_off,
                             label_x_g + x_off, label_y_g + y_off, radius),
           [deco.stroked, deco.filled([color.gray.white])])
    g.insert(c, [trafo.translate(label_x_g, label_y_g)])

def boxed_text(g, x, y, label, border = 1 * unit.v_mm):
    (x_g, y_g) = g.vpos(x, y)
    c = text.text(x_g, y_g, label)
    b = c.bbox().enlarged(all = border)
    g.draw(b.path(), [deco.stroked, deco.filled([color.gray.white])])
    g.insert(c)

def add_color_bar(g, min, max, title, palette, bar_width = 0.5,
		  bar_height_ratio = 0.8, bar_x_offset = 1.8,
                  texter = None):
    colorbar_steps = 1000
    color_d = []
    for i in range(colorbar_steps):
	x0 = float(i) / float(colorbar_steps)
	xh = (float(i) + 0.5) / float(colorbar_steps)
	x1 = float(i + 1) / float(colorbar_steps)
	v0 = x0 * (max - min) + min
	v1 = x1 * (max - min) + min
	color_d.append([0, 1, v0, v1, xh])
    if texter:
        y2axis = graph.axis.linear(
            min = min,
            max = max,
            title = title,
            texter = texter)
    else:
        y2axis = graph.axis.linear(
            min = min,
            max = max,
            title = title)
    gc = g.insert(
	graph.graphxy(
	    width = bar_width,
	    height = bar_height_ratio * g.height,
	    xpos = g.width + bar_x_offset,
	    ypos = (1.0 - bar_height_ratio) / 2.0 * g.height,
	    x = graph.axis.linear(min = 0, max = 1,
				  parter = None),
	    y2 = y2axis))
    gc.plot(graph.data.points(color_d, xmin = 1, xmax = 2,
                            ymin = 3, ymax = 4, color = 5),
	    [hsb_rect(palette)])
    gc.dolayout()
    gc.dobackground()
    gc.dodata()
    gc.doaxes()

def add_color_bar_new(g, min, max, title, bar_width = 0.5,
                      bar_height_ratio = 0.8, bar_x_offset = 1.8,
                      texter = None):
    colorbar_steps = 1000
    value = linspace(0, 1, colorbar_steps).reshape([1,colorbar_steps])
    hue = value_to_hue(value)
    saturation = ones(shape(hue))
    brightness = ones(shape(hue))
    if texter:
        y2axis = graph.axis.linear(
            min = min,
            max = max,
            title = title,
            texter = texter)
    else:
        y2axis = graph.axis.linear(
            min = min,
            max = max,
            title = title)
    gc = g.insert(
	graph.graphxy(
	    width = bar_width,
	    height = bar_height_ratio * g.height,
	    xpos = g.width + bar_x_offset,
	    ypos = (1.0 - bar_height_ratio) / 2.0 * g.height,
	    x = graph.axis.linear(min = 0, max = 1,
				  parter = None),
	    y2 = y2axis))
    pmc_plot_image(gc, hue, saturation, brightness)

def add_horiz_color_bar(g, min, max, title, palette, bar_height = 0.5,
		  bar_width_ratio = 0.8, bar_offset = 1.8):
    colorbar_steps = 1000
    color_d = []
    for i in range(colorbar_steps):
	x0 = float(i) / float(colorbar_steps)
	xh = (float(i) + 0.5) / float(colorbar_steps)
	x1 = float(i + 1) / float(colorbar_steps)
	v0 = x0 * (max - min) + min
	v1 = x1 * (max - min) + min
	color_d.append([0, 1, v0, v1, xh])
    gc = g.insert(
	graph.graphxy(
	    width = bar_width_ratio * g.width,
	    height = bar_height,
	    xpos = (1.0 - bar_width_ratio) / 2.0 * g.width,
	    ypos = g.height + bar_offset,
	    y = graph.axis.linear(min = 0, max = 1,
				  parter = None),
	    x2 = graph.axis.linear(
		min = min,
		max = max,
		title = title)))
    gc.plot(graph.data.points(color_d, ymin = 1, ymax = 2,
                            xmin = 3, xmax = 4, color = 5),
	    [graph.style.rect(palette)])
    gc.dolayout()
    gc.dobackground()
    gc.dodata()
    gc.doaxes()

def add_canvas_color_bar(c, min, max, title, palette, bar_width = 0.5,
		  bar_height_ratio = 0.6, bar_x_offset = -1.2,
                         bar_y_offset = -0.7, min_palette_index = 0.0,
                         max_palette_index = 1.0, texter = None,
                         extra_box_value = None, extra_box_label = None,
                         extra_box_pattern = None):
    colorbar_steps = 1000
    color_d = []
    for i in range(colorbar_steps):
	x0 = float(i) / float(colorbar_steps)
	xh = (float(i) + 0.5) / float(colorbar_steps)
	x1 = float(i + 1) / float(colorbar_steps)
	v0 = x0 * (max - min) + min
	v1 = x1 * (max - min) + min
        pi = (1 - xh) * min_palette_index + xh * max_palette_index
	color_d.append([0, 1, v0, v1, pi])
    if texter:
        y2axis = graph.axis.linear(
            min = min,
            max = max,
            title = title,
            texter = texter)
    else:
        y2axis = graph.axis.linear(
            min = min,
            max = max,
            title = title)
    xpos = c.bbox().width() + bar_x_offset
    ypos = (1.0 - bar_height_ratio) / 2.0 * c.bbox().height() + bar_y_offset
    gc = c.insert(
	graph.graphxy(
            width = bar_width,
	    height = bar_height_ratio * c.bbox().height(),
	    xpos = xpos,
	    ypos = ypos,
	    x = graph.axis.linear(min = 0, max = 1,
				  parter = None),
	    y2 = y2axis))
    gc.plot(graph.data.points(color_d, xmin = 1, xmax = 2,
                            ymin = 3, ymax = 4, color = 5),
	    [graph.style.rect(palette)])
    gc.dolayout()
    gc.dobackground()
    gc.dodata()
    gc.doaxes()
    if extra_box_value != None:
        box_fill_attributes = [palette.getcolor(extra_box_value)]
        if extra_box_pattern != None:
            box_fill_attributes.append(extra_box_pattern)
        c.draw(path.rect(xpos, ypos - 2 * bar_width, bar_width, bar_width),
               [deco.stroked([color.rgb.black]),
                deco.filled(box_fill_attributes)])
        if extra_box_label != None:
            c.text(xpos + bar_width + 0.3, ypos - 1.5 * bar_width,
                   extra_box_label, [text.halign.left, text.valign.middle])

class time_of_day:
    "a texter creating labels of the form 04:50 for 4 hours and 50 minutes"

    def __init__(self, base_time = 0, labelattrs = []):
        """initializes the instance
        - base_time is the offset to add (minutes)
        - labelattrs is a list of attributes to be added to the label
          attributes given in the painter"""
        self.base_time = base_time
        self.labelattrs = labelattrs

    def labels(self, ticks):
        for tick in ticks:
            if tick.label is None and tick.labellevel is not None:
                time = float(tick.num) / float(tick.denom) + self.base_time
                hours, minutes = divmod(time, 60)
                hours = hours % 24
                tick.label = "%02d:%02d" % (hours, minutes)
                tick.labelattrs = tick.labelattrs + self.labelattrs

class hsb_rect(graph.style._style):

    needsdata = ["vrange", "vrangeminmissing", "vrangemaxmissing"]

    def __init__(self, gradient=color.gradient.Grey,
                 fill_pattern = None, do_stroke = True):
        self.gradient = gradient
        self.fill_pattern = fill_pattern
        self.do_stroke = do_stroke

    def columnnames(self, privatedata, sharedata, graph, columnnames):
        if len(graph.axesnames) != 2:
            raise TypeError("arrow style restricted on two-dimensional graphs")
        if len(sharedata.vrangeminmissing) + len(sharedata.vrangemaxmissing):
            raise ValueError("incomplete range")
        ret_names = []
        for name in ["color", "hue", "saturation", "brightness",
                     "stroke_color", "stroke_width"]:
            if name in columnnames:
                ret_names.append(name)
        return ret_names

    def initdrawpoints(self, privatedata, sharedata, graph):
        privatedata.rectcanvas = graph.insert(canvas.canvas())

    def drawpoint(self, privatedata, sharedata, graph, point):
        xvmin = sharedata.vrange[0][0]
        xvmax = sharedata.vrange[0][1]
        yvmin = sharedata.vrange[1][0]
        yvmax = sharedata.vrange[1][1]
        if (xvmin is not None and xvmin < 1 and
            xvmax is not None and xvmax > 0 and
            yvmin is not None and yvmin < 1 and
            yvmax is not None and yvmax > 0):
            if xvmin < 0:
                xvmin = 0
            elif xvmax > 1:
                xvmax = 1
            if yvmin < 0:
                yvmin = 0
            elif yvmax > 1:
                yvmax = 1
            p = graph.vgeodesic(xvmin, yvmin, xvmax, yvmin)
            p.append(graph.vgeodesic_el(xvmax, yvmin, xvmax, yvmax))
            p.append(graph.vgeodesic_el(xvmax, yvmax, xvmin, yvmax))
            p.append(graph.vgeodesic_el(xvmin, yvmax, xvmin, yvmin))
            p.append(path.closepath())
            hue = 0.0
            saturation = 0.0
            brightness = 1.0
            if "color" in point.keys():
                c = self.gradient.getcolor(point["color"]).hsb()
                hue = c.color["h"]
                saturation = c.color["s"]
                brightness = c.color["b"]
            if "hue" in point.keys():
                hue = point["hue"]
            if "saturation" in point.keys():
                saturation = point["saturation"]
            if "brightness" in point.keys():
                brightness = point["brightness"]
            fill_color = color.hsb(hue, saturation, brightness)
            stroke_color = fill_color
            if "stroke_color" in point.keys():
                stroke_color = self.gradient.getcolor(point["stroke_color"])
            stroke_width = 0.0001
            if "stroke_width" in point.keys():
                stroke_width = point["stroke_width"]
            stroke_attributes = [stroke_color,
                                 style.linewidth(stroke_width * unit.v_cm)]
            fill_attributes = [fill_color]
            if self.fill_pattern != None:
                fill_attributes.append(self.fill_pattern)
            if self.do_stroke:
                privatedata.rectcanvas.draw(p, [deco.stroked(stroke_attributes),
                                                deco.filled(fill_attributes)])
            else:
                privatedata.rectcanvas.draw(p, [deco.filled(fill_attributes)])

    def key_pt(self, privatedata, sharedata, graph, x_pt, y_pt, width_pt, height_pt):
        raise RuntimeError("Style currently doesn't provide a graph key")

def pmc_plot_image(g, hue, saturation, brightness):
    if (numpy.shape(hue) != numpy.shape(saturation)) \
       or (numpy.shape(hue) != numpy.shape(brightness)):
        raise Exception("hue, saturation, and brightness must be the same shape")
    (nx, ny) = numpy.shape(hue)
    image_data = numpy.zeros([ny, nx, 3], dtype='uint8')
    for i in range(nx):
        for j in range(ny):
            pixel_hsb = color.hsb(hue[i,j], saturation[i,j], brightness[i,j])
            pixel_rgb = pixel_hsb.rgb()
            image_data[j,i,0] = int(pixel_rgb.color["r"] * 255)
            image_data[j,i,1] = int(pixel_rgb.color["g"] * 255)
            image_data[j,i,2] = int(pixel_rgb.color["b"] * 255)
    image_obj = bitmap.image(nx, ny, "RGB", image_data.tostring())
    bitmap_obj = bitmap.bitmap(0, 0, image_obj, height = 1)
    bitmap_width = float(nx) / float(ny)
    (x0, y0) = g.vpos(0, 0)
    (x1, y1) = g.vpos(1, 1)
    g.insert(bitmap_obj, [trafo.mirror(0), trafo.translate(0, 1),
                          trafo.scale((x1 - x0) / bitmap_width, y1 - y0),
                          trafo.translate(x0, y0)])

class pmc_gradient:

    def __init__(self, map_list):
        self.map_list = map_list

    def map(self, x):
        groups = numpy.digitize(x.ravel(),
                                bins = [val for [val, hue] in self.map_list])
        groups = groups.reshape(shape(x))
        x_arr = numpy.array([xa for [xa, ya] in self.map_list])
        y_arr = numpy.array([ya for [xa, ya] in self.map_list])
        x_arr = numpy.concatenate([array([x_arr[0] - 1.0]),
                                   x_arr, array([x_arr[-1] + 1.0])])
        y_arr = numpy.concatenate([array([y_arr[0]]),
                                   y_arr, array([y_arr[-1]])])
        x0 = x_arr[groups]
        x1 = x_arr[groups + 1]
        y0 = y_arr[groups]
        y1 = y_arr[groups + 1]
        y = (x - x0) / (x1 - x0) * (y1 - y0) + y0
        return y

rainbow_gradient = pmc_gradient([
    [0.0, color.rgb(0, 0, 1).hsb().color["h"]], # blue
    [0.2, color.rgb(0, 1, 1).hsb().color["h"]], # cyan
    [0.4, color.rgb(0, 1, 0).hsb().color["h"]], # green
    [0.6, color.rgb(1, 1, 0).hsb().color["h"]], # yellow
    [0.8, color.rgb(1, 0, 0).hsb().color["h"]], # red
    [1.0, color.rgb(1, 0, 1).hsb().color["h"]], # magenta
    ])

gray_gradient = pmc_gradient([[0.0, 0.8],
                              [1.0, 1.0]])

def value_to_hue(value):
    return (2.0/3.0 - value * 5.0 / 6.0) % 1.0

def hash_pattern(n_lines = 10, line_attributes = [style.linewidth.normal]):
    p = pattern.pattern(painttype = 2,
                        tilingtype = 2,
                        xstep = 1,
                        ystep = 1,
                        bbox = bbox.bbox(0, 0, 1, 1))
    p.stroke(path.line(1, 0, 0, 1), line_attributes)
    for i in range(n_lines):
        x = float(i) / float(n_lines)
        p.stroke(path.line(x, 0, 0, x), line_attributes)
        p.stroke(path.line(1 - x, 1, 1, 1 - x), line_attributes)
    return p

def label_plot_line(g, plot_data, label_time, label, label_pos = [0, 1],
                    label_offset = unit.v_mm, xaxis = None, yaxis = None):
    i = pmc_data_nc.find_nearest_time(plot_data, label_time)
    [label_x, label_y] = plot_data[i]
    [label_pos_h, label_pos_v] = label_pos
    [label_vx, label_vy] = g.pos(label_x, label_y, xaxis, yaxis)
    label_vx += 2.0 * (0.5 - label_pos_h) * label_offset
    label_vy -= 2.0 * (0.5 - label_pos_v) * label_offset
    g.text(label_vx, label_vy, label, [text.halign(label_pos_h, label_pos_h),
                                       text.valign(label_pos_v)])
