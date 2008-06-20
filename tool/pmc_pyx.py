#!/usr/bin/env python
# Copyright (C) 2007-2008 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys
sys.path.append(os.path.expanduser("~/.python"))
from pyx import *

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

rainbow_palette = listpalette([[0, color.rgb(0, 0, 1)],
			       [0.3, color.rgb(0, 1, 1)],
			       [0.5, color.rgb(0, 1, 0)],
			       [0.7, color.rgb(1, 1, 0)],
			       [1, color.rgb(1, 0, 0)]])

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
	    [graph.style.rect(palette)])
    gc.dolayout()
    gc.dobackground()
    gc.dodata()
    gc.doaxes()

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
		  bar_height_ratio = 0.6, bar_x_offset = -1.2, bar_y_offset = -0.7):
    colorbar_steps = 1000
    color_d = []
    for i in range(colorbar_steps):
	x0 = float(i) / float(colorbar_steps)
	xh = (float(i) + 0.5) / float(colorbar_steps)
	x1 = float(i + 1) / float(colorbar_steps)
	v0 = x0 * (max - min) + min
	v1 = x1 * (max - min) + min
	color_d.append([0, 1, v0, v1, xh])
    gc = c.insert(
	graph.graphxy(
            width = bar_width,
	    height = bar_height_ratio * c.bbox().height(),
	    xpos = c.bbox().width() + bar_x_offset,
	    ypos = (1.0 - bar_height_ratio) / 2.0 * c.bbox().height() \
            + bar_y_offset,
	    x = graph.axis.linear(min = 0, max = 1,
				  parter = None),
	    y2 = graph.axis.linear(
		min = min,
		max = max,
		title = title)))
    gc.plot(graph.data.points(color_d, xmin = 1, xmax = 2,
                            ymin = 3, ymax = 4, color = 5),
	    [graph.style.rect(palette)])
    gc.dolayout()
    gc.dobackground()
    gc.dodata()
    gc.doaxes()

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

