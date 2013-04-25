#!/usr/bin/env python

import sys, os
sys.path.append("../../tool")
import mpl_helper
import matplotlib
import matplotlib.pyplot as plt
import partmc, math
import scipy.io, numpy, numpy.ma

def make_fig_array(n_vert=2,
                   n_horiz=2,
                   figure_width=5,
                   axis_ratio=(1 + math.sqrt(5)) / 2, # golden ratio
                   left_margin=0.8,
                   right_margin=0.2,
                   bottom_margin=0.5,
                   top_margin=0.2,
                   horiz_sep=0.4,
                   vert_sep=0.4,
                   colorbar_location="left", # "left", "right", "top", "bottom"
                   colorbar=None, # "individual", "single", "shared"
                   top_colorbar=False,
                   colorbar_width=0.15,
                   colorbar_length_fraction=0.8,
                   colorbar_offset=0.2,
                   share_x_axes=True,
                   share_y_axes=True,
                   ):
    """
    (figure, axes_array) = make_fig_array()
    (figure, axes_array, cbar_axes_array) = make_fig_array(colorbar="individual")
    (figure, axes_array, cbar_axes_array) = make_fig_array(colorbar="shared")
    (figure, axes_array, cbar_axes) = make_fig_array(colorbar="single")

    Numbering convention:
    axes_array[i][j] has (i,j) layout:

          (2,0) (2,1)
          (1,0) (1,1)
          (0,0) (0,1)

    If colorbar is "individual" then cbar_axes_array has the same
    dimensions and format as axes_array.

    If colorbar is "shared" then cbar_axes_array is a 1D list of
    colorbars (one colorbar per row or per column).

    If colorbar is "single" then only one cbar_axes is returned.
    """
    axis_width = (figure_width - left_margin - right_margin - (n_horiz - 1) * horiz_sep) / float(n_horiz)
    axis_height = axis_width / axis_ratio
    figure_height = bottom_margin + axis_height * n_vert + vert_sep * (n_vert - 1) + top_margin
    figure = plt.figure()
    figure.set_figwidth(figure_width)
    figure.set_figheight(figure_height)
    axes_array = []
    for i_vert in range(n_vert):
        axes_array.append([])
        for i_horiz in range(n_horiz):
            x_left = left_margin + i_horiz * (axis_width + horiz_sep)
            y_bottom = bottom_margin + i_vert * (axis_height + vert_sep)
            kwargs = {}
            if i_horiz > 0 and share_y_axes:
                kwargs["sharey"] = last_y_axes
            if i_vert > 0 and share_x_axes:
                kwargs["sharex"] = last_x_axes
            new_axes = figure.add_axes([x_left / figure_width,
                                        y_bottom / figure_height,
                                        axis_width / figure_width,
                                        axis_height / figure_height],
                                       **kwargs)
            axes_array[-1].append(new_axes)
            if i_horiz == 0:
                last_y_axes = new_axes
            if i_vert == 0:
                last_x_axes = new_axes
    if colorbar in ["individual", "shared"]:
        cbar_axes_array = []
        for i_vert in range(n_vert):
            if colorbar == "individual":
                cbar_axes_array.append([])
            for i_horiz in range(1):
                if colorbar_location == "left":
                    x_left = left_margin + i_horiz * (axis_width + horiz_sep) \
                        - colorbar_offset - colorbar_width
                    y_bottom = bottom_margin \
                        + i_vert * (axis_height + vert_sep) \
                        + axis_height * (1.0 - colorbar_length_fraction) / 2.0
                    x_width = colorbar_width
                    y_height = axis_height * colorbar_length_fraction
                elif colorbar_location == "right":
                    x_left = left_margin + (i_horiz + 1) * axis_width \
                        + i_horiz * horiz_sep + colorbar_offset
                    y_bottom = bottom_margin \
                        + i_vert * (axis_height + vert_sep) \
                        + axis_height * (1.0 - colorbar_length_fraction) / 2.0
                    x_width = colorbar_width
                    y_height = axis_height * colorbar_length_fraction
                elif colorbar_location == "bottom":
                    x_left = left_margin + i_horiz * (axis_width + horiz_sep) \
                        + axis_width * (1.0 - colorbar_length_fraction) / 2.0
                    y_bottom = bottom_margin \
                        + i_vert * (axis_height + vert_sep) \
                        - colorbar_offset - colorbar_width
                    x_width = axis_width * colorbar_length_fraction
                    y_height = colorbar_width
                elif colorbar_location == "top":
                    x_left = left_margin + i_horiz * (axis_width + horiz_sep) \
                        + axis_width * (1.0 - colorbar_length_fraction) / 2.0
                    y_bottom = bottom_margin \
                        + (i_vert + 1) * axis_height + i_vert * vert_sep \
                        + colorbar_offset
                    x_width = axis_width * colorbar_length_fraction
                    y_height = colorbar_width
                else:
                    raise Exception("unknown colorbar: %s" % str(colorbar))
                if colorbar == "individual" \
                        or (colorbar == "shared"
                            and ((colorbar_location == "left" and i_horiz == 0)
                                 or (colorbar_location == "right" and i_horiz == n_horiz - 1)
                                 or (colorbar_location == "bottom" and i_vert == 0)
                                 or (colorbar_location == "top" and i_vert == n_vert - 1))):
                        cbar_axes = figure.add_axes([x_left / figure_width,
                                                     y_bottom / figure_height,
                                                     x_width / figure_width,
                                                     y_height / figure_height])
                        if colorbar_location == "left":
                            cbar_axes.xaxis.tick_left()
                            cbar_axes.xaxis.set_label_position('left')
                        if colorbar_location == "top":
                            cbar_axes.xaxis.tick_top()
                            cbar_axes.xaxis.set_label_position('top')
                        if colorbar == "individual":
                            cbar_axes_array[-1].append(cbar_axes)
                        elif colorbar == "shared":
                            cbar_axes_array.append(cbar_axes)
        return (figure, axes_array, cbar_axes_array)
    elif colorbar == "single":
        total_width = n_horiz * axis_width + (n_horiz - 1) * horiz_sep
        total_height = n_vert * axis_height + (n_vert - 1) * vert_sep
        if colorbar_location == "left":
            x_left = left_margin - colorbar_offset - colorbar_width
            y_bottom = bottom_margin \
                + total_height * (1.0 - colorbar_length_fraction) / 2.0
            x_width = colorbar_width
            y_height = total_height * colorbar_length_fraction
        elif colorbar_location == "right":
            x_left = left_margin + total_width + colorbar_offset
            y_bottom = bottom_margin \
                + total_height * (1.0 - colorbar_length_fraction) / 2.0
            x_width = colorbar_width
            y_height = total_height * colorbar_length_fraction
        elif colorbar_location == "bottom":
            x_left = left_margin \
                + total_width * (1.0 - colorbar_length_fraction) / 2.0
            y_bottom = bottom_margin - colorbar_offset - colorbar_width
            x_width = total_width * colorbar_length_fraction
            y_height = colorbar_width
        elif colorbar_location == "top":
            x_left = left_margin \
                + total_width * (1.0 - colorbar_length_fraction) / 2.0
            y_bottom = bottom_margin + total_height + colorbar_offset
            x_width = total_width * colorbar_length_fraction
            y_height = colorbar_width
        else:
            raise Exception("unknown colorbar: %s" % str(colorbar))
        cbar_axes = figure.add_axes([x_left / figure_width,
                                     y_bottom / figure_height,
                                     x_width / figure_width,
                                     y_height / figure_height])
        if colorbar_location == "left":
            cbar_axes.xaxis.tick_left()
            cbar_axes.xaxis.set_label_position('left')
        if colorbar_location == "top":
            cbar_axes.xaxis.tick_top()
            cbar_axes.xaxis.set_label_position('top')
        return (figure, axes_array, cbar_axes)
    elif colorbar is not None:
        raise Exception("unknown colorbar: %s" % str(colorbar))
    return (figure, axes_array)

(figure, axes_array, cbar_axes_array) \
    = make_fig_array(1,2,figure_width=6.9,
                     colorbar_offset=0.1,
                                top_margin=0.25, bottom_margin=0.45,
                                left_margin=1, right_margin=0.6,
                                vert_sep=0.3, horiz_sep=1.6,
                                colorbar="individual",colorbar_location="right",
                                share_y_axes=False)

######## first row ###############
filename = 'urban_plume/out/urban_plume2_process.nc'
ncf = scipy.io.netcdf_file(filename)
time_grid_edges = ncf.variables["time_grid_edges"].data
diversity_edges = ncf.variables["diversity_edges"].data
time_diversity_dist = ncf.variables["time_diversity_dist"].data

time = ncf.variables["time"].data / 3600
avg_part_entropy = numpy.exp(ncf.variables["avg_part_entropy"].data)
entropy_of_avg_part = numpy.exp(ncf.variables["entropy_of_avg_part"].data)
#tot_entropy_ratio = ncf.variables["tot_entropy_ratio"].data
tot_entropy_ratio = (avg_part_entropy-1) / (entropy_of_avg_part-1)
ncf.close()

axes = axes_array[0][0]
cbar_axes = cbar_axes_array[0][0]
d = numpy.ma.masked_less_equal(time_diversity_dist, 0)
vmin = 10**math.floor(math.log10(d.min()))
vmax = 10**math.ceil(math.log10(d.max()))
p = axes.imshow(numpy.flipud(d.transpose()), interpolation='nearest',
                extent=[time_grid_edges.min(), time_grid_edges.max(),
                        diversity_edges.min(), diversity_edges.max()],
                norm = matplotlib.colors.LogNorm(), aspect='auto')

axes.set_xscale("linear")
axes.set_xlabel(r"time $t$ / h")
axes.set_xlim([0,48])
axes.set_xticks([0, 12, 24, 36, 48])

axes.set_yscale("linear")
axes.set_ylabel(r"particle diversity $D_i$")
axes.set_ylim(0, 10)

axes.grid(True)

cbar = figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext(),
                       orientation='vertical')
cbar.solids.set_edgecolor("face")
cbar.solids.set_rasterized(True)
cbar_axes.xaxis.set_label_position('top')
cbar.set_label(r"$n(t, D_i)$ / $\rm m^{-3}$")
mpl_helper.add_boxed_text(axes, "(a)")

axes = axes_array[0][1]
axes.plot(time, avg_part_entropy, "b-")
axes.plot(time, entropy_of_avg_part, "k:")
axes.set_xlabel(r"time $t$ / h")

axes.set_ylabel(r"diversity $D_{\alpha}, D{\gamma}$")
axes.set_ylim([0,10])

axes2 =  axes.twinx()
axes2.plot(time, tot_entropy_ratio, "r--", markersize = 2)
axes2.set_ylabel(r"mixing state index $\chi$")
axes2.set_ylim([0.4,0.9])
axes2.set_xlim([0,48])
axes2.set_xticks([0, 12, 24, 36, 48])

axes.annotate(r"$D_{\gamma}$", (12,7.8),
              verticalalignment="bottom", horizontalalignment="right",
              bbox = dict(facecolor='white', edgecolor='white'),
              xytext=(0, 5), textcoords='offset points')
axes.annotate(r"$D_{\alpha}$", (36,3),
              verticalalignment="bottom", horizontalalignment="right",
              bbox = dict(facecolor='white', edgecolor='white'),
              xytext=(0, 5), textcoords='offset points')
axes.annotate(r"$\chi$", (36,7),
              verticalalignment="bottom", horizontalalignment="right",
              bbox = dict(facecolor='white', edgecolor='white'),
              xytext=(0, 5), textcoords='offset points')

axes.grid(True)
mpl_helper.add_boxed_text(axes, "(b)")
mpl_helper.remove_fig_array_axes(axes_array, remove_y_axes=False)

out_filename = "urban_plume_all_d.pdf"
figure.savefig(out_filename, dpi=1200)
print out_filename
