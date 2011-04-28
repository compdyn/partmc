#!/usr/bin/env python

import os, sys
import numpy

sys.path.append("../../tool")
import partmc
import mpl_helper
import matplotlib

if not os.path.isdir("figs"):
    os.mkdir("figs")

def plot_kernel(prefix):
    kernel_data = numpy.loadtxt("data/" + prefix + "_kernel.txt")
    counts_data = numpy.loadtxt("data/" + prefix + "_count.txt")

    diam_edges = numpy.concatenate((counts_data[:,2], counts_data[-1:,3])) * 1e6
    diam_centers = counts_data[:,1] * 1e6
    n_parts = counts_data[:,4]
    min_bin = numpy.flatnonzero(n_parts).min()
    max_bin = numpy.flatnonzero(n_parts).max()
    diam_edges_nonzero = diam_edges[min_bin:(max_bin + 2)]
    diam_centers_nonzero = diam_centers[min_bin:(max_bin + 1)]
    n_parts_nonzero = n_parts[min_bin:(max_bin + 1)]
    n_bin = counts_data.shape[0]

    k_min = kernel_data[:,2].reshape((n_bin, n_bin))
    k_max = kernel_data[:,3].reshape((n_bin, n_bin))
    k_min = numpy.ma.masked_less_equal(k_min, 0)
    k_max = numpy.ma.masked_less_equal(k_max, 0)
    k_ratio = k_max / k_min
    k_min_nonzero = k_min[min_bin:(max_bin + 1), min_bin:(max_bin + 1)]
    k_max_nonzero = k_max[min_bin:(max_bin + 1), min_bin:(max_bin + 1)]
    k_ratio_nonzero = k_ratio[min_bin:(max_bin + 1), min_bin:(max_bin + 1)]

    n_prod = numpy.zeros_like(k_max)
    for i in range(n_bin):
        for j in range(n_bin):
            if i == j:
                n_prod[i,j] = 0.5 * float(n_parts[i]) * (float(n_parts[j]) - 1)
            else:
                n_prod[i,j] = 0.5 * float(n_parts[i]) * float(n_parts[j])
    n_prod_nonzero = n_prod[min_bin:(max_bin + 1), min_bin:(max_bin + 1)]

    (figure, axes) = mpl_helper.make_fig()
    axes.loglog(diam_centers, n_parts)
    axes.set_xlim(min(diam_edges), max(diam_edges))
    axes.set_xlabel(r"diameter / $\rm \mu m$")
    axes.set_ylabel(r"number of particles")
    figure.savefig("figs/" + prefix + "_count.pdf")
    
    (figure, axes, cbar_axes) = mpl_helper.make_fig(axis_ratio=1, colorbar=True, right_margin=1.2)
    p = axes.pcolor(diam_edges, diam_edges, n_prod, norm=matplotlib.colors.LogNorm(), linewidths=0.1)
    axes.set_xscale("log")
    axes.set_yscale("log")
    axes.set_xlim(min(diam_edges), max(diam_edges))
    axes.set_ylim(min(diam_edges), max(diam_edges))
    axes.set_xlabel(r"diameter $D_1$ / $\rm \mu m$")
    axes.set_ylabel(r"diameter $D_2$ / $\rm \mu m$")
    cbar = figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext(), orientation='vertical')
    cbar_axes.xaxis.set_label_position('top')
    cbar.set_label(r"number of reaction pairs $n(D_1) \big( n(D_2) - \delta(D_1,D_2) \big)$")
    figure.savefig("figs/" + prefix + "_reactions.pdf")
    
    (figure, axes, cbar_axes) = mpl_helper.make_fig(axis_ratio=1, colorbar=True, right_margin=1.2)
    p = axes.pcolor(diam_edges, diam_edges, n_prod * k_max, norm=matplotlib.colors.LogNorm(), linewidths=0.1)
    axes.set_xscale("log")
    axes.set_yscale("log")
    axes.set_xlim(min(diam_edges), max(diam_edges))
    axes.set_ylim(min(diam_edges), max(diam_edges))
    axes.set_xlabel(r"diameter $D_1$ / $\rm \mu m$")
    axes.set_ylabel(r"diameter $D_2$ / $\rm \mu m$")
    cbar = figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext(), orientation='vertical')
    cbar_axes.xaxis.set_label_position('top')
    cbar.set_label(r"number of events $K(D_1, D_2) n(D_1) \big( n(D_2) - \delta(D_1,D_2) \big)$")
    figure.savefig("figs/" + prefix + "_events.pdf")
    
    (figure, axes, cbar_axes) = mpl_helper.make_fig(axis_ratio=1, colorbar=True, right_margin=1.2)
    p = axes.pcolor(diam_edges, diam_edges, k_max, norm=matplotlib.colors.LogNorm(), linewidths=0.1)
    axes.set_xscale("log")
    axes.set_yscale("log")
    axes.set_xlim(min(diam_edges), max(diam_edges))
    axes.set_ylim(min(diam_edges), max(diam_edges))
    axes.set_xlabel(r"diameter $D_1$ / $\rm \mu m$")
    axes.set_ylabel(r"diameter $D_2$ / $\rm \mu m$")
    cbar = figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext(), orientation='vertical')
    cbar_axes.xaxis.set_label_position('top')
    cbar.set_label(r"kernel $K(D_1,D_2)$ / $(\rm m^3\ s^{-1})$")
    figure.savefig("figs/" + prefix + "_kernel_max.pdf")

    (figure, axes, cbar_axes) = mpl_helper.make_fig(axis_ratio=1, colorbar=True, right_margin=1.2)
    p = axes.pcolor(diam_edges, diam_edges, k_ratio, linewidths=0.1)
    axes.set_xscale("log")
    axes.set_yscale("log")
    axes.set_xlim(min(diam_edges), max(diam_edges))
    axes.set_ylim(min(diam_edges), max(diam_edges))
    axes.set_xlabel(r"diameter $D_1$ / $\rm \mu m$")
    axes.set_ylabel(r"diameter $D_2$ / $\rm \mu m$")
    cbar = figure.colorbar(p, cax=cbar_axes, orientation='vertical')
    cbar_axes.xaxis.set_label_position('top')
    cbar.set_label(r"kernel ratio $K_{\rm max}(D_1,D_2) / K_{\rm min}(D_1,D_2)$")
    figure.savefig("figs/" + prefix + "_kernel_ratio.pdf")

    (figure, axes) = mpl_helper.make_fig()
    axes.loglog(diam_centers_nonzero, n_parts_nonzero)
    axes.set_xlim(min(diam_edges_nonzero), max(diam_edges_nonzero))
    axes.set_xlabel(r"diameter / $\rm \mu m$")
    axes.set_ylabel(r"number of particles")
    figure.savefig("figs/" + prefix + "_count_nonzero.pdf")
    
    (figure, axes, cbar_axes) = mpl_helper.make_fig(axis_ratio=1, colorbar=True, right_margin=1.2)
    p = axes.pcolor(diam_edges_nonzero, diam_edges_nonzero, n_prod_nonzero, norm=matplotlib.colors.LogNorm(), linewidths=0.1)
    axes.set_xscale("log")
    axes.set_yscale("log")
    axes.set_xlim(min(diam_edges_nonzero), max(diam_edges_nonzero))
    axes.set_ylim(min(diam_edges_nonzero), max(diam_edges_nonzero))
    axes.set_xlabel(r"diameter $D_1$ / $\rm \mu m$")
    axes.set_ylabel(r"diameter $D_2$ / $\rm \mu m$")
    cbar = figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext(), orientation='vertical')
    cbar_axes.xaxis.set_label_position('top')
    cbar.set_label(r"number of reaction pairs $n(D_1) \big( n(D_2) - \delta(D_1,D_2) \big)$")
    figure.savefig("figs/" + prefix + "_reactions_nonzero.pdf")
    
    (figure, axes, cbar_axes) = mpl_helper.make_fig(axis_ratio=1, colorbar=True, right_margin=1.2)
    p = axes.pcolor(diam_edges_nonzero, diam_edges_nonzero, n_prod_nonzero * k_max_nonzero, norm=matplotlib.colors.LogNorm(), linewidths=0.1)
    axes.set_xscale("log")
    axes.set_yscale("log")
    axes.set_xlim(min(diam_edges_nonzero), max(diam_edges_nonzero))
    axes.set_ylim(min(diam_edges_nonzero), max(diam_edges_nonzero))
    axes.set_xlabel(r"diameter $D_1$ / $\rm \mu m$")
    axes.set_ylabel(r"diameter $D_2$ / $\rm \mu m$")
    cbar = figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext(), orientation='vertical')
    cbar_axes.xaxis.set_label_position('top')
    cbar.set_label(r"number of events $K(D_1, D_2) n(D_1) \big( n(D_2) - \delta(D_1,D_2) \big)$")
    figure.savefig("figs/" + prefix + "_events_nonzero.pdf")
    
    (figure, axes, cbar_axes) = mpl_helper.make_fig(axis_ratio=1, colorbar=True, right_margin=1.2)
    p = axes.pcolor(diam_edges_nonzero, diam_edges_nonzero, k_max_nonzero, norm=matplotlib.colors.LogNorm(), linewidths=0.1)
    axes.set_xscale("log")
    axes.set_yscale("log")
    axes.set_xlim(min(diam_edges_nonzero), max(diam_edges_nonzero))
    axes.set_ylim(min(diam_edges_nonzero), max(diam_edges_nonzero))
    axes.set_xlabel(r"diameter $D_1$ / $\rm \mu m$")
    axes.set_ylabel(r"diameter $D_2$ / $\rm \mu m$")
    cbar = figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext(), orientation='vertical')
    cbar_axes.xaxis.set_label_position('top')
    cbar.set_label(r"kernel $K(D_1,D_2)$ / $(\rm m^3\ s^{-1})$")
    figure.savefig("figs/" + prefix + "_kernel_max_nonzero.pdf")
    
    (figure, axes, cbar_axes) = mpl_helper.make_fig(axis_ratio=1, colorbar=True, right_margin=1.2)
    p = axes.pcolor(diam_edges_nonzero, diam_edges_nonzero, k_ratio_nonzero, linewidths=0.1)
    axes.set_xscale("log")
    axes.set_yscale("log")
    axes.set_xlim(min(diam_edges_nonzero), max(diam_edges_nonzero))
    axes.set_ylim(min(diam_edges_nonzero), max(diam_edges_nonzero))
    axes.set_xlabel(r"diameter $D_1$ / $\rm \mu m$")
    axes.set_ylabel(r"diameter $D_2$ / $\rm \mu m$")
    cbar = figure.colorbar(p, cax=cbar_axes, orientation='vertical')
    cbar_axes.xaxis.set_label_position('top')
    cbar.set_label(r"kernel ratio $K_{\rm max}(D_1,D_2) / K_{\rm min}(D_1,D_2)$")
    figure.savefig("figs/" + prefix + "_kernel_ratio_nonzero.pdf")

plot_kernel("brownian_w0")
plot_kernel("brownian_w-3")
plot_kernel("sedi_w0")
plot_kernel("sedi_w-1")
