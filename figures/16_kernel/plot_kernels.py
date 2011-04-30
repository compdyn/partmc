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

    n_bin = counts_data.shape[0]

    k_min = kernel_data[:,2].reshape((n_bin, n_bin))
    k_max = kernel_data[:,3].reshape((n_bin, n_bin))
    n_possible = kernel_data[:,4].reshape((n_bin, n_bin))
    r_samp = kernel_data[:,5].reshape((n_bin, n_bin))
    n_samp_mean = kernel_data[:,6].reshape((n_bin, n_bin))
    r_samp_rate = kernel_data[:,7].reshape((n_bin, n_bin))
    n_samp_mean_rate = kernel_data[:,8].reshape((n_bin, n_bin))

    diam_centers = counts_data[:,1] * 1e6
    diam_edges = numpy.concatenate((counts_data[:,2], counts_data[-1:,3])) * 1e6
    n_parts = counts_data[:,4]

    min_bin = numpy.flatnonzero(n_parts).min()
    max_bin = numpy.flatnonzero(n_parts).max()
    diam_edges_nonzero = diam_edges[min_bin:(max_bin + 2)]
    diam_centers_nonzero = diam_centers[min_bin:(max_bin + 1)]
    n_parts_nonzero = n_parts[min_bin:(max_bin + 1)]

    k_min = numpy.ma.masked_less_equal(k_min, 0)
    k_max = numpy.ma.masked_less_equal(k_max, 0)
    n_possible = numpy.ma.masked_less_equal(n_possible, 0)
    r_samp = numpy.ma.masked_less_equal(r_samp, 0)
    n_samp_mean = numpy.ma.masked_less_equal(n_samp_mean, 0)
    r_samp_rate = numpy.ma.masked_less_equal(r_samp_rate, 0)
    n_samp_mean_rate = numpy.ma.masked_less_equal(n_samp_mean_rate, 0)

    k_ratio = k_max / k_min

    k_min_nonzero = k_min[min_bin:(max_bin + 1), min_bin:(max_bin + 1)]
    k_max_nonzero = k_max[min_bin:(max_bin + 1), min_bin:(max_bin + 1)]
    k_ratio_nonzero = k_ratio[min_bin:(max_bin + 1), min_bin:(max_bin + 1)]
    n_possible_nonzero = n_possible[min_bin:(max_bin + 1), min_bin:(max_bin + 1)]
    r_samp_nonzero = r_samp[min_bin:(max_bin + 1), min_bin:(max_bin + 1)]
    n_samp_mean_nonzero = n_samp_mean[min_bin:(max_bin + 1), min_bin:(max_bin + 1)]
    r_samp_rate_nonzero = r_samp_rate[min_bin:(max_bin + 1), min_bin:(max_bin + 1)]
    n_samp_mean_rate_nonzero = n_samp_mean_rate[min_bin:(max_bin + 1), min_bin:(max_bin + 1)]

    n_parts_x = []
    n_parts_y = []
    for i in range(len(n_parts)):
        n_parts_x.append(diam_edges[i])
        n_parts_y.append(n_parts[i])
        n_parts_x.append(diam_edges[i + 1])
        n_parts_y.append(n_parts[i])

    (figure, axes) = mpl_helper.make_fig()
    axes.loglog(n_parts_x, n_parts_y)
    axes.set_xlim(min(diam_edges), max(diam_edges))
    axes.set_ylim(bottom=1e-1)
    axes.set_xlabel(r"diameter / $\rm \mu m$")
    axes.set_ylabel(r"number of particles")
    figure.savefig("figs/" + prefix + "_count.pdf")
    
    (figure, axes, cbar_axes) = mpl_helper.make_fig(axis_ratio=1, colorbar=True, right_margin=1.2)
    p = axes.pcolor(diam_edges, diam_edges, n_possible, norm=matplotlib.colors.LogNorm(), linewidths=0.1)
    axes.set_xscale("log")
    axes.set_yscale("log")
    axes.set_xlim(min(diam_edges), max(diam_edges))
    axes.set_ylim(min(diam_edges), max(diam_edges))
    axes.set_xlabel(r"diameter $D_1$ / $\rm \mu m$")
    axes.set_ylabel(r"diameter $D_2$ / $\rm \mu m$")
    cbar = figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext(), orientation='vertical')
    cbar_axes.xaxis.set_label_position('top')
    cbar.set_label(r"number of reaction pairs $n_i(n_j - \delta_{ij})$")
    figure.savefig("figs/" + prefix + "_n_possible.pdf")
    
    (figure, axes, cbar_axes) = mpl_helper.make_fig(axis_ratio=1, colorbar=True, right_margin=1.2)
    p = axes.pcolor(diam_edges, diam_edges, r_samp, norm=matplotlib.colors.LogNorm(), linewidths=0.1)
    axes.set_xscale("log")
    axes.set_yscale("log")
    axes.set_xlim(min(diam_edges), max(diam_edges))
    axes.set_ylim(min(diam_edges), max(diam_edges))
    axes.set_xlabel(r"diameter $D_1$ / $\rm \mu m$")
    axes.set_ylabel(r"diameter $D_2$ / $\rm \mu m$")
    cbar = figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext(), orientation='vertical')
    cbar_axes.xaxis.set_label_position('top')
    cbar.set_label(r"number of samples per pair $K_{ij} \Delta t / V_{\rm comp}$")
    figure.savefig("figs/" + prefix + "_r_samp.pdf")
    
    (figure, axes, cbar_axes) = mpl_helper.make_fig(axis_ratio=1, colorbar=True, right_margin=1.2)
    p = axes.pcolor(diam_edges, diam_edges, r_samp_rate, norm=matplotlib.colors.LogNorm(), linewidths=0.1)
    axes.set_xscale("log")
    axes.set_yscale("log")
    axes.set_xlim(min(diam_edges), max(diam_edges))
    axes.set_ylim(min(diam_edges), max(diam_edges))
    axes.set_xlabel(r"diameter $D_1$ / $\rm \mu m$")
    axes.set_ylabel(r"diameter $D_2$ / $\rm \mu m$")
    cbar = figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext(), orientation='vertical')
    cbar_axes.xaxis.set_label_position('top')
    cbar.set_label(r"number of samples per pair rate $K_{ij} / V_{\rm comp}$")
    figure.savefig("figs/" + prefix + "_r_samp_rate.pdf")
    
    (figure, axes, cbar_axes) = mpl_helper.make_fig(axis_ratio=1, colorbar=True, right_margin=1.2)
    p = axes.pcolor(diam_edges, diam_edges, n_samp_mean, norm=matplotlib.colors.LogNorm(), linewidths=0.1)
    axes.set_xscale("log")
    axes.set_yscale("log")
    axes.set_xlim(min(diam_edges), max(diam_edges))
    axes.set_ylim(min(diam_edges), max(diam_edges))
    axes.set_xlabel(r"diameter $D_1$ / $\rm \mu m$")
    axes.set_ylabel(r"diameter $D_2$ / $\rm \mu m$")
    cbar = figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext(), orientation='vertical')
    cbar_axes.xaxis.set_label_position('top')
    cbar.set_label(r"number of events $(K_{ij} \Delta t / V) n_i(n_j - \delta_{ij})$")
    figure.savefig("figs/" + prefix + "_n_samp_mean.pdf")
    
    (figure, axes, cbar_axes) = mpl_helper.make_fig(axis_ratio=1, colorbar=True, right_margin=1.2)
    p = axes.pcolor(diam_edges, diam_edges, n_samp_mean_rate, norm=matplotlib.colors.LogNorm(), linewidths=0.1)
    axes.set_xscale("log")
    axes.set_yscale("log")
    axes.set_xlim(min(diam_edges), max(diam_edges))
    axes.set_ylim(min(diam_edges), max(diam_edges))
    axes.set_xlabel(r"diameter $D_1$ / $\rm \mu m$")
    axes.set_ylabel(r"diameter $D_2$ / $\rm \mu m$")
    cbar = figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext(), orientation='vertical')
    cbar_axes.xaxis.set_label_position('top')
    cbar.set_label(r"number of events rate $(K_{ij} / V) n_i(n_j - \delta_{ij})$")
    figure.savefig("figs/" + prefix + "_n_samp_mean_rate.pdf")
    
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
    cbar.set_label(r"kernel $K_{ij}^{\rm max}$ / $(\rm m^3\ s^{-1})$")
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
    cbar.set_label(r"kernel ratio $K_{ij}^{\rm max} / K_{ij}^{\rm min}$")
    figure.savefig("figs/" + prefix + "_kernel_ratio.pdf")

    (figure, axes) = mpl_helper.make_fig()
    axes.loglog(n_parts_x, n_parts_y)
    axes.set_xlim(min(diam_edges_nonzero), max(diam_edges_nonzero))
    axes.set_ylim(bottom=1e-1)
    axes.set_xlabel(r"diameter / $\rm \mu m$")
    axes.set_ylabel(r"number of particles")
    figure.savefig("figs/" + prefix + "_count_nonzero.pdf")
    
    (figure, axes, cbar_axes) = mpl_helper.make_fig(axis_ratio=1, colorbar=True, right_margin=1.2)
    p = axes.pcolor(diam_edges_nonzero, diam_edges_nonzero, n_possible_nonzero, norm=matplotlib.colors.LogNorm(), linewidths=0.1)
    axes.set_xscale("log")
    axes.set_yscale("log")
    axes.set_xlim(min(diam_edges_nonzero), max(diam_edges_nonzero))
    axes.set_ylim(min(diam_edges_nonzero), max(diam_edges_nonzero))
    axes.set_xlabel(r"diameter $D_1$ / $\rm \mu m$")
    axes.set_ylabel(r"diameter $D_2$ / $\rm \mu m$")
    cbar = figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext(), orientation='vertical')
    cbar_axes.xaxis.set_label_position('top')
    cbar.set_label(r"number of reaction pairs $n_i(n_j - \delta_{ij})$")
    figure.savefig("figs/" + prefix + "_n_possible_nonzero.pdf")
    
    (figure, axes, cbar_axes) = mpl_helper.make_fig(axis_ratio=1, colorbar=True, right_margin=1.2)
    p = axes.pcolor(diam_edges_nonzero, diam_edges_nonzero, r_samp_nonzero, norm=matplotlib.colors.LogNorm(), linewidths=0.1)
    axes.set_xscale("log")
    axes.set_yscale("log")
    axes.set_xlim(min(diam_edges_nonzero), max(diam_edges_nonzero))
    axes.set_ylim(min(diam_edges_nonzero), max(diam_edges_nonzero))
    axes.set_xlabel(r"diameter $D_1$ / $\rm \mu m$")
    axes.set_ylabel(r"diameter $D_2$ / $\rm \mu m$")
    cbar = figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext(), orientation='vertical')
    cbar_axes.xaxis.set_label_position('top')
    cbar.set_label(r"number of samples per pair $K_{ij} \Delta t / V_{\rm comp}$")
    figure.savefig("figs/" + prefix + "_r_samp_nonzero.pdf")
    
    (figure, axes, cbar_axes) = mpl_helper.make_fig(axis_ratio=1, colorbar=True, right_margin=1.2)
    p = axes.pcolor(diam_edges_nonzero, diam_edges_nonzero, r_samp_rate_nonzero, norm=matplotlib.colors.LogNorm(), linewidths=0.1)
    axes.set_xscale("log")
    axes.set_yscale("log")
    axes.set_xlim(min(diam_edges_nonzero), max(diam_edges_nonzero))
    axes.set_ylim(min(diam_edges_nonzero), max(diam_edges_nonzero))
    axes.set_xlabel(r"diameter $D_1$ / $\rm \mu m$")
    axes.set_ylabel(r"diameter $D_2$ / $\rm \mu m$")
    cbar = figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext(), orientation='vertical')
    cbar_axes.xaxis.set_label_position('top')
    cbar.set_label(r"number of samples per pair rate $K_{ij} / V_{\rm comp}$")
    figure.savefig("figs/" + prefix + "_r_samp_rate_nonzero.pdf")
    
    (figure, axes, cbar_axes) = mpl_helper.make_fig(axis_ratio=1, colorbar=True, right_margin=1.2)
    p = axes.pcolor(diam_edges_nonzero, diam_edges_nonzero, n_samp_mean_nonzero, norm=matplotlib.colors.LogNorm(), linewidths=0.1)
    axes.set_xscale("log")
    axes.set_yscale("log")
    axes.set_xlim(min(diam_edges_nonzero), max(diam_edges_nonzero))
    axes.set_ylim(min(diam_edges_nonzero), max(diam_edges_nonzero))
    axes.set_xlabel(r"diameter $D_1$ / $\rm \mu m$")
    axes.set_ylabel(r"diameter $D_2$ / $\rm \mu m$")
    cbar = figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext(), orientation='vertical')
    cbar_axes.xaxis.set_label_position('top')
    cbar.set_label(r"number of events $(K_{ij} \Delta t / V) n_i(n_j - \delta_{ij})$")
    figure.savefig("figs/" + prefix + "_n_samp_mean_nonzero.pdf")
    
    (figure, axes, cbar_axes) = mpl_helper.make_fig(axis_ratio=1, colorbar=True, right_margin=1.2)
    p = axes.pcolor(diam_edges_nonzero, diam_edges_nonzero, n_samp_mean_rate_nonzero, norm=matplotlib.colors.LogNorm(), linewidths=0.1)
    axes.set_xscale("log")
    axes.set_yscale("log")
    axes.set_xlim(min(diam_edges_nonzero), max(diam_edges_nonzero))
    axes.set_ylim(min(diam_edges_nonzero), max(diam_edges_nonzero))
    axes.set_xlabel(r"diameter $D_1$ / $\rm \mu m$")
    axes.set_ylabel(r"diameter $D_2$ / $\rm \mu m$")
    cbar = figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext(), orientation='vertical')
    cbar_axes.xaxis.set_label_position('top')
    cbar.set_label(r"number of events rate $(K_{ij} / V) n_i(n_j - \delta_{ij})$")
    figure.savefig("figs/" + prefix + "_n_samp_mean_rate_nonzero.pdf")
    
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
    cbar.set_label(r"kernel $K_{ij}^{\rm max}$ / $(\rm m^3\ s^{-1})$")
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
    cbar.set_label(r"kernel ratio $K_{ij}^{\rm max} / K_{ij}^{\rm min}$")
    figure.savefig("figs/" + prefix + "_kernel_ratio_nonzero.pdf")

plot_kernel("brownian_w0")
plot_kernel("brownian_w-3")
plot_kernel("sedi_w0")
plot_kernel("sedi_w-1")
plot_kernel("sedi_w-3")
