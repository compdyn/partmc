#!/usr/bin/env python
import numpy as np

def plot_bar(hist, edges, ax, coef = 1e9, **kwargs):
    min_interval =  np.min(edges[1:] - edges[:-1])
    ax.bar(edges[:-1], hist / (edges[1:] - edges[:-1]) / coef, align = "edge", width = min_interval * 0.85, **kwargs)

def bar_minmax(hist, edges, coef = 1e9):
    min_interval =  np.min(edges[1:] - edges[:-1])
    bar_data = hist / (edges[1:] - edges[:-1]) / coef
    bar_min, bar_max = bar_data.min(), bar_data.max()
    return bar_min, bar_max

def mod_xlabel_to_log(ax_hist):
    x_ticks = ax_hist.get_xticks()
    x_ticks = [int(xtick) for xtick in x_ticks]
    x_ticks_label = ["10$^{" + "%.1f" % (xtick) + "}$" for xtick in x_ticks]
    ax_hist.set_xticks(x_ticks, x_ticks_label)
    return ax_hist
    
def mod_ylabel_to_log(ax_hist):
    y_ticks = ax_hist.get_yticks()
    y_ticks = [int(ytick) for ytick in y_ticks]
    y_ticks_label = ["10$^{" + "%.1f" % (ytick) + "}$" for ytick in y_ticks]
    ax_hist.set_yticks(y_ticks, y_ticks_label)
    return ax_hist
