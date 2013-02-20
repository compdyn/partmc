#!/usr/bin/env python
import sys
sys.path.append('../../tool/')
import subprocess
import os
import math
import numpy as np
import mpl_helper
import matplotlib.pyplot as plt

# plot number distribution
barrel_num_1 = np.loadtxt("ref_0212/ref_aero_size_num_regrid.txt")
barrel_num_2 = np.loadtxt("ref_0909/ref_aero_size_num_regrid.txt")
barrel_num_3 = np.loadtxt("ref_0925/ref_aero_size_num_regrid.txt")
(figure, axes) = mpl_helper.make_fig(colorbar=False)
b1 = axes.semilogx(barrel_num_1[:,0], barrel_num_1[:,1], color='r')
b2 = axes.semilogx(barrel_num_2[:,0], barrel_num_2[:,1], color='b')
b3 = axes.semilogx(barrel_num_3[:,0], barrel_num_3[:,1], color='g')
axes.set_xlabel("Diameter (m)")
axes.set_ylabel(r"Number concentration ($\mathrm{m}^{-3}$)")
axes.grid()
axes.set_ylim(0, 1.8e12)
#axes.legend([b1,b2,b3], ('02/12', '09/09', '09/25'),loc='upper right')
axes.annotate("(1) high dil.", color='r',
            xy=(2e-7,1.25e12), xycoords='data',
            xytext=(0.88,0.8), textcoords='axes fraction',
            size=12, va="center", ha="center",
            arrowprops=dict(arrowstyle="fancy",
                            fc="0.6", ec="r",
                                connectionstyle="angle3,angleA=0,angleB=-90"), 
            )
axes.annotate("(2) low dil.,\n high conc.",color='b',
            xy=(0.4e-7,1.35e12), xycoords='data',
            xytext=(0.12,0.8), textcoords='axes fraction',
            size=12, va="center", ha="center",
            arrowprops=dict(arrowstyle="fancy",
                            fc="0.6", ec="b",
                                connectionstyle="angle3,angleA=0,angleB=-90"),
            )
axes.annotate("(3) low dil.,\n low conc.",color='g',
            xy=(1.1e-7,0.25e12), xycoords='data',
            xytext=(0.88,0.4), textcoords='axes fraction',
            size=12, va="center", ha="center",
            arrowprops=dict(arrowstyle="fancy",
                            fc="0.6", ec="g",
                                connectionstyle="angle3,angleA=0,angleB=-90"),
            )
filename_out = "aero_num_size_initial.pdf"
figure.savefig(filename_out)

barrel_num_1 = np.loadtxt("ref_0322/ref_aero_size_num_regrid.txt")
barrel_num_2 = np.loadtxt("ref_0908/ref_aero_size_num_regrid.txt")
(figure, axes) = mpl_helper.make_fig(colorbar=False)
b1 = axes.semilogx(barrel_num_1[:,0], barrel_num_1[:,1], color='r')
b2 = axes.semilogx(barrel_num_2[:,0], barrel_num_2[:,1], color='b')
axes.set_xlabel("Diameter (m)")
axes.set_ylabel(r"Number concentration ($\mathrm{m}^{-3}$)")
axes.grid()
axes.set_ylim(0, 1.7e12)
axes.annotate("(4)", color='r',
            xy=(1.2e-7,1.2e12), xycoords='data',
            xytext=(0.7,0.8), textcoords='axes fraction',
            size=12, va="center", ha="center",
            arrowprops=dict(arrowstyle="fancy",
                            fc="0.6", ec="r",
                                connectionstyle="angle3,angleA=0,angleB=-90"),
            )
axes.annotate("(5)",color='b',
            xy=(0.5e-7,0.5e12), xycoords='data',
            xytext=(0.12,0.5), textcoords='axes fraction',
            size=12, va="center", ha="center",
            arrowprops=dict(arrowstyle="fancy",
                            fc="0.6", ec="b",
                                connectionstyle="angle3,angleA=0,angleB=-90"),
            )
filename_out = "aero_num_size_initial_1.pdf"
figure.savefig(filename_out)
