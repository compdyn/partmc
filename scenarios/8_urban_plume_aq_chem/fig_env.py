#!/usr/bin/env python

import scipy.io
import sys
sys.path.append("../../tool")
import numpy as np
import mpl_helper
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
import partmc
 
in_dir = "out/"
in_file_pattern = "urban_plume_no_aq_chem_w_0001_.*.nc"
env_state_history = partmc.read_history(partmc.env_state_t, in_dir, in_file_pattern)
time = np.array([env_state_history[i][0] for i in range(len(env_state_history))])
rh = [env_state_history[i][1].relative_humidity for i in range(len(env_state_history))]
temp = [env_state_history[i][1].temperature for i in range(len(env_state_history))]


(figure, axes) = mpl_helper.make_fig(figure_width=5,
                                     top_margin=0.5, bottom_margin=0.45,
                                     left_margin=0.65, right_margin=1.2)

print time/60.-time[0]/60
print rh
axes.set_xscale("linear")
axes.set_yscale("linear")
axes.plot(time/60.-time[0]/60., rh)

axes.set_ylabel(r"relative humidity")
axes.set_xlabel(r"time in minutes")
axes.set_xlim(0,50)
axes.set_xticks([0, 10, 20, 30, 40, 50])
axes.grid(True)
figure.savefig("figs/time_rh_w.pdf")
plt.close()

(figure, axes) = mpl_helper.make_fig(figure_width=5,
                                     top_margin=0.5, bottom_margin=0.45,
                                     left_margin=0.65, right_margin=1.2)

axes.set_xscale("linear")
axes.set_yscale("linear")
axes.plot(time/60.-time[0]/60., temp)

axes.set_ylabel(r"temperature in K")
axes.set_xlabel(r"time in minutes")
axes.set_ylim(274,295)
axes.set_xlim(0,50)
axes.set_xticks([0, 10, 20, 30, 40, 50])
axes.grid(True)
figure.savefig("figs/time_temp_w.pdf")
plt.close()

