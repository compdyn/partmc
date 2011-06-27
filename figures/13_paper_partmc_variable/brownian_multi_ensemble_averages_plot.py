#!/usr/bin/env python

import os, sys
import scipy.io
import scipy.stats
import numpy as np

sys.path.append("../../tool")
import partmc
import mpl_helper
import matplotlib
import config

top_label_width = 30
top_label_height = 6
top_label_voffset = 12.2

c = config.c_value
n = config.i_ens_max
r = scipy.stats.t.ppf((1 + c) / 2, n - 1)
conf_factor = r / np.sqrt(n)


f1 = "data/avg_ensemble_norm_num_brownian_600s.txt"
f2 = "data/avg_ensemble_norm_mass_brownian_600s.txt" 
f3 = "data/std_ensemble_norm_num_brownian_600s.txt"
f4 = "data/std_ensemble_norm_mass_brownian_600s.txt" 

f5 = "data/avg_ensemble_norm_num_brownian_60s.txt"
f6 = "data/avg_ensemble_norm_mass_brownian_60s.txt" 
f7 = "data/std_ensemble_norm_num_brownian_60s.txt"
f8 = "data/std_ensemble_norm_mass_brownian_60s.txt" 

f9 = "data/avg_ensemble_norm_num_brownian_6s.txt"
f10 = "data/avg_ensemble_norm_mass_brownian_6s.txt" 
f11 = "data/std_ensemble_norm_num_brownian_6s.txt"
f12 = "data/std_ensemble_norm_mass_brownian_6s.txt" 

f13 = "data/avg_ensemble_norm_num_brownian_1200s.txt"
f14 = "data/avg_ensemble_norm_mass_brownian_1200s.txt" 
f15 = "data/std_ensemble_norm_num_brownian_1200s.txt"
f16 = "data/std_ensemble_norm_mass_brownian_1200s.txt" 

f17 = "data/avg_ensemble_norm_num_brownian_1800s.txt"
f18 = "data/avg_ensemble_norm_mass_brownian_1800s.txt"
f19 = "data/std_ensemble_norm_num_brownian_1800s.txt"
f20 = "data/std_ensemble_norm_mass_brownian_1800s.txt"

f21 = "data/avg_ensemble_norm_num_brownian_3600s.txt"
f22 = "data/avg_ensemble_norm_mass_brownian_3600s.txt"
f23 = "data/std_ensemble_norm_num_brownian_3600s.txt"
f24 = "data/std_ensemble_norm_mass_brownian_3600s.txt"

f25 = "data/avg_ensemble_norm_num_brownian_1000p.txt"
f26 = "data/avg_ensemble_norm_mass_brownian_1000p.txt"
f27 = "data/std_ensemble_norm_num_brownian_1000p.txt"
f28 = "data/std_ensemble_norm_mass_brownian_1000p.txt"

f29 = "data/avg_ensemble_norm_num_brownian_100p.txt"
f30 = "data/avg_ensemble_norm_mass_brownian_100p.txt"
f31 = "data/std_ensemble_norm_num_brownian_100p.txt"
f32 = "data/std_ensemble_norm_mass_brownian_100p.txt"

f33 = "data/avg_ensemble_norm_num_brownian_10p.txt"
f34 = "data/avg_ensemble_norm_mass_brownian_10p.txt"
f35 = "data/std_ensemble_norm_num_brownian_10p.txt"
f36 = "data/std_ensemble_norm_mass_brownian_10p.txt"

num_avg_overall_600s = np.loadtxt(f1) / 1e6
mass_avg_overall_600s = np.loadtxt(f2)*1e9
num_std_overall_600s = np.loadtxt(f3) / 1e6
mass_std_overall_600s = np.loadtxt(f4)*1e9

num_avg_overall_60s = np.loadtxt(f5) / 1e6
mass_avg_overall_60s = np.loadtxt(f6)*1e9
num_std_overall_60s = np.loadtxt(f7) / 1e6
mass_std_overall_60s = np.loadtxt(f8)*1e9

num_avg_overall_6s = np.loadtxt(f9) / 1e6
mass_avg_overall_6s = np.loadtxt(f10)*1e9
num_std_overall_6s = np.loadtxt(f11) / 1e6
mass_std_overall_6s = np.loadtxt(f12)*1e9

num_avg_overall_1200s = np.loadtxt(f13) / 1e6
mass_avg_overall_1200s = np.loadtxt(f14)*1e9
num_std_overall_1200s = np.loadtxt(f15) / 1e6
mass_std_overall_1200s = np.loadtxt(f16)*1e9

num_avg_overall_1800s = np.loadtxt(f17) / 1e6
mass_avg_overall_1800s = np.loadtxt(f18)*1e9
num_std_overall_1800s = np.loadtxt(f19) / 1e6
mass_std_overall_1800s = np.loadtxt(f20)*1e9

num_avg_overall_3600s = np.loadtxt(f21) / 1e6
mass_avg_overall_3600s = np.loadtxt(f22)*1e9
num_std_overall_3600s = np.loadtxt(f23) / 1e6
mass_std_overall_3600s = np.loadtxt(f24)*1e9

num_avg_overall_1000p = np.loadtxt(f25) / 1e6
mass_avg_overall_1000p = np.loadtxt(f26)*1e9
num_std_overall_1000p = np.loadtxt(f27) / 1e6
mass_std_overall_1000p = np.loadtxt(f28)*1e9

num_avg_overall_100p = np.loadtxt(f29) / 1e6
mass_avg_overall_100p = np.loadtxt(f30)*1e9
num_std_overall_100p = np.loadtxt(f31) / 1e6
mass_std_overall_100p = np.loadtxt(f32)*1e9

num_avg_overall_10p = np.loadtxt(f33) / 1e6
mass_avg_overall_10p = np.loadtxt(f34)*1e9
num_std_overall_10p = np.loadtxt(f35) / 1e6
mass_std_overall_10p = np.loadtxt(f36)*1e9

x_array = range(1,101)

(figure, axes_array) = mpl_helper.make_fig_array(2,2, figure_width=config.figure_width_double, 
						 top_margin=0.41,
                                                 left_margin=0.65, right_margin = 1.2, horiz_sep=0.6, vert_sep = 0.3, share_y_axes = False)

axes = axes_array[1][0]

#time_step_6s = axes.errorbar(x_array, num_avg_overall_6s, conf_factor*num_std_overall_6s, fmt = 'g-')
#time_step_60s = axes.errorbar(x_array, num_avg_overall_60s, conf_factor*num_std_overall_60s, fmt = 'b-')
#time_step_600s = axes.errorbar(x_array, num_avg_overall_600s, conf_factor*num_std_overall_600s, fmt = 'r-')
#time_step_1200s = axes.errorbar(x_array, num_avg_overall_1200s, conf_factor*num_std_overall_1200s, fmt = 'm-')

time_step_6s = axes.plot(x_array, num_avg_overall_6s,  'g-')
time_step_6s_high = axes.plot(x_array, num_avg_overall_6s+conf_factor*num_std_overall_6s,  'g:', linewidth=0.2)
time_step_6s_low = axes.plot(x_array, num_avg_overall_6s-conf_factor*num_std_overall_6s,  'g:', linewidth=0.2)
time_step_60s = axes.plot(x_array, num_avg_overall_60s, 'b-')
time_step_60s_high = axes.plot(x_array, num_avg_overall_60s+conf_factor*num_std_overall_60s,  'b:', linewidth=0.2)
time_step_60s_low = axes.plot(x_array, num_avg_overall_60s-conf_factor*num_std_overall_60s,  'b:', linewidth=0.2)
time_step_600s = axes.plot(x_array, num_avg_overall_600s, 'r-')
time_step_600s_high = axes.plot(x_array, num_avg_overall_600s+conf_factor*num_std_overall_600s,  'r:', linewidth=0.2)
time_step_600s_low = axes.plot(x_array, num_avg_overall_600s-conf_factor*num_std_overall_600s,  'r:', linewidth=0.2)
time_step_1200s = axes.plot(x_array, num_avg_overall_1200s, 'm-')
time_step_1200s_high = axes.plot(x_array, num_avg_overall_1200s+conf_factor*num_std_overall_1200s,  'm:', linewidth=0.2)
time_step_1200s_low = axes.plot(x_array, num_avg_overall_1200s-conf_factor*num_std_overall_1200s,  'm:', linewidth=0.2)
time_step_1800s = axes.plot(x_array, num_avg_overall_1800s, 'k-')
time_step_1800s_high = axes.plot(x_array, num_avg_overall_1800s+conf_factor*num_std_overall_1800s,  'k:', linewidth=0.2)
time_step_1800s_low = axes.plot(x_array, num_avg_overall_1800s-conf_factor*num_std_overall_1800s,  'k:', linewidth=0.2)
time_step_3600s = axes.plot(x_array, num_avg_overall_3600s, 'c-')
time_step_3600s_high = axes.plot(x_array, num_avg_overall_3600s+conf_factor*num_std_overall_3600s,  'c:', linewidth=0.2)
time_step_3600s_low = axes.plot(x_array, num_avg_overall_3600s-conf_factor*num_std_overall_3600s,  'c:', linewidth=0.2)


axes.set_xscale("log")
axes.set_yscale("log")

axes.set_xlim([0, 100])
axes.set_ylim([num_avg_overall_60s.min(), num_avg_overall_3600s.max()])
axes.grid(True)
axes.text(0.5, 1.15, r'number', horizontalalignment='center',
           verticalalignment='baseline', transform=axes.transAxes, )
pos = axes.transAxes.transform_point((0.5, 1)) * 72.0 / figure.get_dpi()
axes.add_artist(matplotlib.patches.FancyBboxPatch((pos[0] - top_label_width / 2.0,
                                                   pos[1] + top_label_voffset),
                                                  top_label_width, top_label_height,
                                                  boxstyle=matplotlib.patches.BoxStyle.Round(pad=3),
                                                  facecolor='white', edgecolor='black',
                                                  transform=None, clip_on=False))

#axes.set_xlabel(r"ensemble size $R$")
axes.set_ylabel(r" E$\| \langle  n_{\rm p}(11\, {\rm h}) - n_{\rm sect}(11 \, {\rm h}) \rangle \|_2$")
#axes.legend((time_step_6s, time_step_60s, time_step_600s, time_step_1200s, time_step_1800s, time_step_3600s), (r"6 s", r"60 s", r"600 s", r"1200 s", r"1800 s", r"3600"))

axes = axes_array[1][1]

time_step_6s = axes.plot(x_array, mass_avg_overall_6s,  'g-')
time_step_6s_high = axes.plot(x_array, mass_avg_overall_6s+conf_factor*mass_std_overall_6s,  'g:', linewidth=0.2)
time_step_6s_low = axes.plot(x_array, mass_avg_overall_6s-conf_factor*mass_std_overall_6s,  'g:', linewidth=0.2)
time_step_60s = axes.plot(x_array, mass_avg_overall_60s, 'b-')
time_step_60s_high = axes.plot(x_array, mass_avg_overall_60s+conf_factor*mass_std_overall_60s,  'b:', linewidth=0.2)
time_step_60s_low = axes.plot(x_array, mass_avg_overall_60s-conf_factor*mass_std_overall_60s,  'b:', linewidth=0.2)
time_step_600s = axes.plot(x_array, mass_avg_overall_600s, 'r-')
time_step_600s_high = axes.plot(x_array, mass_avg_overall_600s+conf_factor*mass_std_overall_600s,  'r:', linewidth=0.2)
time_step_600s_low = axes.plot(x_array, mass_avg_overall_600s-conf_factor*mass_std_overall_600s,  'r:', linewidth=0.2)
time_step_1200s = axes.plot(x_array, mass_avg_overall_1200s, 'm-')
time_step_1200s_high = axes.plot(x_array, mass_avg_overall_1200s+conf_factor*mass_std_overall_1200s,  'm:', linewidth=0.2)
time_step_1200s_low = axes.plot(x_array, mass_avg_overall_1200s-conf_factor*mass_std_overall_1200s,  'm:', linewidth=0.2)
time_step_1800s = axes.plot(x_array, mass_avg_overall_1800s, 'k-')
time_step_1800s_high = axes.plot(x_array, mass_avg_overall_1800s+conf_factor*mass_std_overall_1800s,  'k:', linewidth=0.2)
time_step_1800s_low = axes.plot(x_array, mass_avg_overall_1800s-conf_factor*mass_std_overall_1800s,  'k:', linewidth=0.2)
time_step_3600s = axes.plot(x_array, mass_avg_overall_3600s, 'c-')
time_step_3600s_high = axes.plot(x_array, mass_avg_overall_3600s+conf_factor*mass_std_overall_3600s,  'c:', linewidth=0.2)
time_step_3600s_low = axes.plot(x_array, mass_avg_overall_3600s-conf_factor*mass_std_overall_3600s,  'c:', linewidth=0.2)

axes.set_xscale("log")
axes.set_yscale("log")

axes.set_xlim([0, 100])
axes.set_ylim([mass_avg_overall_60s.min(), mass_avg_overall_60s.max()])
axes.grid(True)
axes.text(0.5, 1.15, r'mass', horizontalalignment='center',
           verticalalignment='baseline', transform=axes.transAxes)
pos = axes.transAxes.transform_point((0.5, 1)) * 72.0 / figure.get_dpi()
axes.add_artist(matplotlib.patches.FancyBboxPatch((pos[0] - top_label_width / 2.0,
                                                   pos[1] + top_label_voffset),
                                                  top_label_width, top_label_height,
                                                  boxstyle=matplotlib.patches.BoxStyle.Round(pad=3),
                                                  facecolor='white', edgecolor='black',
                                                  transform=None, clip_on=False))

#axes.set_xlabel(r"ensemble size $R$")
axes.set_ylabel(r" E$\| \langle  m_{\rm p}(11\, {\rm h}) - m_{\rm sect}(11 \, {\rm h}) \rangle \|_2$")
axes.legend((time_step_3600s, time_step_1800s, time_step_1200s, time_step_600s, time_step_60s, time_step_6s), (r"$\Delta t=3600$ s", r"$\Delta t=1800$ s", r"$\Delta t=1200$ s", r"$\Delta t=600$ s", r"$\Delta t=60$ s", r"$\Delta t=6$ s"), loc='center left', bbox_to_anchor = (1.0,0.5),ncol=1)


axes = axes_array[0][0]

time_step_60s = axes.plot(x_array, num_avg_overall_60s, 'b-')
time_step_60s_high = axes.plot(x_array, num_avg_overall_60s+conf_factor*num_std_overall_60s,  'b:', linewidth=0.2)
time_step_60s_low = axes.plot(x_array, num_avg_overall_60s-conf_factor*num_std_overall_60s,  'b:', linewidth=0.2)
time_step_1000p = axes.plot(x_array, num_avg_overall_1000p, 'r-')
time_step_1000p_high = axes.plot(x_array, num_avg_overall_1000p+conf_factor*num_std_overall_1000p,  'r:', linewidth=0.2)
time_step_1000p_low = axes.plot(x_array, num_avg_overall_1000p-conf_factor*num_std_overall_1000p,  'r:', linewidth=0.2)
time_step_100p = axes.plot(x_array, num_avg_overall_100p, 'g-')
time_step_100p_high = axes.plot(x_array, num_avg_overall_100p+conf_factor*num_std_overall_100p,  'g:', linewidth=0.2)
time_step_100p_low = axes.plot(x_array, num_avg_overall_100p-conf_factor*num_std_overall_100p,  'g:', linewidth=0.2)
time_step_10p = axes.plot(x_array, num_avg_overall_10p, 'c-')
time_step_10p_high = axes.plot(x_array, num_avg_overall_10p+conf_factor*num_std_overall_10p,  'c:', linewidth=0.2)
time_step_10p_low = axes.plot(x_array, num_avg_overall_10p-conf_factor*num_std_overall_10p,  'c:', linewidth=0.2)

axes.set_xscale("log")
axes.set_yscale("log")

axes.set_xlim([0, 100])
axes.set_ylim([num_avg_overall_60s.min(), num_avg_overall_10p.max()])
axes.grid(True)
axes.set_xlabel(r"ensemble size $R$")
axes.set_ylabel(r" E$\| \langle  n_{\rm p}(11\, {\rm h}) - n_{\rm sect}(11 \, {\rm h}) \rangle \|_2$")


axes = axes_array[0][1]

time_step_60s = axes.plot(x_array, mass_avg_overall_60s, 'b-')
time_step_60s_high = axes.plot(x_array, mass_avg_overall_60s+conf_factor*mass_std_overall_60s,  'b:', linewidth=0.2)
time_step_60s_low = axes.plot(x_array, mass_avg_overall_60s-conf_factor*mass_std_overall_60s,  'b:', linewidth=0.2)
time_step_1000p = axes.plot(x_array, mass_avg_overall_1000p, 'r-')
time_step_1000p_high = axes.plot(x_array, mass_avg_overall_1000p+conf_factor*mass_std_overall_1000p,  'r:', linewidth=0.2)
time_step_1000p_low = axes.plot(x_array, mass_avg_overall_1000p-conf_factor*mass_std_overall_1000p,  'r:', linewidth=0.2)
time_step_100p = axes.plot(x_array, mass_avg_overall_100p, 'g-')
time_step_100p_high = axes.plot(x_array, mass_avg_overall_100p+conf_factor*mass_std_overall_100p,  'g:', linewidth=0.2)
time_step_100p_low = axes.plot(x_array, mass_avg_overall_100p-conf_factor*mass_std_overall_100p,  'g:', linewidth=0.2)
time_step_10p = axes.plot(x_array, mass_avg_overall_10p, 'c-')
time_step_10p_high = axes.plot(x_array, mass_avg_overall_10p+conf_factor*mass_std_overall_10p,  'c:', linewidth=0.2)
time_step_10p_low = axes.plot(x_array, mass_avg_overall_10p-conf_factor*mass_std_overall_10p,  'c:', linewidth=0.2)

axes.set_xscale("log")
axes.set_yscale("log")

axes.set_xlim([0, 100])
axes.set_ylim([mass_avg_overall_60s.min(), mass_avg_overall_10p.max()])
axes.grid(True)
axes.set_xlabel(r"ensemble size $R$")
axes.set_ylabel(r" E$\| \langle  m_{\rm p}(11\, {\rm h}) - m_{\rm sect}(11 \, {\rm h}) \rangle \|_2$")
axes.legend((time_step_10p, time_step_100p, time_step_1000p, time_step_60s), (r"$N_{\rm p} = 10$ ", r"$N_{\rm p} = 100$", r"$N_{\rm p} = 10^3$", r"$N_{\rm p} = 10^4$"), loc='center left', bbox_to_anchor = (1.0,0.5),ncol=1)

mpl_helper.remove_fig_array_axes(axes_array, remove_y_axes=False)

figure.savefig("figs/brownian_multi_ensemble_averages.pdf")

