#!/usr/bin/env python
import sys
sys.path.append('../../tool/')
import math
import numpy
import os
import mpl_helper
import matplotlib.pyplot as plt

def air_mean_free_path(temp, press):
    boltzmann = 1.3806505e-23
    avogad = 6.02214179e23
    air_molec_weight = 2.89644e-2
    univ_gas_const = 8.314472

    boltz = boltzmann * 1e7
    mwair = air_molec_weight * 1e3
    rgas = univ_gas_const * 1e-2
    rhoair = 0.001 * ((press/1.01325e5) * mwair / (rgas * temp))

    viscosd = (1.8325e-4 * (296.16 + 120) / (temp + 120)) * (temp / 296.16)**1.5
    viscosk = viscosd / rhoair
    gasspeed = math.sqrt(8 * boltz * temp * avogad / (math.pi * mwair))
    free_path = 2 * viscosk / gasspeed * 1e-2
    return free_path

def slip_correct(diameter, temp, press):
    a_slip = 1.142
    q_slip = 0.588
    b_slip = 0.999
    free_path = air_mean_free_path(temp, press)

    slip_correct = 1 + a_slip * free_path / diameter \
                     + q_slip * free_path / diameter \
                     * math.exp(-b_slip * diameter / free_path)
    return slip_correct

def diff_coef(diameter, temp, press):

    boltzmann = 1.3806505e-23
    air_dyn_visc = 1.78e-5

    diff_coef = boltzmann * temp * slip_correct(diameter, temp, press) \
                     / (6 * math.pi * air_dyn_visc * diameter)
    return diff_coef

def diffus_BL_thick(diameter, temp, press, kD, a):
    delta_D = kD * diff_coef(diameter, temp, press)**a
    return delta_D

def wall_loss(diameter, temp, press, kD, a):
    A_D = 1.988
    V = 0.2093
    wall_loss_rate = diff_coef(diameter, temp, press) * A_D \
                     / diffus_BL_thick(diameter, temp, press, kD, a) / V
    return wall_loss_rate

# make particle diameter array (m)
Dp = numpy.logspace(-8, -6, 100)

temp = 293.4
press = 1e5
list_diff_coef = []
for diameter in Dp:
    list_diff_coef.append(diff_coef(diameter, temp, press))

(figure, axes) = mpl_helper.make_fig(colorbar=False)
axes.semilogx(Dp, list_diff_coef)
axes.set_xlabel("Diameter (m)")
axes.set_ylabel(r"Diffusional coefficient ($\mathrm{m}^{2}$ / s)")
axes.grid()
filename_out = "plot_diff_coef.pdf"
figure.savefig(filename_out)

if not os.path.exists("plots_wall_loss"):
   os.mkdir("plots_wall_loss")

i = 0
list_wall_loss = []
for prefactor in numpy.arange(0.005,0.055,0.005):
    for exponent in numpy.arange(0.1,0.31,0.01):
        i = i + 1
        for diameter in Dp:
            list_wall_loss.append(wall_loss(diameter, temp, press, prefactor, exponent))
        (figure, axes) = mpl_helper.make_fig(colorbar=False)
        axes.semilogx(Dp, list_wall_loss)
        axes.set_xlabel("Diameter (m)")
        axes.set_ylabel(r"Wall loss rate ($\mathrm{s}^{-1}$)")
        axes.set_title("prefactor = %.3f, exponent = %.3f" %(prefactor, exponent))
        axes.grid()
        axes.set_ylim(0,0.0062)
        filename_out = "plots_wall_loss/plot_wall_loss_%04d.png" %(i)
        figure.set_dpi(600)
        figure.savefig(filename_out)
        list_wall_loss = []
