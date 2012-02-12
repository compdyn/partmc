#!/usr/bin/env python

import math

H = 0.95
T = 287.5
Tdot = -5.0 / (10 * 60)
P0 = 611.2 * math.exp(7.45 * math.log(10) * (T - 273.15) / (T - 38))
dP0_dT = P0 * 7.45 * math.log(10) * (273.15 - 38) / (T - 38)**2
Hdot_env = - 1 / P0 * dP0_dT * Tdot * H
t_sat = (1 - H) / Hdot_env

print "H = ", H
print "T = ", T
print "Tdot = ", Tdot
print "P0 = ", P0
print "dP0_dT = ", dP0_dT
print "Hdot_env = ", Hdot_env
print "t_sat = ", t_sat
