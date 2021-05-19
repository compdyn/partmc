#!/usr/bin/env python

import numpy as np

print('#' * 78)

print('\ndensities of species A, B, C, D')
rho = np.array([1000, 1500, 2000, 2500])
print(f'rho = {rho} kg m^{-3}')

print('\nmass fractions in populations P1 and P2')
w1 = np.array([0.2, 0.3, 0.1, 0.4])
w2 = np.array([0.5, 0.2, 0.2, 0.1])
print(f'w1 = {w1}')
print(f'w2 = {w2}')

print('\ndiameters of particles in populations P1 and P2')
D1 = 1e-7
D2 = 1e-6
print(f'D1 = {D1} m')
print(f'D2 = {D2} m')

print('\nnumber concentrations of populations P1 and P2')
N1 = 1e13
N2 = 1e10
print(f'N1 = {N1:e} m^{-3}')
print(f'N2 = {N2:e} m^{-3}')

print('\naverage densities in each population')
rho1 = 1/np.sum(w1/rho)
rho2 = 1/np.sum(w2/rho)
print(f'rho1 = {rho1} kg m^{-3}')
print(f'rho2 = {rho2} kg m^{-3}')

print('\nper-particle mass in each population')
m1 = np.pi/6 * D1**3 * rho1
m2 = np.pi/6 * D2**3 * rho2
print(f'm1 = {m1} kg')
print(f'm2 = {m2} kg')

print('\nper-species masses in each population')
mu1 = m1 * w1
mu2 = m2 * w2
print(f'mu1 = {mu1} kg')
print(f'mu2 = {mu2} kg')

def compute_chi(name, mu1, mu2, weighting):
    """Compute chi for the population.
    name is a string describing the population
    mu1, mu2 are the per-particle species masses in each population
    weighting is a string that can be one of:
        'mass'
    """

    print('\nper-species masses in each population')
    print(f'mu1 = {mu1} kg')
    print(f'mu2 = {mu2} kg')

    print('\ntotal per-particle mass in each population')
    m1 = np.sum(mu1)
    m2 = np.sum(mu2)
    print(f'm1 = {m1}')
    print(f'm2 = {m2}')

    print('\nmass concentrations of each population')
    M1 = m1 * N1
    M2 = m2 * N2
    print(f'M1 = {M1} kg m^{-3}')
    print(f'M2 = {M2} kg m^{-3}')

    print(f'\nmass fractions in each population')
    p1 = mu1 / np.sum(mu1)
    p2 = mu2 / np.sum(mu2)
    print(f'p1 = {p1}')
    print(f'p2 = {p2}')

    print('\nentropy of particles in each population')
    Halpha1 = -np.sum(p1 * np.log(p1))
    Halpha2 = -np.sum(p2 * np.log(p2))
    print(f'Halpha1 = {Halpha1}')
    print(f'Halpha2 = {Halpha2}')

    print('\ndiversity of particles in each population')
    Dalpha1 = np.exp(Halpha1)
    Dalpha2 = np.exp(Halpha2)
    print(f'Dalpha1 = {Dalpha1}')
    print(f'Dalpha2 = {Dalpha2}')

    print(f'\nweighting = {weighting}')
    if (weighting == 'mass'):
        pw1 = M1
        pw2 = M2
    else:
        raise Exception(f'unknown weighting: {weighting}')
    print(f'pw1 = {pw1}')
    print(f'pw2 = {pw2}')

    print('\naverage particle entropy over entire population')
    Halpha = (pw1*Halpha1 + pw2*Halpha2) / (pw1 + pw2)
    print(f'Halpha = {Halpha}')

    print('\naverage particle diversity over entire population')
    Dalpha = np.exp(Halpha)
    print(f'Dalpha = {Dalpha}')

    print('\naverage particle mass fractions')
    p = (pw1*p1 + pw2*p2) / (pw1 + pw2)
    print(f'p = {p}')

    print('\nentropy of average particle')
    Hgamma = -np.sum(p * np.log(p))
    print(f'Hgamma = {Hgamma}')

    print('\ndiversity of average particle')
    Dgamma = np.exp(Hgamma)
    print(f'Dgamma = {Dgamma}')

    print('\nmixing state index')
    chi = (Dalpha - 1) / (Dgamma - 1)
    print(f'chi = {chi}')

    output_name = f'ref_{name}.txt'
    print(f'\nwriting output to {output_name}')
    with open(output_name, 'w') as f:
        f.write(f'{Dalpha} {Dgamma} {chi}\n')

print('\n' + '#' * 78)
use_mu1 = mu1
use_mu2 = mu2
compute_chi('all_species', use_mu1, use_mu2, 'mass')

print('\n' + '#' * 78)
use_mu1 = np.array([mu1[0], mu1[2] + mu1[3]])
use_mu2 = np.array([mu2[0], mu2[2] + mu2[3]])
compute_chi('groups', use_mu1, use_mu2, 'mass')

print('\n' + '#' * 78)

