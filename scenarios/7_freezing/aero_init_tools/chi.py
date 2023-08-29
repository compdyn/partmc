#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from pdb import set_trace

class aero_init(object):

    def __init__(self, mean = 1e-6, logstd = 0.326, Npart = 1000, species = ["OIN", "ILT", "KLN"]):
        self.Npart = Npart
        self.species = species
        self.Ns = len(species)

        self.aero_dens = {
            "H2O": 1000,
            "OIN": 2600,
            "ILT": 2750,
            "KLN": 2650,
        }
        self.species_dens = []
        for spec in self.species:
            self.species_dens.append(self.aero_dens[spec])
        self.species_dens = np.array(self.species_dens)
        self.species_dens = self.species_dens.reshape(self.Ns, 1)
        self.logDp = np.random.normal(np.log10(mean), logstd, size = (self.Npart))
        self.Dp = 10 **self.logDp
        self.Vp = 1.0/6 * np.pi * self.Dp **3


    def init_random(self):

        self.species_volume_ratio  = np.random.rand(self.Ns, self.Npart)
        self.species_volume_ratio_one_parm = np.sum(self.species_volume_ratio, axis = 0)
        self.species_volume_ratio = self.species_volume_ratio / self.species_volume_ratio_one_parm
        self.species_volume = self.species_volume_ratio * self.Vp
        self.update()

    def init_external(self, species_total_volume_ratio):
        self.species_volume_ratio = np.full((self.Ns, self.Npart), 0.0)
        species_total_volume_ratio = np.array(species_total_volume_ratio)
        species_total_volume_ratio = species_total_volume_ratio / np.sum(species_total_volume_ratio)
        aero_ind = np.arange(self.Npart)
        np.random.shuffle(aero_ind)
        species_volume_cumsum = np.cumsum(self.Vp[aero_ind])
        species_volume_cumratio = species_volume_cumsum / species_volume_cumsum[-1]
        threshold_ratio = np.cumsum(species_total_volume_ratio)
        threshold = [0]
        for ratio in threshold_ratio:
            threshold.append(np.argmin(np.abs(species_volume_cumratio - ratio)) + 1)
        for i_spec in range(self.Ns):
            spec_index = aero_ind[threshold[i_spec]:threshold[i_spec+1]]
            self.species_volume_ratio[i_spec, spec_index] = 1.0
        self.species_volume = self.species_volume_ratio * self.Vp

        #species_total_volume = (self.species_volume_ratio * self.Vp).sum(axis=1)
        #species_total_volume_ratio = species_total_volume / np.sum(species_total_volume)
        self.update()

    def update(self):
        self.species_volume_ratio = self.species_volume / self.Vp
        self.species_mass = self.species_volume * self.species_dens
        self.aero_mass = np.sum(self.species_mass, axis = 0)
        self.species_total_mass = np.sum(self.species_mass, axis = 1)
        self.total_mass = np.sum(self.aero_mass)
        self.P_ai = self.species_mass / self.aero_mass
        self.P_i = self.aero_mass / self.total_mass
        self.P_a = self.species_total_mass / self.total_mass
        self.compute_chi()

    def compute_chi(self):
        self.Hi = -np.sum(np.where(self.P_ai == 0.0, 0.0, self.P_ai * np.log(self.P_ai)), axis = 0)
        self.Ha = np.sum(self.P_i * self.Hi)
        self.Hy = -np.sum(np.where(self.P_a == 0.0, 0.0, self.P_a * np.log(self.P_a)))
        self.Di = np.exp(self.Hi)
        self.Da = np.exp(self.Ha)
        self.Dy = np.exp(self.Hy)
        self.Db = self.Dy / self.Da
        self.chi = (self.Da - 1) / (self.Dy - 1)

    def mixing(self, Npars, chi_obj = None):
        assert Npars * 2 <= self.Npart
        aero_index = np.arange(self.Npart)
        np.random.shuffle(aero_index)
        parA_ind = aero_index[:Npars]
        parB_ind = aero_index[-Npars:]
        parA_vol = self.Vp[parA_ind]
        parB_vol = self.Vp[parB_ind]
        min_vol = np.min(np.array([parA_vol, parB_vol]), axis = 0)
        exchange_ratio = np.random.rand(len(min_vol)) 
        if not(chi_obj is None):
            exchange_ratio *= np.abs(self.chi - chi_obj) / np.abs(chi_obj)
        exchange_vol = exchange_ratio * min_vol
        exchange_species_vol_A = self.species_volume_ratio[:, parA_ind] * exchange_vol
        exchange_species_vol_B = self.species_volume_ratio[:, parB_ind] * exchange_vol
        self.species_volume[:, parB_ind] = self.species_volume[:, parB_ind] - exchange_species_vol_B + exchange_species_vol_A
        self.species_volume[:, parA_ind] = self.species_volume[:, parA_ind] - exchange_species_vol_A + exchange_species_vol_B
        self.update()
        #set_trace()
        #print(np.sum(self.species_volume, axis = 0))






if __name__ == "__main__":
    chi_list = []
    a = aero_init()
#a.init_random()
#a.compute_chi()
#plt.hist(a.logDp, bins = 100)
#plt.show()
    a.init_external([0.3, 0.5, 0.2])
#a.compute_chi()
    chi_obj = 1
    while True:
        a.mixing(Npars = 100)#, chi_obj = chi_obj)
#a.compute_chi()
        print("chi = ", a.chi)
        chi_list.append(a.chi)
        if np.abs(a.chi - chi_obj) < 0.0001:
            break
    print(a.species_volume_ratio)
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(chi_list)
    ax.set_xlabel("Iteration")
    ax.set_ylabel("χ")
    ax.grid()
    plt.show()


